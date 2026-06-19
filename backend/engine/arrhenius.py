"""
Arrheniusパラメータ逆解析エンジン（断熱n次反応・n=1固定）

設計書§5.1-5.4準拠。

モデル（§5.1）：
  状態 y = [T_degC, α]  （α ∈ [0,1]: 反応進行度）
  dα/dt = A · (1-α) · exp(-Ea / (R·T_K))       n=1固定
  dT/dt = (ΔH / (φ·Cp)) · dα/dt

最適化変数 θ = [Ea, log₁₀(A), ΔH, Cp]（§5.2）
CVP出力：A = 10^log_A（物理量1/s。log_Aは内部変数・顧客成果物に残さない）
NRMSE：range正規化 sqrt(mean((T_meas−T_sim)²)) / (T_max−T_min)（§5.2）
       HWSギャップリセット点はNRMSE算定から除外（強制正解点による過小評価防止）
HWSギャップ：dt>30分区間はスキップ・T実測値リセット・α引き継ぎ（§5.3）
CI：LM再フィット→conf_interval → bootstrap(n=200)フォールバック（§5.4）
"""
from __future__ import annotations

import math
import time
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import minimize

from .schema import CalorInput
from .solver import rk4_substep

R: float = 8.314            # 気体定数 [J/(mol·K)]
DT_GAP_THRESH_MIN: float = 30.0   # HWSギャップ判定閾値 [分]（超過区間は断熱自己発熱でない）
N_SUBSTEP: int = 40         # 各測定間隔内のRK4サブステップ数（§5.3）


@dataclass
class ArrheniusTheta:
    """Arrheniusモデルパラメータ（内部最適化用）"""
    Ea: float      # 活性化エネルギー [J/mol]
    log_A: float   # log₁₀(A)  ← 最適化変数（ODE内でA=10^log_Aに変換）
    dH: float      # 反応エンタルピー（総発熱量） [J/kg]
    Cp: float      # 比熱容量 [J/(kg·K)]
    phi: float     # 熱容量補正係数（≥1）


@dataclass
class ArrheniusFitResult:
    """逆解析結果（CVP params.json に格納するフィールドを含む）"""
    Ea: float              # 活性化エネルギー [J/mol]
    A: float               # 前指数因子 [1/s]（= 10^log_A・物理量）
    dH: float              # 反応エンタルピー [J/kg]
    Cp: float              # 比熱容量 [J/(kg·K)]
    n: int = 1             # 反応次数（v1固定）
    ci_lower_95: dict = field(default_factory=dict)   # 95%CI下限 {Ea, A, dH, Cp}（A空間）
    ci_upper_95: dict = field(default_factory=dict)   # 95%CI上限
    nrmse: float = 0.0
    method: str = ""
    n_evals: int = 0
    converged: bool = False
    elapsed_sec: float = 0.0


# ── 境界制約（§5.2 物理的妥当範囲） ─────────────────────────────────────────
_BOUNDS: list[tuple[float, float]] = [
    (50_000.0,    200_000.0),   # Ea [J/mol]（LIB熱暴走文献値）
    (6.0,         20.0),         # log₁₀(A)（反応速度定数対数スケール）
    (100_000.0,  2_000_000.0),  # ΔH [J/kg]（LIB化学系別発熱量範囲）
    (800.0,       1_500.0),     # Cp [J/(kg·K)]（LIB比熱容量文献値）
]

# 最適化初期値（境界の中央付近。合成テストではtheta_trueから離れた値を使用）
_X0: np.ndarray = np.array([120_000.0, 14.0, 500_000.0, 1_000.0])


def _ode_rhs(theta: ArrheniusTheta):
    """
    Arrhenius ODEの右辺 f(t, y) を返す（クロージャ）。

    y = [T_degC, alpha]
    f(t, y) = [dT/dt (°C/s), dα/dt (1/s)]
    """
    Ea = theta.Ea
    log_A = theta.log_A
    dH_phi_Cp = theta.dH / (theta.phi * theta.Cp)

    def f(t: float, y: np.ndarray) -> np.ndarray:
        T_K = y[0] + 273.15   # 摂氏→ケルビン変換必須（§5.1）
        # RK4サブステップ内での一時的オーバーシュートに対するガード
        if not math.isfinite(T_K) or T_K < 200.0 or T_K > 6000.0:
            return np.zeros(2)
        # αをODE内部でもクリップ（サブステップ間オーバーシュートで(1-α)^1が負になるのを防ぐ）
        alpha_safe = min(max(float(y[1]), 0.0), 1.0)
        if alpha_safe >= 1.0:
            return np.zeros(2)  # 反応完了：発熱ゼロ
        A = 10.0 ** log_A    # log→線形変換（§5.2）
        rate = A * (1.0 - alpha_safe) * math.exp(-Ea / (R * T_K))   # n=1固定
        if not math.isfinite(rate):
            return np.zeros(2)
        dT = dH_phi_Cp * rate
        return np.array([dT, rate])

    return f


def simulate_arc(
    t_arr: np.ndarray,
    T_meas_arr: np.ndarray,
    T_onset_degC: float,
    phi: float,
    theta: ArrheniusTheta,
) -> tuple[np.ndarray, np.ndarray]:
    """
    断熱ARCモデルで T_sim(t) を順方向計算する（§5.3）。

    HWSギャップ（dt > DT_GAP_THRESH_MIN）は積分スキップ・T実測値リセット・α引き継ぎ。
    ギャップリセット点は gap_mask に記録し、NRMSE算定から除外する
    （強制的に T_meas と一致させた点を正解として計上しないため）。

    Args:
        t_arr:        測定時刻配列 [min]
        T_meas_arr:   実測温度配列 [°C]（HWSギャップリセット用）
        T_onset_degC: 発熱開始温度 [°C]（初期条件 T(t_0)=T_meas[0]≈T_onset, α(t_0)=0）
        phi:          熱容量補正係数（≥1）
        theta:        Arrheniusパラメータ

    Returns:
        (T_sim, gap_mask):
          T_sim    - 模擬温度配列 [°C]（len(t_arr) と同一長）
          gap_mask - HWSギャップリセット点マスク（True=NRMSE除外対象）
    """
    f = _ode_rhs(theta)
    n = len(t_arr)
    T_sim = np.empty(n)
    T_sim[0] = T_onset_degC
    gap_mask = np.zeros(n, dtype=bool)

    y = np.array([T_onset_degC, 0.0])   # 初期状態 [T_degC, alpha=0]

    for i in range(n - 1):
        dt_min = t_arr[i + 1] - t_arr[i]

        if dt_min > DT_GAP_THRESH_MIN:
            # HWSギャップ：ARC装置強制加熱区間。積分スキップ・T実測値リセット・αは引き継ぎ
            y[0] = T_meas_arr[i + 1]
            gap_mask[i + 1] = True   # このリセット点はNRMSEから除外
        else:
            dt_sec = dt_min * 60.0
            try:
                y, _ = rk4_substep(f, t_arr[i] * 60.0, y, dt_sec, n_substep=N_SUBSTEP)
            except (OverflowError, ValueError, FloatingPointError):
                # 数値発散：残りを高温値で埋めて打ち切り
                T_sim[i + 1 :] = T_meas_arr.max() * 2.0
                return T_sim, gap_mask
            y[1] = min(max(y[1], 0.0), 1.0)   # αを[0,1]にクリップ（§5.3）

        T_sim[i + 1] = y[0]

    return T_sim, gap_mask


def compute_nrmse(
    T_meas: np.ndarray,
    T_sim: np.ndarray,
    gap_mask: Optional[np.ndarray] = None,
) -> float:
    """
    Range正規化RMSE（§5.2）。

    NRMSE = sqrt(mean((T_meas − T_sim)²)) / (T_max − T_min)

    分母は T_meas の断熱昇温幅（range）。平均値による正規化は使用しない
    （温度は原点任意の量ゆえ mean 正規化は誤差を過小評価する）。

    gap_mask が指定された場合はHWSギャップリセット点を除外する
    （T_meas[i]=T_sim[i] が保証されるリセット点を算定に含めると精度が過小評価される）。
    """
    if gap_mask is not None and gap_mask.any():
        keep = ~gap_mask
        T_meas = T_meas[keep]
        T_sim = T_sim[keep]
    if len(T_meas) == 0:
        return float("inf")
    denom = float(T_meas.max() - T_meas.min())
    if denom < 1e-6:
        return float("inf")
    rmse = float(np.sqrt(np.mean((T_meas - T_sim) ** 2)))
    return rmse / denom


def _make_objective(
    t_arr: np.ndarray,
    T_meas_arr: np.ndarray,
    T_onset_degC: float,
    phi: float,
    eval_counter: list[int],
    T_reset_arr: Optional[np.ndarray] = None,
):
    """
    NRMSE目的関数を生成する（scipy.optimize.minimize 用スカラー関数）。

    T_reset_arr: HWSギャップリセット用温度配列（Noneの場合はT_meas_arrを使用）。
                 Bootstrapでは原実測値をリセット基準として保持するため分離する。
    """
    if T_reset_arr is None:
        T_reset_arr = T_meas_arr

    def objective(x: np.ndarray) -> float:
        Ea, log_A, dH, Cp = float(x[0]), float(x[1]), float(x[2]), float(x[3])
        theta = ArrheniusTheta(Ea=Ea, log_A=log_A, dH=dH, Cp=Cp, phi=phi)
        T_sim, gap_mask = simulate_arc(t_arr, T_reset_arr, T_onset_degC, phi, theta)
        eval_counter[0] += 1
        nrmse = compute_nrmse(T_meas_arr, T_sim, gap_mask)
        return nrmse if math.isfinite(nrmse) else 1e6

    return objective


def fit_arrhenius(calor_input: CalorInput) -> ArrheniusFitResult:
    """
    calor-input v1 JSON から Arrhenius パラメータを逆解析する。

    3段階最適化（§5.2）：
      Stage 1: L-BFGS-B（大局的勾配降下）
      Stage 2: Nelder-Mead（勾配不要・局所精密化）
      Stage 3: デフォルト値フォールバック

    信頼区間（§5.4）：
      LM再フィット → conf_interval → bootstrap(n=200) フォールバック
    """
    t_start = time.monotonic()

    cond = calor_input.arc_conditions
    T_onset = cond.T_onset_degC
    phi = cond.phi

    # ── データ準備: T_onset 以降の発熱相データのみ使用 ────────────────────────
    t_all = np.array([p.time_min for p in calor_input.arc_data])
    T_all = np.array([p.T_degC for p in calor_input.arc_data])

    # np.searchsorted は単調ソートを前提とするため使用禁止
    # 実ARCデータはHWS待機相の揺らぎで非単調になりうる（設計書§5.3 onset_idx注記）
    candidates = np.nonzero(T_all >= T_onset)[0]
    onset_idx = int(candidates[0]) if len(candidates) > 0 else 0
    if onset_idx >= len(t_all) - 1:
        onset_idx = 0

    t_arr = t_all[onset_idx:]
    T_meas_arr = T_all[onset_idx:]
    T_onset_actual = float(T_meas_arr[0])   # 実データの第1点（≈T_onset）

    eval_counter: list[int] = [0]
    objective = _make_objective(t_arr, T_meas_arr, T_onset_actual, phi, eval_counter)

    best_x: np.ndarray = _X0.copy()
    best_cost: float = float("inf")
    best_method: str = "fallback"
    converged: bool = False

    # ── Stage 1: L-BFGS-B ──────────────────────────────────────────────────
    try:
        res1 = minimize(
            objective,
            _X0.copy(),
            method="L-BFGS-B",
            bounds=_BOUNDS,
            options={"maxiter": 500, "ftol": 1e-10, "gtol": 1e-8},
        )
        if math.isfinite(res1.fun) and res1.fun < best_cost:
            best_cost = float(res1.fun)
            best_x = res1.x.copy()
            best_method = "L-BFGS-B"
            converged = bool(res1.success)
    except Exception:
        pass

    # ── Stage 2: Nelder-Mead（L-BFGS-B 結果を初期値として精密化） ───────────
    try:
        res2 = minimize(
            objective,
            best_x,
            method="Nelder-Mead",
            options={"maxiter": 3000, "xatol": 1e-6, "fatol": 1e-9, "adaptive": True},
        )
        if math.isfinite(res2.fun) and res2.fun < best_cost:
            # 境界制約の事後確認（Nelder-Mead は bounds を無視する）
            x2 = res2.x
            in_bounds = all(
                lo <= float(xi) <= hi
                for xi, (lo, hi) in zip(x2, _BOUNDS)
            )
            if in_bounds:
                best_cost = float(res2.fun)
                best_x = x2.copy()
                best_method = "Nelder-Mead"
                converged = bool(res2.success)
    except Exception:
        pass

    Ea_opt, log_A_opt, dH_opt, Cp_opt = (
        float(best_x[0]), float(best_x[1]),
        float(best_x[2]), float(best_x[3]),
    )
    A_opt = 10.0 ** log_A_opt

    theta_opt = ArrheniusTheta(
        Ea=Ea_opt, log_A=log_A_opt, dH=dH_opt, Cp=Cp_opt, phi=phi
    )
    T_sim_opt, gap_mask_opt = simulate_arc(
        t_arr, T_meas_arr, T_onset_actual, phi, theta_opt
    )
    nrmse_opt = compute_nrmse(T_meas_arr, T_sim_opt, gap_mask_opt)

    # ── Stage 3+: LM再フィット → 信頼区間（§5.4） ─────────────────────────
    ci_lower, ci_upper = _compute_ci(
        t_arr, T_meas_arr, T_onset_actual, phi, best_x
    )

    return ArrheniusFitResult(
        Ea=Ea_opt,
        A=A_opt,
        dH=dH_opt,
        Cp=Cp_opt,
        n=1,
        ci_lower_95=ci_lower,
        ci_upper_95=ci_upper,
        nrmse=nrmse_opt,
        method=best_method,
        n_evals=eval_counter[0],
        converged=converged,
        elapsed_sec=time.monotonic() - t_start,
    )


# ── 信頼区間計算 ────────────────────────────────────────────────────────────

def _compute_ci(
    t_arr: np.ndarray,
    T_meas_arr: np.ndarray,
    T_onset_degC: float,
    phi: float,
    x_opt: np.ndarray,
) -> tuple[dict, dict]:
    """LM再フィット → conf_interval、失敗時 bootstrap(n=200) フォールバック。"""
    try:
        return _ci_via_lmfit(t_arr, T_meas_arr, T_onset_degC, phi, x_opt)
    except Exception as exc:
        warnings.warn(
            f"LM信頼区間計算失敗 ({exc!r})。bootstrap(n=200)にフォールバック。"
            " CI幅が過小評価される可能性あり（§5.4注記）。",
            RuntimeWarning,
            stacklevel=3,
        )
        return _ci_via_bootstrap(t_arr, T_meas_arr, T_onset_degC, phi, x_opt, n=200)


def _ci_via_lmfit(
    t_arr: np.ndarray,
    T_meas_arr: np.ndarray,
    T_onset_degC: float,
    phi: float,
    x_opt: np.ndarray,
) -> tuple[dict, dict]:
    """
    lmfit Levenberg-Marquardt + conf_interval で 95% CI を計算する（§5.4）。

    conf_interval() は leastsq 法後にのみ適用可能（L-BFGS-B/NM後には不可）。
    HWSギャップリセット点は残差から除外する。
    """
    from lmfit import Parameters, Minimizer, conf_interval as lm_ci

    def residuals(params) -> np.ndarray:
        theta = ArrheniusTheta(
            Ea=float(params["Ea"].value),
            log_A=float(params["log_A"].value),
            dH=float(params["dH"].value),
            Cp=float(params["Cp"].value),
            phi=phi,
        )
        T_sim, gap_mask = simulate_arc(t_arr, T_meas_arr, T_onset_degC, phi, theta)
        res_arr = T_meas_arr - T_sim
        return res_arr[~gap_mask]   # ギャップリセット点を残差から除外

    lm_params = Parameters()
    Ea0, log_A0, dH0, Cp0 = x_opt
    lm_params.add("Ea",    value=Ea0,    min=50_000.0,    max=200_000.0)
    lm_params.add("log_A", value=log_A0, min=6.0,         max=20.0)
    lm_params.add("dH",    value=dH0,    min=100_000.0,   max=2_000_000.0)
    lm_params.add("Cp",    value=Cp0,    min=800.0,       max=1_500.0)

    minner = Minimizer(residuals, lm_params)
    lm_result = minner.minimize(method="leastsq")

    # sigmas=[1, 2] → conf_interval returns entries sorted by sigma
    # entries[0][1] = −2σ value, entries[-1][1] = +2σ value（~95% CI）
    ci = lm_ci(minner, lm_result, sigmas=[1, 2])

    def lo(name: str) -> float:
        return float(ci[name][0][1])

    def hi(name: str) -> float:
        return float(ci[name][-1][1])

    # A空間への逆変換（log_A CI → A CI・§5.4）
    ci_lower = {"Ea": lo("Ea"), "A": 10.0 ** lo("log_A"), "dH": lo("dH"), "Cp": lo("Cp")}
    ci_upper = {"Ea": hi("Ea"), "A": 10.0 ** hi("log_A"), "dH": hi("dH"), "Cp": hi("Cp")}
    return ci_lower, ci_upper


def _ci_via_bootstrap(
    t_arr: np.ndarray,
    T_meas_arr: np.ndarray,
    T_onset_degC: float,
    phi: float,
    x_opt: np.ndarray,
    n: int = 200,
) -> tuple[dict, dict]:
    """
    Bootstrap法（n=200サンプル）で 95% CI を計算する（LM失敗時フォールバック・§5.4）。

    残差リサンプリング法：T_boot[非ギャップ点] = T_sim(θ*) + resample(T_meas − T_sim(θ*))
    HWSギャップリセット点は T_meas_arr の実測値を保持する。

    注意：各サンプルをL-BFGS-Bのみ・θ*始点で再フィットするため、
    多峰性景観では CI が過小評価される可能性がある（既知の制約）。
    """
    theta_opt = ArrheniusTheta(
        Ea=float(x_opt[0]), log_A=float(x_opt[1]),
        dH=float(x_opt[2]), Cp=float(x_opt[3]), phi=phi,
    )
    T_sim_opt, gap_mask_opt = simulate_arc(
        t_arr, T_meas_arr, T_onset_degC, phi, theta_opt
    )
    # ギャップリセット点を除いた残差（リセット点は常に残差=0ゆえリサンプリング不要）
    non_gap_residuals = (T_meas_arr - T_sim_opt)[~gap_mask_opt]

    rng = np.random.default_rng(seed=42)
    boot_params: list[np.ndarray] = []

    for _ in range(n):
        idx = rng.integers(0, len(non_gap_residuals), size=len(non_gap_residuals))
        # T_boot: 非ギャップ点のみリサンプル、ギャップ点は実測値のまま
        T_boot = T_meas_arr.copy()
        T_boot[~gap_mask_opt] = T_sim_opt[~gap_mask_opt] + non_gap_residuals[idx]

        ec_b: list[int] = [0]
        obj_b = _make_objective(
            t_arr, T_boot, T_onset_degC, phi, ec_b,
            T_reset_arr=T_meas_arr,   # ギャップリセットは原実測値で固定
        )
        try:
            res_b = minimize(
                obj_b,
                x_opt,
                method="L-BFGS-B",
                bounds=_BOUNDS,
                options={"maxiter": 200},
            )
            if math.isfinite(res_b.fun):
                boot_params.append(res_b.x.copy())
        except Exception:
            pass

    if len(boot_params) < 10:
        # サンプル不足：点推定値をそのまま返す
        pt = {"Ea": float(x_opt[0]), "A": 10.0 ** float(x_opt[1]),
              "dH": float(x_opt[2]), "Cp": float(x_opt[3])}
        return pt, pt

    boot_arr = np.array(boot_params)
    lo = np.percentile(boot_arr, 2.5, axis=0)
    hi = np.percentile(boot_arr, 97.5, axis=0)

    ci_lower = {"Ea": float(lo[0]), "A": 10.0 ** float(lo[1]),
                "dH": float(lo[2]), "Cp": float(lo[3])}
    ci_upper = {"Ea": float(hi[0]), "A": 10.0 ** float(hi[1]),
                "dH": float(hi[2]), "Cp": float(hi[3])}
    return ci_lower, ci_upper
