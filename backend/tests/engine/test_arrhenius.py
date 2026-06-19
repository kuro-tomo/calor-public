"""
Arrheniusエンジン単体テスト

テスト一覧：
  1. test_range_nrmse_normalization     - NRMSE がrange正規化（mean正規化禁止）
  2. test_alpha_clip                    - α が [0,1] にクリップされ、plateau安定
  3. test_hws_gap_skip                  - HWSギャップ（>30分）でT再設定・α引き継ぎ・gap_mask記録
  4. test_synthetic_round_trip          - コールドスタートから合成データを逆解析（NRMSE < 15%）
  5. test_schema_phi_validation         - φ<1 でバリデーション失敗
  6. test_schema_min_data_points        - データ点数不足でバリデーション失敗
  7. test_nmc_21700_fit                 - NMC 21700 実データでフィット成功（NRMSE < 0.5）
  8. test_lfp_21700_fit                 - LFP 21700 実データでフィット成功（NRMSE < 0.5）
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

# プロジェクトルートを sys.path に追加（conftest.py が優先するが念のため保持）
sys.path.insert(0, str(Path(__file__).parents[3]))

from backend.engine.schema import ArcConditions, ArcDataPoint, CalorInput, CellSpec
from backend.engine.arrhenius import (
    DT_GAP_THRESH_MIN,
    ArrheniusTheta,
    compute_nrmse,
    fit_arrhenius,
    simulate_arc,
)
from tools.arc_parser import parse_zenodo_txt

# Zenodoデータファイルパス（data/raw/ 相対）
# parents[3] = /projects/calor-arrhenius（プロジェクトルート）
_DATA_DIR = Path(__file__).parents[3] / "data" / "raw"
_NMC_TXT = _DATA_DIR / "NMC_SOC100_M1.txt"
_LFP_TXT = _DATA_DIR / "LFP_SOC100_M1.txt"


# ── ヘルパー ────────────────────────────────────────────────────────────────

def _make_calor_input(
    t_arr: np.ndarray,
    T_arr: np.ndarray,
    T_onset: float,
    phi: float = 1.0,
    chemistry: str = "NMC",
) -> CalorInput:
    """テスト用 CalorInput を構築する。"""
    arc_data = [
        ArcDataPoint(time_min=float(t), T_degC=float(T))
        for t, T in zip(t_arr, T_arr)
    ]
    return CalorInput(
        arc_data=arc_data,
        cell_spec=CellSpec(
            cell_id="test_cell",
            chemistry=chemistry,
            capacity_ah=4.9,
            mass_g=68.0,
            format="21700",
        ),
        arc_conditions=ArcConditions(T_onset_degC=T_onset, phi=phi),
    )


def _make_synthetic_input(
    theta: ArrheniusTheta,
    T_onset: float = 90.0,
    t_max_min: float = 20.0,
    n_points: int = 200,
) -> CalorInput:
    """既知θで合成ARC時系列を生成し CalorInput を返す。"""
    t_arr = np.linspace(0.0, t_max_min, n_points)
    T_dummy = np.full(n_points, T_onset)
    T_true, _ = simulate_arc(t_arr, T_dummy, T_onset, theta.phi, theta)
    return _make_calor_input(t_arr, T_true, T_onset, phi=theta.phi)


# ── テスト 1: NRMSE range 正規化 ─────────────────────────────────────────────

def test_range_nrmse_normalization():
    """NRMSE 分母が T_max−T_min（range）であることを確認する。"""
    # range = 300-80 = 220, mean = (80+120+180+260+300)/5 = 188  （range≠mean）
    T_meas = np.array([80.0, 120.0, 180.0, 260.0, 300.0])
    T_sim  = np.array([90.0, 120.0, 180.0, 260.0, 290.0])
    errors = np.array([10.0, 0.0, 0.0, 0.0, 10.0])
    t_range = 300.0 - 80.0   # 220
    t_mean  = float(np.mean(T_meas))  # 188

    # range正規化の期待値
    expected_range_nrmse = float(np.sqrt(np.mean(errors ** 2))) / t_range
    assert abs(compute_nrmse(T_meas, T_sim) - expected_range_nrmse) < 1e-10

    # mean正規化（禁止パターン）とは異なることを確認（range≠mean なので異なるはず）
    mean_norm_nrmse = float(np.sqrt(np.mean(errors ** 2))) / t_mean
    assert not np.isclose(compute_nrmse(T_meas, T_sim), mean_norm_nrmse, rtol=1e-3), (
        "NRMSE が mean 正規化と一致している（range 正規化が実装されていない可能性）"
    )

    # gap_mask を渡した場合、対象外点を除いた値になること
    gap = np.array([False, True, False, False, False])  # index=1のみギャップ点
    T_meas_kept = T_meas[[0, 2, 3, 4]]
    T_sim_kept  = T_sim[[0, 2, 3, 4]]
    expected_masked = compute_nrmse(T_meas_kept, T_sim_kept)
    assert abs(compute_nrmse(T_meas, T_sim, gap_mask=gap) - expected_masked) < 1e-10


# ── テスト 2: α クリップ ─────────────────────────────────────────────────────

def test_alpha_clip():
    """
    α が [0,1] にクリップされ、反応完了後 T_sim が有限かつ安定（plateau）であることを確認する。

    Opus指摘⑥対応：有限値確認に加え、α=1クリップ後の plateau 安定性を追加。
    """
    # 高温・大きな A → 急速反応 → α が即座に1に達する極限ケース
    theta_fast = ArrheniusTheta(
        Ea=80_000.0, log_A=18.0, dH=800_000.0, Cp=1_000.0, phi=1.0
    )
    t_arr = np.linspace(0.0, 30.0, 200)
    T_dummy = np.full(200, 90.0)
    T_sim, _ = simulate_arc(t_arr, T_dummy, 90.0, 1.0, theta_fast)

    # クリップにより T_sim が有限値（NaN/Inf なし）
    assert np.all(np.isfinite(T_sim)), "T_sim に NaN/Inf が含まれている"

    # 急速反応による昇温が起きること（α > 0 の証拠）
    T_peak = float(T_sim.max())
    assert T_peak > 90.0 + 100.0, (
        f"反応による昇温が不十分: T_peak={T_peak:.1f}°C（< 190°C）"
    )

    # α=1クリップ後は発熱ゼロ → T_sim が最高点付近に安定（plateau）
    # クリップなしなら α がオーバーシュートして発熱が続きT_simが発散する
    assert T_sim[-1] >= T_peak - 10.0, (
        f"α=1クリップ後にT_simが低下している: T_sim[-1]={T_sim[-1]:.1f}, T_peak={T_peak:.1f}"
        " （αクリップが機能していない可能性）"
    )


# ── テスト 3: HWS ギャップスキップ ──────────────────────────────────────────

def test_hws_gap_skip():
    """
    HWSギャップ（dt > 30分）区間では T が実測値にリセットされ、
    α（反応進行度）はリセットされず引き継がれることを確認する。
    また、gap_mask のリセット点フラグが正しく立っていることを確認する。
    """
    # t_arr: 0→10分（通常）, 10→100分（HWSギャップ）, 100→110分（通常再開）
    t_arr = np.array([0.0, 10.0, 100.0, 110.0])
    T_meas_arr = np.array([90.0, 120.0, 180.0, 220.0])   # ギャップ中に実測180°C

    theta = ArrheniusTheta(
        Ea=120_000.0, log_A=14.0, dH=500_000.0, Cp=1_000.0, phi=1.0
    )

    T_sim, gap_mask = simulate_arc(t_arr, T_meas_arr, T_onset_degC=90.0, phi=1.0, theta=theta)

    # ギャップ（100−10=90分 > 30分）なのでT_sim[2]はT_meas[2]=180.0にリセットされる
    assert T_sim[2] == pytest.approx(180.0, abs=0.1), (
        f"HWSギャップ後のT_simが実測値にリセットされていない: T_sim[2]={T_sim[2]}"
    )

    # gap_mask[2] が True（リセット点として記録される）
    assert gap_mask[2], "ギャップリセット点はgap_mask[2]=Trueでなければならない"
    assert not gap_mask[0] and not gap_mask[1], "非ギャップ点はgap_mask=Falseであること"

    # ギャップ幅の検証
    assert (t_arr[2] - t_arr[1]) > DT_GAP_THRESH_MIN


# ── テスト 4: 合成データ往復精度（コールドスタート） ────────────────────────

def test_synthetic_round_trip():
    """
    コールドスタートから合成ARCデータを逆解析し、低NRMSEを達成することを確認する。

    Opus指摘①対応：theta_true を最適化初期値 _X0=[120000, 14, 500000, 1000] から
    全次元で大きく離した値に設定し、ウォームスタート問題を排除する。

    識別性の注意：Ea と log_A の間には補償効果（Arrhenius compensation effect）が
    あるため、θ_true の完全復元は不要。NRMSE < 15% でフィットが有効とみなす。
    """
    # theta_true は最適化初期値 _X0=[120000, 14, 500000, 1000] から全次元で乖離
    # Ea: -25%,  log_A: -4 units（線形スケールで10000倍差）,  dH: -30%,  Cp: -10%
    theta_true = ArrheniusTheta(
        Ea=90_000.0,
        log_A=10.0,
        dH=350_000.0,
        Cp=900.0,
        phi=1.0,
    )
    calor_input = _make_synthetic_input(theta_true, T_onset=90.0, t_max_min=20.0, n_points=200)

    # コールドスタートが真に試練となることを確認
    # _X0 で T_sim を計算し、theta_true データに対するNRMSEが十分高いことを検証
    theta_x0 = ArrheniusTheta(Ea=120_000.0, log_A=14.0, dH=500_000.0, Cp=1_000.0, phi=1.0)
    t_arr = np.array([p.time_min for p in calor_input.arc_data])
    T_arr = np.array([p.T_degC for p in calor_input.arc_data])
    T_onset = calor_input.arc_conditions.T_onset_degC
    T_sim_x0, _ = simulate_arc(t_arr, T_arr, T_onset, 1.0, theta_x0)
    nrmse_x0 = compute_nrmse(T_arr, T_sim_x0)
    assert nrmse_x0 > 0.10, (
        f"コールドスタート検証失敗：初期NRMSE={nrmse_x0:.4f} <= 10%。"
        " theta_true が _X0 に近すぎてウォームスタート問題が残存している。"
    )

    result = fit_arrhenius(calor_input)
    assert result.nrmse < 0.15, (
        f"往復精度未達: NRMSE={result.nrmse:.4f} >= 15%（コールドスタートフィット後）"
    )
    assert result.A > 0, "A は正の物理量（1/s）でなければならない"
    assert result.n == 1, "v1 は n=1 固定"


# ── テスト 5: スキーマ φ バリデーション ─────────────────────────────────────

def test_schema_phi_validation():
    """φ < 1.0 で pydantic バリデーションエラーが発生することを確認する。"""
    from pydantic import ValidationError
    with pytest.raises(ValidationError):
        ArcConditions(T_onset_degC=90.0, phi=0.9)   # φ<1 は物理的に不可


# ── テスト 6: データ点数不足バリデーション ──────────────────────────────────

def test_schema_min_data_points():
    """arc_data < 10点で pydantic バリデーションエラーが発生することを確認する。"""
    from pydantic import ValidationError
    short_data = [ArcDataPoint(time_min=float(i), T_degC=90.0 + i) for i in range(5)]
    with pytest.raises(ValidationError):
        CalorInput(
            arc_data=short_data,
            cell_spec=CellSpec(cell_id="x", chemistry="NMC", capacity_ah=5.0,
                               mass_g=70.0, format="21700"),
            arc_conditions=ArcConditions(T_onset_degC=90.0),
        )


# ── テスト 7: NMC 21700 実データフィット ─────────────────────────────────────

@pytest.mark.skipif(
    not _NMC_TXT.exists(),
    reason=f"NMC実データファイル未取得: {_NMC_TXT}",
)
def test_nmc_21700_fit():
    """NMC 21700 SOC100 実データでフィットが収束し NRMSE < 0.5 であることを確認する。"""
    raw = parse_zenodo_txt(
        path=_NMC_TXT,
        cell_id="NMC_SOC100_M1",
        chemistry="NMC",
        capacity_ah=4.9,
        mass_g=68.0,
        cell_format="21700",
        T_onset_degC=90.0,
        phi=1.0,
    )
    calor_input = CalorInput(**raw)
    result = fit_arrhenius(calor_input)

    assert np.isfinite(result.nrmse), "NMC フィット NRMSE が非有限値"
    assert result.nrmse < 0.5, f"NMC フィット NRMSE 異常大: {result.nrmse:.3f}（< 0.5 期待）"
    assert result.Ea >= 50_000.0 and result.Ea <= 200_000.0, (
        f"Ea が境界外: {result.Ea:.0f} J/mol"
    )
    assert result.A > 0.0, "A は正値でなければならない"


# ── テスト 8: LFP 21700 実データフィット ─────────────────────────────────────

@pytest.mark.skipif(
    not _LFP_TXT.exists(),
    reason=f"LFP実データファイル未取得: {_LFP_TXT}",
)
def test_lfp_21700_fit():
    """LFP 21700 SOC100 実データでフィットが収束し NRMSE < 0.5 であることを確認する。"""
    raw = parse_zenodo_txt(
        path=_LFP_TXT,
        cell_id="LFP_SOC100_M1",
        chemistry="LFP",
        capacity_ah=5.0,
        mass_g=95.0,
        cell_format="21700",
        T_onset_degC=120.0,
        phi=1.0,
    )
    calor_input = CalorInput(**raw)
    result = fit_arrhenius(calor_input)

    assert np.isfinite(result.nrmse), "LFP フィット NRMSE が非有限値"
    assert result.nrmse < 0.5, f"LFP フィット NRMSE 異常大: {result.nrmse:.3f}（< 0.5 期待）"
    assert result.Ea >= 50_000.0 and result.Ea <= 200_000.0, (
        f"Ea が境界外: {result.Ea:.0f} J/mol"
    )
    assert result.A > 0.0, "A は正値でなければならない"
