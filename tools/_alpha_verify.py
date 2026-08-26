"""
白書#1 改訂用：αf・補正NRMSE・絶対RMSE一括算出
既存フィットパラメータを使い順方向積分のみ実行（最適化不要・数秒で完了）
"""
import json, math, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))

from backend.engine.arrhenius import (
    ArrheniusTheta, _ode_rhs, compute_nrmse, DT_GAP_THRESH_MIN, N_SUBSTEP,
)
from backend.engine.solver import rk4_substep
from backend.engine.schema import CalorInput
from tools.arc_parser import parse_zenodo_txt

DATA_DIR  = Path(__file__).parents[1] / "data" / "raw"
JSON_PATH = Path(__file__).parents[1] / "data" / "processed" / "w4_benchmark.json"

CELLS = [
    {"file": "NMC_SOC100_M1.txt",      "chem": "NMC",     "cah": 4.9, "mg": 68.0, "To": 90.0},
    {"file": "LFP_SOC100_M1.txt",      "chem": "LFP",     "cah": 5.0, "mg": 95.0, "To": 120.0},
    {"file": "NCA_HEI_SOC100_M1.txt",  "chem": "NCA-HEI", "cah": 4.8, "mg": 68.0, "To": 95.0},
    {"file": "NCA_HEII_SOC100_M1.txt", "chem": "NCA-HEII","cah": 4.8, "mg": 68.0, "To": 90.0},
    {"file": "NCA_HP_SOC100_M1.txt",   "chem": "NCA-HP",  "cah": 4.0, "mg": 65.0, "To": 100.0},
]

with open(JSON_PATH) as f:
    w4 = json.load(f)
params_map = {c["chemistry"]: c for c in w4["cells"]}

HDR = (f"{'Chem':<10} {'αf':>5} {'ΔTad_K':>8} {'αf·ΔTad':>8} "
       f"{'obs_rise':>9} {'nrmse_old':>10} {'nrmse_new':>10} "
       f"{'absRMSE':>8} {'range_K':>8} {'n_pts':>6}")
print(HDR)
print("-" * len(HDR))

for c in CELLS:
    p = params_map[c["chem"]]
    path = DATA_DIR / c["file"]

    raw = parse_zenodo_txt(
        path=path, cell_id=Path(c["file"]).stem,
        chemistry=c["chem"], capacity_ah=c["cah"],
        mass_g=c["mg"], cell_format="21700",
        T_onset_degC=c["To"], phi=1.0,
    )
    ci = CalorInput(**raw)

    t_all = np.array([pt.time_min for pt in ci.arc_data])
    T_all = np.array([pt.T_degC    for pt in ci.arc_data])

    candidates = np.nonzero(T_all >= c["To"])[0]
    onset_idx  = int(candidates[0]) if len(candidates) > 0 else 0
    t_arr  = t_all[onset_idx:]
    T_meas = T_all[onset_idx:]
    T0     = float(T_meas[0])

    theta = ArrheniusTheta(
        Ea=p["Ea_J_mol"], log_A=math.log10(p["A_1_s"]),
        dH=p["dH_J_kg"], Cp=p["Cp_J_kgK"], phi=1.0,
    )
    f_ode = _ode_rhs(theta)

    n = len(t_arr)
    T_sim    = np.empty(n)
    T_sim[0] = T0
    gap_mask = np.zeros(n, dtype=bool)   # [0] は False のまま（現状の挙動）
    y = np.array([T0, 0.0])

    for i in range(n - 1):
        dt_min = t_arr[i + 1] - t_arr[i]
        if dt_min > DT_GAP_THRESH_MIN:
            y[0] = T_meas[i + 1]
            gap_mask[i + 1] = True
        else:
            dt_sec = dt_min * 60.0
            try:
                y, _ = rk4_substep(f_ode, t_arr[i] * 60.0, y, dt_sec, n_substep=N_SUBSTEP)
            except Exception:
                pass
            y[1] = min(max(float(y[1]), 0.0), 1.0)
        T_sim[i + 1] = y[0]

    alpha_final = float(y[1])

    # 旧 NRMSE（初期化点込み）
    nrmse_old = compute_nrmse(T_meas, T_sim, gap_mask)

    # 新 NRMSE（初期化点除外）
    gm_new = gap_mask.copy()
    gm_new[0] = True
    nrmse_new = compute_nrmse(T_meas, T_sim, gm_new)

    keep_new = ~gm_new
    abs_rmse = float(np.sqrt(np.mean((T_meas[keep_new] - T_sim[keep_new]) ** 2)))
    n_pts    = int(keep_new.sum())

    T_range  = float(T_meas[~gap_mask].max() - T_meas[~gap_mask].min())
    obs_rise = T_range   # ≈ T_onset〜T_peak（断熱自己発熱幅）

    dTad    = p["dH_J_kg"] / p["Cp_J_kgK"]
    af_dTad = alpha_final * dTad

    print(
        f"{c['chem']:<10} {alpha_final:>5.3f} {dTad:>8.0f} {af_dTad:>8.0f} "
        f"{obs_rise:>9.0f} {nrmse_old:>10.4f} {nrmse_new:>10.4f} "
        f"{abs_rmse:>8.1f} {T_range:>8.1f} {n_pts:>6d}"
    )
