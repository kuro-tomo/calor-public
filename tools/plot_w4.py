"""
W4 ベンチマーク図生成：T_meas vs T_sim（5化学系）

データソース：Zenodo DOI 10.5281/zenodo.7707929
出力：docs/figures/fig1_benchmark_fits.png

使用法：
  cd /Users/MBP/projects/calor
  python -m tools.plot_w4
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))

from backend.engine.arrhenius import simulate_arc, ArrheniusTheta
from tools.arc_parser import parse_zenodo_txt

_DATA_DIR  = Path(__file__).parents[1] / "data" / "raw"
_JSON_PATH = Path(__file__).parents[1] / "data" / "processed" / "w4_benchmark.json"
_FIG_DIR   = Path(__file__).parents[1] / "docs" / "figures"
_FIG_DIR.mkdir(parents=True, exist_ok=True)

CELLS = [
    {"file": "NMC_SOC100_M1.txt",    "chemistry": "NMC",     "capacity_ah": 4.9, "mass_g": 68.0, "T_onset": 90.0,  "color": "#2166ac"},
    {"file": "LFP_SOC100_M1.txt",    "chemistry": "LFP",     "capacity_ah": 5.0, "mass_g": 95.0, "T_onset": 120.0, "color": "#4dac26"},
    {"file": "NCA_HEI_SOC100_M1.txt","chemistry": "NCA-HEI", "capacity_ah": 4.8, "mass_g": 68.0, "T_onset": 95.0,  "color": "#d01c8b"},
    {"file": "NCA_HEII_SOC100_M1.txt","chemistry":"NCA-HEII","capacity_ah": 4.8, "mass_g": 68.0, "T_onset": 90.0,  "color": "#f1a340"},
    {"file": "NCA_HP_SOC100_M1.txt", "chemistry": "NCA-HP",  "capacity_ah": 4.0, "mass_g": 65.0, "T_onset": 100.0, "color": "#01665e"},
]


def main() -> None:
    # JSONから最適パラメータを読み込む
    with open(_JSON_PATH, encoding="utf-8") as fh:
        bench = json.load(fh)

    params_by_chem = {c["chemistry"]: c for c in bench["cells"]}

    fig = plt.figure(figsize=(14, 9))
    fig.suptitle(
        "Fig. 1 — ARC Calibration Fits: Measured vs. Simulated Temperature\n"
        "for Five 21700 Cell Chemistries (Ohneseit et al. [4], Zenodo 10.5281/zenodo.7707929)",
        fontsize=11, y=0.98,
    )

    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.32)
    axes_pos = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
    # 6番目のセルは凡例用（row1, col2）
    ax_legend = fig.add_subplot(gs[1, 2])
    ax_legend.axis("off")

    for idx, cell in enumerate(CELLS):
        row, col = axes_pos[idx]
        ax = fig.add_subplot(gs[row, col])

        chem = cell["chemistry"]
        p = params_by_chem.get(chem)
        if p is None:
            ax.set_title(f"{chem} — データ未取得")
            continue

        # ARC生データ読み込み
        raw = parse_zenodo_txt(
            path=_DATA_DIR / cell["file"],
            cell_id=Path(cell["file"]).stem,
            chemistry=chem,
            capacity_ah=cell["capacity_ah"],
            mass_g=cell["mass_g"],
            cell_format="21700",
            T_onset_degC=cell["T_onset"],
            phi=1.0,
        )
        t_arr   = np.array([pt["time_min"] for pt in raw["arc_data"]])
        T_meas  = np.array([pt["T_degC"]   for pt in raw["arc_data"]])

        # 最適パラメータで順方向シミュレーション
        log_A = math.log10(p["A_1_s"])
        theta = ArrheniusTheta(
            Ea=p["Ea_J_mol"], log_A=log_A,
            dH=p["dH_J_kg"], Cp=p["Cp_J_kgK"], phi=1.0,
        )
        T_sim, gap_mask = simulate_arc(t_arr, T_meas, cell["T_onset"], 1.0, theta)

        # 実測値（グレー丸）
        ax.plot(t_arr, T_meas, "o", markersize=3, color="gray",
                alpha=0.55, label="$T_{\\mathrm{meas}}$", zorder=2)
        # HWSギャップリセット点（赤抜き丸で判別可能に強調・NRMSE除外）
        if gap_mask.any():
            ax.scatter(
                t_arr[gap_mask], T_meas[gap_mask],
                facecolors="none", edgecolors="red", s=55, linewidths=1.3,
                zorder=4, label="HWS gap-reset (excluded)",
            )
        # シミュレーション値
        ax.plot(t_arr, T_sim, "-", linewidth=2.0, color=cell["color"],
                label=f"$T_{{\\mathrm{{sim}}}}$ (NRMSE={p['nrmse']:.3f})")

        Ea_kJ = p["Ea_J_mol"] / 1000.0
        dTad  = p["dH_J_kg"] / p["Cp_J_kgK"]
        ax.set_title(
            f"{chem} (21700)\n"
            f"$E_a$={Ea_kJ:.0f} kJ/mol, $A$={p['A_1_s']:.2e} s⁻¹\n"
            f"$\\Delta T_{{ad}}$={dTad:.0f} K, NRMSE={p['nrmse']:.3f}",
            fontsize=8.5,
        )
        ax.set_xlabel("Time [min]", fontsize=8)
        ax.set_ylabel("Temperature [°C]", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7, loc="upper left")
        ax.grid(True, alpha=0.3, linewidth=0.5)

    # 凡例エリアに注記
    ax_legend.text(
        0.05, 0.85,
        "Notes:\n"
        "• Gray circles: measured temperature ($T_{\\mathrm{meas}}$)\n"
        "• Colored line: ODE simulation ($T_{\\mathrm{sim}}$)\n"
        "• Open red circles: HWS gap-reset points\n"
        "  (excluded from NRMSE)\n\n"
        "Model: adiabatic 2-state ODE\n"
        "  $d\\alpha/dt = A(1-\\alpha)\\exp(-E_a/RT)$\n"
        "  $dT/dt = (\\Delta H/C_p)\\cdot d\\alpha/dt$\n\n"
        "Optimizer: L-BFGS-B → Nelder-Mead\n"
        "CI: bootstrap ($n$=200, fallback)\n"
        "Data: Ohneseit et al. (2023)\n"
        "  Zenodo 10.5281/zenodo.7707929",
        transform=ax_legend.transAxes,
        fontsize=8, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f5f5f5", edgecolor="#cccccc"),
    )

    out_path = _FIG_DIR / "fig1_benchmark_fits.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"図を保存仕りました：{out_path}")


if __name__ == "__main__":
    main()
