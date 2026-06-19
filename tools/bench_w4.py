"""
W4 ベンチマーク：5化学系 21700セル Arrheniusパラメータ一括フィット

データソース：Zenodo DOI 10.5281/zenodo.7707929（CC BY 4.0）
出力：
  - コンソール：サマリーテーブル（Ea, A, ΔH, Cp, NRMSE, 経過時間）
  - data/processed/w4_benchmark.json：機械可読結果（白書#1用）

使用法：
  cd /Users/MBP/projects/calor
  python -m tools.bench_w4
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

from backend.engine.arrhenius import fit_arrhenius
from backend.engine.schema import CalorInput
from tools.arc_parser import parse_zenodo_txt

_DATA_DIR = Path(__file__).parents[1] / "data" / "raw"
_OUT_DIR = Path(__file__).parents[1] / "data" / "processed"

# ── ベンチマーク対象定義 ──────────────────────────────────────────────────────
# 注：capacity_ah・mass_g は公称メタ情報であり、Arrheniusフィットの入力ではない
#     （フィットは T(t) 時系列を使用、ΔH/Cp は per-kg 単位）。Table 1 結果には影響しない。
CELLS = [
    {
        "file":       "NMC_SOC100_M1.txt",
        "label":      "NMC（21700）",
        "chemistry":  "NMC",
        "capacity_ah": 4.9,
        "mass_g":     68.0,
        "T_onset":    90.0,
    },
    {
        "file":       "LFP_SOC100_M1.txt",
        "label":      "LFP（21700）",
        "chemistry":  "LFP",
        "capacity_ah": 5.0,
        "mass_g":     95.0,
        "T_onset":    120.0,
    },
    {
        "file":       "NCA_HEI_SOC100_M1.txt",
        "label":      "NCA-HEI（21700）",
        "chemistry":  "NCA-HEI",
        "capacity_ah": 4.8,
        "mass_g":     68.0,
        "T_onset":    95.0,
    },
    {
        "file":       "NCA_HEII_SOC100_M1.txt",
        "label":      "NCA-HEII（21700）",
        "chemistry":  "NCA-HEII",
        "capacity_ah": 4.8,
        "mass_g":     68.0,
        "T_onset":    90.0,
    },
    {
        "file":       "NCA_HP_SOC100_M1.txt",
        "label":      "NCA-HP（21700）",
        "chemistry":  "NCA-HP",
        "capacity_ah": 4.0,
        "mass_g":     65.0,
        "T_onset":    100.0,
    },
]


def _fmt(val: float, fmt: str = ".3g") -> str:
    return f"{val:{fmt}}"


def main() -> None:
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []

    print("=" * 80)
    print("W4 ベンチマーク：5化学系 21700セル Arrheniusパラメータ逆解析")
    print(f"データソース：Zenodo DOI 10.5281/zenodo.7707929")
    print("=" * 80)
    print(
        f"{'化学系':<18}  {'Ea[kJ/mol]':>10}  {'A[1/s]':>12}  "
        f"{'ΔH[J/kg]':>12}  {'Cp[J/kgK]':>10}  {'NRMSE':>8}  {'収束':>4}  {'時間[s]':>7}"
    )
    print("-" * 80)

    total_start = time.monotonic()

    for cell in CELLS:
        path = _DATA_DIR / cell["file"]
        if not path.exists():
            print(f"  ⚠ スキップ（ファイル未取得）: {cell['file']}")
            continue

        label = cell["label"]
        print(f"  {label:<16}  フィット中...", end="", flush=True)
        t0 = time.monotonic()

        try:
            raw = parse_zenodo_txt(
                path=path,
                cell_id=Path(cell["file"]).stem,
                chemistry=cell["chemistry"],
                capacity_ah=cell["capacity_ah"],
                mass_g=cell["mass_g"],
                cell_format="21700",
                T_onset_degC=cell["T_onset"],
                phi=1.0,
            )
            calor_input = CalorInput(**raw)
            result = fit_arrhenius(calor_input)
            elapsed = time.monotonic() - t0

            row = {
                "label":      label,
                "file":       cell["file"],
                "chemistry":  cell["chemistry"],
                "T_onset_degC": cell["T_onset"],
                "Ea_J_mol":   result.Ea,
                "A_1_s":      result.A,
                "dH_J_kg":    result.dH,
                "Cp_J_kgK":   result.Cp,
                "nrmse":      result.nrmse,
                "method":     result.method,
                "converged":  result.converged,
                "n_evals":    result.n_evals,
                "elapsed_s":  elapsed,
                "ci_lower_95": result.ci_lower_95,
                "ci_upper_95": result.ci_upper_95,
            }
            results.append(row)

            # コンソール出力（\rで上書き）
            print(
                f"\r  {label:<16}  "
                f"{result.Ea/1000:>10.1f}  "
                f"{result.A:>12.3e}  "
                f"{result.dH:>12.0f}  "
                f"{result.Cp:>10.0f}  "
                f"{result.nrmse:>8.4f}  "
                f"{'✓' if result.converged else '△':>4}  "
                f"{elapsed:>7.1f}"
            )

        except Exception as exc:
            elapsed = time.monotonic() - t0
            print(f"\r  {label:<16}  ERROR: {exc}")
            results.append({
                "label": label,
                "file": cell["file"],
                "chemistry": cell["chemistry"],
                "error": str(exc),
                "elapsed_s": elapsed,
            })

    total_elapsed = time.monotonic() - total_start
    print("-" * 80)
    print(f"  {'合計':>16}  {total_elapsed/60:.1f}分")
    print("=" * 80)

    # JSON出力
    out_path = _OUT_DIR / "w4_benchmark.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "schema": "calor-w4-benchmark-v1",
                "dataset": "Zenodo DOI 10.5281/zenodo.7707929",
                "license": "CC BY 4.0",
                "cells": results,
            },
            fh,
            ensure_ascii=False,
            indent=2,
        )
    print(f"\n結果保存完了：{out_path}")
    print(f"（白書#1 Table 1 のソースデータとして使用）")


if __name__ == "__main__":
    main()
