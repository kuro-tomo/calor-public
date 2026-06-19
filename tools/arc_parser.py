"""
Zenodo ARC .txt データ → calor-input v1 JSON 変換ツール

対応フォーマット（Zenodo 21700/18650 データセット）：
  ヘッダ行: "time / min\ttemperature / °C\ttemperature rate / (°C/min)"
  データ行: <float>\t<float>\t<float or blank>

使用法（CLIモジュール実行）：
  python -m tools.arc_parser data/raw/NMC_SOC100_M1.txt \\
      --cell-id NMC_SOC100_M1 \\
      --chemistry NMC \\
      --capacity-ah 4.9 \\
      --mass-g 68.0 \\
      --format 21700 \\
      --T-onset 90.0 \\
      --phi 1.0 \\
      --output data/processed/NMC_SOC100_M1.json

T-onset が未指定の場合、最初のデータ点温度を自動採用する。
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Optional


def _infer_chemistry_from_filename(stem: str) -> str:
    """ファイル名から化学系を推定する（例：NMC_SOC100_M1 → NMC）。"""
    for chem in ("NMC811", "NCA-HE", "NCA-HP", "NCA", "NMC", "LFP"):
        if chem.upper() in stem.upper():
            return chem
    return "UNKNOWN"


def _infer_format_from_path(path: Path) -> str:
    """パスから外形規格を推定する（例：zenodo.7707929/ → 21700）。"""
    parts = str(path).upper()
    if "21700" in parts:
        return "21700"
    if "18650" in parts:
        return "18650"
    return "UNKNOWN"


def parse_zenodo_txt(
    path: Path,
    cell_id: str,
    chemistry: Optional[str] = None,
    capacity_ah: float = 5.0,
    mass_g: float = 70.0,
    cell_format: str = "21700",
    T_onset_degC: Optional[float] = None,
    phi: float = 1.0,
    test_date: Optional[str] = None,
) -> dict:
    """
    Zenodo ARC .txt ファイルを calor-input v1 JSON（dict）に変換する。

    Args:
        path:          入力ファイルパス
        cell_id:       セル識別子
        chemistry:     電池化学系（None なら filename から自動推定）
        capacity_ah:   定格容量 [Ah]
        mass_g:        セル質量 [g]
        cell_format:   外形規格（18650・21700 等）
        T_onset_degC:  発熱開始温度 [°C]（None なら最初のデータ点温度を採用）
        phi:           熱容量補正係数（≥1）
        test_date:     試験日（YYYY-MM-DD・任意）

    Returns:
        calor-input v1 JSON dict
    """
    if chemistry is None:
        chemistry = _infer_chemistry_from_filename(path.stem)

    arc_data: list[dict] = []

    with open(path, encoding="latin-1") as fh:
        lines = fh.readlines()

    # ヘッダ行をスキップ（"time" を含む行）
    data_lines: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if re.search(r"time", stripped, re.IGNORECASE):
            continue
        data_lines.append(stripped)

    for line in data_lines:
        parts = line.split("\t")
        if len(parts) < 2:
            parts = line.split()
        if len(parts) < 2:
            continue

        try:
            t_min = float(parts[0])
            T_degC = float(parts[1])
        except ValueError:
            continue

        dTdt: Optional[float] = None
        if len(parts) >= 3 and parts[2].strip():
            try:
                dTdt = float(parts[2])
            except ValueError:
                pass

        arc_data.append(
            {
                "time_min": t_min,
                "T_degC": T_degC,
                "dTdt_degC_min": dTdt,
            }
        )

    if len(arc_data) < 10:
        raise ValueError(
            f"データ点数不足（{len(arc_data)}点）。calor-input v1 は最低10点必要。"
        )

    if T_onset_degC is None:
        T_onset_degC = arc_data[0]["T_degC"]

    result = {
        "schema_version": "calor-input-v1",
        "arc_data": arc_data,
        "vent_gas": None,
        "cell_spec": {
            "cell_id": cell_id,
            "chemistry": chemistry,
            "capacity_ah": capacity_ah,
            "mass_g": mass_g,
            "format": cell_format,
        },
        "arc_conditions": {
            "T_onset_degC": T_onset_degC,
            "phi": phi,
            "test_date": test_date,
        },
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Zenodo ARC .txt → calor-input v1 JSON 変換"
    )
    parser.add_argument("input", type=Path, help="入力 .txt ファイルパス")
    parser.add_argument("--cell-id", default=None, help="セル識別子（省略時：ファイル名）")
    parser.add_argument("--chemistry", default=None, help="電池化学系（省略時：自動推定）")
    parser.add_argument("--capacity-ah", type=float, default=5.0, help="定格容量 [Ah]")
    parser.add_argument("--mass-g", type=float, default=70.0, help="セル質量 [g]")
    parser.add_argument("--format", dest="cell_format", default=None, help="外形規格")
    parser.add_argument("--T-onset", type=float, default=None, help="発熱開始温度 [°C]")
    parser.add_argument("--phi", type=float, default=1.0, help="熱容量補正係数φ（≥1）")
    parser.add_argument("--test-date", default=None, help="試験日（YYYY-MM-DD）")
    parser.add_argument(
        "--output", "-o", type=Path, default=None, help="出力 JSON ファイルパス（省略時：stdout）"
    )
    args = parser.parse_args()

    input_path: Path = args.input
    if not input_path.exists():
        print(f"エラー：ファイルが見つかりません: {input_path}", file=sys.stderr)
        sys.exit(1)

    cell_id = args.cell_id or input_path.stem
    cell_format = args.cell_format or _infer_format_from_path(input_path)

    try:
        result = parse_zenodo_txt(
            path=input_path,
            cell_id=cell_id,
            chemistry=args.chemistry,
            capacity_ah=args.capacity_ah,
            mass_g=args.mass_g,
            cell_format=cell_format,
            T_onset_degC=args.T_onset,
            phi=args.phi,
            test_date=args.test_date,
        )
    except (ValueError, OSError) as e:
        print(f"エラー：{e}", file=sys.stderr)
        sys.exit(1)

    json_str = json.dumps(result, ensure_ascii=False, indent=2)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json_str, encoding="utf-8")
        print(f"出力完了：{args.output}（{len(result['arc_data'])}点）")
    else:
        print(json_str)


if __name__ == "__main__":
    main()
