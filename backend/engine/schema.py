"""
calor-input v1 JSONスキーマ（pydantic v2）

設計書§2モジュール構成・§8セキュリティ設計準拠。
物理制約バリデーション：φ≥1・容量>0・質量>0。
"""
from __future__ import annotations

from typing import List, Literal, Optional

from pydantic import BaseModel, Field


class ArcDataPoint(BaseModel):
    """ARC実測データ1点"""
    time_min: float = Field(..., description="経過時間 [min]")
    T_degC: float = Field(..., description="セル表面温度 [°C]")
    dTdt_degC_min: Optional[float] = Field(
        default=None,
        description="昇温速度 [°C/min]（来歴保存用・逆解析では未使用）",
    )


class VentGas(BaseModel):
    """ベントガス組成（UL9540A:2025試験報告書より・機能C用）"""
    H2_vol_frac: float = Field(..., ge=0.0, le=1.0, description="H2体積分率")
    CO_vol_frac: float = Field(..., ge=0.0, le=1.0, description="CO体積分率")
    CO2_vol_frac: float = Field(..., ge=0.0, le=1.0, description="CO2体積分率")
    CH4_vol_frac: float = Field(..., ge=0.0, le=1.0, description="CH4体積分率")
    total_mass_g: float = Field(..., gt=0.0, description="総ベントガス質量 [g]")


class CellSpec(BaseModel):
    """セル仕様メタデータ"""
    cell_id: str = Field(..., description="セル識別子（例：NMC_SOC100_M1）")
    chemistry: str = Field(..., description="電池化学系（例：NMC811・LFP・NCA-HE）")
    capacity_ah: float = Field(..., gt=0.0, description="定格容量 [Ah]")
    mass_g: float = Field(..., gt=0.0, description="セル質量 [g]")
    format: str = Field(..., description="外形規格（例：18650・21700・pouch）")


class ArcConditions(BaseModel):
    """ARC試験条件"""
    T_onset_degC: float = Field(
        ..., description="発熱開始温度 [°C]（α=0の初期条件）"
    )
    phi: float = Field(
        default=1.0,
        ge=1.0,
        description="熱容量補正係数φ（≥1。理想断熱ARC=1.0）",
    )
    test_date: Optional[str] = Field(default=None, description="試験日（YYYY-MM-DD）")


class CalorInput(BaseModel):
    """calor-input v1 JSONスキーマ（トップレベル）"""
    schema_version: Literal["calor-input-v1"] = "calor-input-v1"
    arc_data: List[ArcDataPoint] = Field(
        ...,
        min_length=10,
        description="ARC実測時系列データ（発熱相のみ・最低10点）",
    )
    vent_gas: Optional[VentGas] = Field(
        default=None,
        description="ベントガス組成（機能C用）。欠落時は機能CスキップしA成果物のみCVP収録",
    )
    cell_spec: CellSpec
    arc_conditions: ArcConditions
