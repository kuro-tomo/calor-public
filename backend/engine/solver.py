"""
固定ステップ RK4 ソルバー（Fornax solver.py から転用・設計書§5.3準拠）

状態ベクトルの次元に依存しない実装のため、Calor 2元状態 y=[T, α] にそのまま適用可能。
"""
import time
from typing import Callable

import numpy as np


def rk4_step(
    f: Callable[[float, np.ndarray], np.ndarray],
    t: float,
    x: np.ndarray,
    dt: float,
) -> np.ndarray:
    """
    1ステップ固定刻み RK4 積分。

    Args:
        f:  状態方程式 dx/dt = f(t, x)
        t:  現在時刻 [s]
        x:  状態ベクトル
        dt: 刻み幅 [s]

    Returns:
        x_next: 次ステップの状態ベクトル
    """
    k1 = f(t,          x)
    k2 = f(t + dt / 2, x + (dt / 2) * k1)
    k3 = f(t + dt / 2, x + (dt / 2) * k2)
    k4 = f(t + dt,     x + dt       * k3)
    return x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def rk4_substep(
    f: Callable[[float, np.ndarray], np.ndarray],
    t: float,
    x: np.ndarray,
    dt: float,
    n_substep: int = 4,
) -> tuple[np.ndarray, float]:
    """
    dt を n_substep に分割して RK4 積分する（精度・安定性向上版）。

    Calor用途では n_substep=40 を指定（設計書§5.3）。

    Args:
        f:         状態方程式 dx/dt = f(t, x)
        t:         現在時刻 [s]
        x:         状態ベクトル
        dt:        総積分時間 [s]
        n_substep: 分割数

    Returns:
        (x_next, elapsed_sec)
    """
    t_wall_start = time.monotonic()
    dt_sub = dt / n_substep
    x_curr = x.copy()
    t_curr = t
    for _ in range(n_substep):
        x_curr = rk4_step(f, t_curr, x_curr, dt_sub)
        t_curr += dt_sub
    elapsed = time.monotonic() - t_wall_start
    return x_curr, elapsed
