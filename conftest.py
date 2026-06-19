"""pytest ルート設定：プロジェクトルートを sys.path に追加する。"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
