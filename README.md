# Calor — Open Arrhenius Thermal Runaway Calibration

Open inverse-analysis workflow for fitting Arrhenius kinetic parameters
(Ea, A, ΔH, Cp) to Accelerating Rate Calorimetry (ARC) data.

**Preprint:** [whitepaper_1_draft.md](docs/whitepaper_1_draft.md)
(submitted to engrXiv)

**Repository:** https://github.com/kuro-tomo/calor-public

---

## What this is

Calor calibrates the two-state adiabatic ODE:

```
dα/dt = A · (1 − α) · exp(−Ea / (R · T_K))
dT/dt = (ΔH / (φ · Cp)) · dα/dt
```

against ARC time-series data using a three-stage optimizer
(L-BFGS-B → Nelder-Mead → default-value fallback).

Applied to five 21700 cell chemistries from the public Zenodo ARC dataset
(DOI: [10.5281/zenodo.7707929](https://zenodo.org/records/7707929), CC BY 4.0),
the workflow achieves NRMSE ≤ 0.263 for all five chemistries.

---

## Benchmark results (Table 1)

| Chemistry | Ea [kJ/mol] | A [s⁻¹] | ΔH [J/kg] | Cp [J/(kg·K)] | ΔTad [K] | NRMSE |
|-----------|------------|---------|-----------|---------------|----------|-------|
| NMC (21700) | 96.3 | 4.55 × 10⁷ | 617,500 | 903 | **684** | **0.075** |
| LFP (21700) | 120.0 | 6.00 × 10⁸ | 500,000† | 1,000† | 500† | 0.263 |
| NCA-HEI (21700) | 121.1 | 1.02 × 10⁶ | 500,000† | 1,500† | 333† | 0.254 |
| NCA-HEII (21700) | 120.5 | 2.99 × 10¹⁰ | 492,700 | 1,084 | **455** | 0.100 |
| NCA-HP (21700) | 120.0 | 2.25 × 10¹⁰ | 500,000† | 1,000† | 500† | 0.089 |

† ΔH/Cp pinned at start value (unconstrained by data under the Arrhenius compensation effect).
Bold ΔTad values are data-constrained; others are placeholder estimates.

Full machine-readable results: [`data/processed/w4_benchmark.json`](data/processed/w4_benchmark.json)

---

## Reproduce

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Run the benchmark (fits all 5 cells, writes data/processed/w4_benchmark.json)
cd calor-public
python -m tools.bench_w4

# 3. Generate Figure 1
python -m tools.plot_w4

# 4. Run tests
pytest
```

Requirements: Python 3.11+, scipy 1.17.1, lmfit 1.3.4, numpy 2.4.4

---

## Data provenance

ARC raw data (`data/raw/*.txt`) are sourced from:

> Ohneseit, S.; Finster, P.; Floras, C.; Lubenau, N.; Uhlmann, N.; Seifert, H. J.; Ziebert, C.
> *Exothermal data from thermal safety assessment of type 21700 lithium-ion batteries
> with NMC, NCA and LFP cathodes by means of Accelerating Rate Calorimetry (ARC).*
> Karlsruhe Institute of Technology (KIT) / Queen's University, Zenodo, 2023.
> DOI: [10.5281/zenodo.7707929](https://zenodo.org/records/7707929) — CC BY 4.0

---

## Repository structure

```
calor-public/
├── backend/
│   └── engine/
│       ├── arrhenius.py      # core ODE + optimizer
│       ├── schema.py         # calor-input v1 JSON schema (Pydantic)
│       └── solver.py
├── tools/
│   ├── arc_parser.py         # Zenodo .txt → calor-input JSON
│   ├── bench_w4.py           # batch benchmark runner
│   └── plot_w4.py            # Figure 1 generator
├── data/
│   ├── raw/                  # ARC .txt files (Ohneseit et al., CC BY 4.0)
│   └── processed/
│       └── w4_benchmark.json # Table 1 source data
├── docs/
│   ├── whitepaper_1_draft.md
│   └── figures/
│       └── fig1_benchmark_fits.png
└── requirements.txt
```

---

## License

Code: MIT License — © 2026 Shinonome Engineering LLC

ARC data (`data/raw/`): CC BY 4.0 — Ohneseit et al. / KIT / Queen's University

---

## Citation

If you use this code or benchmark results, please cite:

> Shinonome Engineering LLC. *Open Validation of Arrhenius Thermal Runaway Calibration
> for Lithium-Ion Cells — A Benchmark Against Public ARC Datasets.*
> engrXiv, 2026. [DOI to be added upon publication]

---

*Shinonome Engineering LLC — mikasa2564@gmail.com*

> **Note:** This repository will be transferred to the Shinonome Engineering LLC GitHub organization upon its creation.
