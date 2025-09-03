# VIX Chaos Analysis - Quantitative Study of Financial Volatility

A comprehensive quantitative analysis of the VIX (Volatility Index) using chaos theory and nonlinear dynamics to distinguish between stochastic financial processes and deterministic chaotic systems.

## 🎯 Overview

This project implements rigorous scientific methods to analyze VIX data over 20 years (2004-2024), applying:
- **Takens Embedding** for phase space reconstruction
- **Lyapunov Exponent** calculation (Rosenstein algorithm)
- **Hurst Exponent** analysis (R/S method)
- **Correlation Dimension** (Grassberger-Procaccia)
- **Surrogate Data Testing** (IAAFT method)
- **Lorenz Attractor** comparison for chaos benchmarking

## 📊 Key Results

**VIX vs Lorenz System Comparison:**
- **VIX Lyapunov**: 0.03 ± 0.006 (stochastic process)
- **Lorenz Lyapunov**: 0.923 (deterministic chaos)
- **Surrogate Test**: VIX does not reject stochasticity hypothesis
- **Conclusion**: VIX exhibits stochastic behavior, not deterministic chaos

## 🚀 Usage

```bash
julia vix_robust_analysis.jl
```

The analysis will:
1. Download real VIX data (20 years)
2. Perform quantitative chaos analysis
3. Generate scientific visualization
4. Export detailed results

## 📁 Output Files

- `vix_historic_complete.pdf/.png` - Main scientific visualization
- `VIX_Chaos_Analysis_Results_*.txt` - Quantitative results and parameters

## 🔬 Scientific Methodology

Based on established literature:
- Rosenstein et al. (1993): Lyapunov exponent calculation
- Takens (1981): Embedding theorem for phase space reconstruction
- Theiler et al. (1992): Surrogate data testing methodology
- Grassberger & Procaccia (1983): Correlation dimension analysis

## 📦 Dependencies

```julia
using YFinance, Plots, StatsBase, Statistics, Random, Dates, Printf, Distributions
```

## 🎨 Features

- **Real VIX Data**: 20 years of market volatility data
- **High-Quality Visualization**: Publication-ready figures (300 DPI)
- **Regime Detection**: Automatic volatility threshold detection
- **Phase Space Analysis**: 3D embedding visualization
- **Comparative Analysis**: VIX vs Lorenz system benchmarking
- **Statistical Rigor**: Error bars, confidence intervals, p-values

## 📈 Analysis Panels

1. **Historical Series**: Complete 20-year VIX timeline with regime zones
2. **Bimodal Distribution**: Probability density with KDE smoothing
3. **Phase Space 3D**: VIX embedding with Takens reconstruction
4. **Lorenz Attractor**: Reference chaotic system for comparison

## 🧮 Mathematical Framework

- **Embedding Dimension**: 3D for VIX, 5D for Lorenz
- **Time Delay**: τ=5 days (VIX), τ=2 steps (Lorenz)
- **Evolution Time**: 50 steps (VIX), 200 steps (Lorenz)
- **Integration Method**: Runge-Kutta 4th order (Lorenz)
- **Sampling Rate**: Daily (VIX), dt=0.001 (Lorenz)

---

*This analysis provides scientific evidence that VIX behaves as a stochastic financial process rather than a deterministic chaotic system.*