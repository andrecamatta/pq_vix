# Análise de Caos do VIX - Estudo Quantitativo da Volatilidade Financeira

Uma análise quantitativa abrangente do VIX (Índice de Volatilidade) usando teoria do caos e dinâmica não-linear para distinguir entre processos financeiros estocásticos e sistemas caóticos determinísticos.

## 🎯 Visão Geral

Este projeto implementa métodos científicos rigorosos para analisar dados do VIX ao longo de 20 anos (2004-2024), aplicando:
- **Embedding de Takens** para reconstrução do espaço de fase
- **Expoente de Lyapunov** (algoritmo de Rosenstein)
- **Expoente de Hurst** (método R/S)
- **Dimensão de Correlação** (Grassberger-Procaccia)
- **Teste de Dados Substitutos** (método IAAFT)
- **Atrator de Lorenz** como referência para comparação caótica

## 📊 Resultados Principais

**Comparação VIX vs Sistema de Lorenz:**
- **VIX Lyapunov**: 0.03 ± 0.006 (processo estocástico)
- **Lorenz Lyapunov**: 0.923 (caos determinístico)
- **Teste Substituto**: VIX não rejeita hipótese de estocasticidade
- **Conclusão**: VIX exibe comportamento estocástico, não caótico determinístico

## 🚀 Como Usar

```bash
julia vix_robust_analysis.jl
```

A análise irá:
1. Baixar dados reais do VIX (20 anos)
2. Executar análise quantitativa de caos
3. Gerar visualização científica
4. Exportar resultados detalhados

## 📁 Arquivos de Saída

- `vix_historic_complete.pdf/.png` - Visualização científica principal
- `VIX_Chaos_Analysis_Results_*.txt` - Resultados quantitativos e parâmetros

## 🔬 Metodologia Científica

Baseada na literatura estabelecida:
- Rosenstein et al. (1993): Cálculo do expoente de Lyapunov
- Takens (1981): Teorema de embedding para reconstrução do espaço de fase
- Theiler et al. (1992): Metodologia de teste com dados substitutos
- Grassberger & Procaccia (1983): Análise da dimensão de correlação

## 📦 Dependências

```julia
using YFinance, Plots, StatsBase, Statistics, Random, Dates, Printf, Distributions
```

## 🎨 Características

- **Dados Reais do VIX**: 20 anos de dados de volatilidade do mercado
- **Visualização de Alta Qualidade**: Figuras prontas para publicação (300 DPI)
- **Detecção de Regimes**: Detecção automática de limiares de volatilidade
- **Análise do Espaço de Fase**: Visualização de embedding 3D
- **Análise Comparativa**: Benchmarking VIX vs sistema de Lorenz
- **Rigor Estatístico**: Barras de erro, intervalos de confiança, p-valores

## 📈 Painéis de Análise

1. **Série Histórica**: Timeline completa de 20 anos do VIX com zonas de regime
2. **Distribuição Bimodal**: Densidade de probabilidade com suavização KDE
3. **Espaço de Fase 3D**: Embedding do VIX com reconstrução de Takens
4. **Atrator de Lorenz**: Sistema caótico de referência para comparação

## 🧮 Framework Matemático

- **Dimensão de Embedding**: 3D para VIX, 5D para Lorenz
- **Atraso Temporal**: τ=5 dias (VIX), τ=2 passos (Lorenz)
- **Tempo de Evolução**: 50 passos (VIX), 200 passos (Lorenz)
- **Método de Integração**: Runge-Kutta 4ª ordem (Lorenz)
- **Taxa de Amostragem**: Diária (VIX), dt=0.001 (Lorenz)

---

*Esta análise fornece evidência científica de que o VIX se comporta como um processo financeiro estocástico em vez de um sistema caótico determinístico.*