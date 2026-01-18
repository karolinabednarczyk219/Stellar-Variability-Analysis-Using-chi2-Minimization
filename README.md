# Stellar Variability Analysis using χ² Minimization

This project implements a simple and transparent method for detecting periodic variability in stellar light curves using χ² minimization.

## Overview
The analysis is based on fitting a sinusoidal model to photometric time-series data for a range of trial periods. For each period, the model parameters are obtained by minimizing the χ² statistic, and the best-fit period is selected as the global minimum of χ²(P).

The project focuses on:
- implementing χ²-based period search from scratch,
- handling heteroscedastic uncertainties via weighted fitting,
- validating detected periods using simple statistical criteria,
- emphasizing methodological clarity rather than optimization.

## Methods
- Weighted mean subtraction
- Linear least-squares fitting of sine and cosine components
- Period grid search with χ² minimization
- Basic variability selection criteria based on reduced χ² and contrast

## Contents
- `analyze_stellar_variability.py` — main analysis script
- `Bednarczyk_Karolina_StellarVariability_Chi2.pdf` — short project report describing the method and results

## Notes
This project was designed as a methodological exercise to understand period-search techniques and their limitations, rather than to produce a high-performance or production-ready pipeline.
