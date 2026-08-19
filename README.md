# BARIA Motor - Ballistics

Experimental data reduction and internal-ballistics analysis of the BARIA
solid rocket motor: 27 pressure traces (9 batches x 3 pressure levels),
Bayern-Chemie data reduction, Vieille's law fit, characteristic velocity,
and Monte Carlo uncertainty propagation through a 0D ballistic model.

## Contributors

- [Matteo Mafrici](https://github.com/matteomafrici)
- [Mario Guida](https://github.com/marioguida27)
- Martina Lucia Magarelli
- Margherita Palitta

## Results

| Quantity                  | Value                     |
|---------------------------|---------------------------|
| Vieille `a`               | 1.744 ± 0.021 mm s^-1 bar^-n |
| Vieille `n`               | 0.380 ± 0.003             |
| Fit R^2                   | 0.9984                    |
| Characteristic velocity c*| 1511.6 ± 11.5 m/s         |
| 0D model RMSE (mid P)     | 2.35 bar (~5 % P_eff)     |

Full derivation and discussion are in
[report/baria-report.pdf](report/baria-report.pdf).

## Usage

Requires MATLAB (R2026a or later) and a dataset file in `data/`. The
pipeline auto-detects every numeric `N x 3` matrix in the dataset file
(one per firing batch, 1 kHz sampling). See
[data/README.md](data/README.md) for the expected format.

Run from the repository root:

```matlab
matlab -batch "run('matlab/baria_ballistics.m')"
```

The script runs the complete analysis in one pass:

1. BC data reduction on all traces, per-regime statistics
2. Vieille's law OLS fit + characteristic velocity
3. Monte Carlo uncertainty propagation (N = 5000, fixed seed)
4. 0D model validation vs experimental firings, RMSE

Figures are written to `figures/` with a light theme at 200 dpi.

## Repository layout

```
matlab/    analysis script + helper functions
data/      dataset placement + format notes (dataset not distributed)
figures/   generated figures
report/    technical report (PDF)
```

## License

MIT - see [LICENSE](LICENSE).