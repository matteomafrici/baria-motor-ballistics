# Data

The plots and results in this repository are produced from a dataset that
is not published.

To run the analysis, place your dataset here as a MATLAB `.mat` file, then
run from the repository root:

```matlab
matlab -batch "run('matlab/baria_ballistics.m')"
```

## Expected dataset format

The `.mat` file must contain one variable per firing batch, each an
`N x 3` matrix sampled at 1 kHz (1 ms per row):

| Column | Regime   | D_t [mm] |
|--------|----------|----------|
| 1      | Low P    | 28.80    |
| 2      | Mid P    | 25.26    |
| 3      | High P   | 21.81    |

Batches are processed in alphabetical order of their variable names. The
pipeline auto-detects every numeric `N x 3` variable in the file, so the
number of batches is not fixed.

For a different motor, adjust the constants in the CONFIG section at the
top of `matlab/baria_ballistics.m` (grain geometry, web, throat
diameters, propellant composition).

The generated figures and the technical report are the published outputs
of the analysis.