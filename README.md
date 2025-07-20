# Diode Fitting Workflow

This repository contains MATLAB/Octave scripts for fitting diode I-V characteristics. The core program `main.m` loads example data, estimates starting parameters, performs a multi-stage optimisation and optionally allows interactive parameter refinement.

## Setup

1. Clone or download this repository.
2. Start **MATLAB** or **GNU Octave** and add the repository folder to your path, e.g.
   ```matlab
   addpath(genpath('path/to/repo'));
   ```
3. Ensure the optimisation functions are available:
   - In MATLAB, the **Optimization Toolbox** is required for `lsqnonlin`.
   - In Octave, install the `optim` package.

## Running

Run the main routine from the command line:
```matlab
>> main
```
The script prompts whether to load parameters saved from previous runs. It then fits the built-in measurement data, plots the results and asks if you want to save them. Saved files include:
- `fit_results_<timestamp>.mat` – fitting data and parameters
- `fit_plot_<timestamp>.png` – generated figure
- `fit_data_<timestamp>.csv` – exported data table
- `fit_params_<timestamp>.txt` – summary of fitted parameters
Interactive adjustments, when chosen, are written to `adjusted_params_<timestamp>.(mat|txt)`.

## Fitting workflow

`main.m` performs the following steps:
1. Load configuration constants and example I-V data.
2. Optionally initialise parameters from previous results; otherwise use a Lambert‑W based estimate.
3. Call `performFitting` which optimises different voltage regions before a global fit using Levenberg–Marquardt and trust‑region algorithms.
4. Compute diode, ohmic and non‑ohmic current components with `calculateCurrents` and plot them via `plotResults`.
5. Save the results and display final parameters. Interactive refinement can further tweak parameters with real‑time plots.
6. After each global fit the mean and maximum relative errors are checked.
   If either exceeds their thresholds the optimisation is repeated using the
   current parameters as the starting point with increased iteration limits
   until the criteria are met or `config.optimization.max_attempts` is reached.
## Parameter regularisation

`errorFunction` and its partial variants now accept a vector of prior parameter values. When `config.regularization.lambda` is greater than zero an L2 penalty
`lambda * ((x - prior).^2)` is appended to the residuals. Set `config.regularization.prior` and `config.regularization.lambda` in `loadConfig.m` to bias
the optimisation toward expected parameter values.

## Additional physical parameters

The model now includes a second saturation current `J02` describing tunnelling
effects. Typical magnitudes lie between `1e-12` and `1e-3` A depending on the
device. Its ideality factor `n2` is chosen at run time within
`config.physics.n2_range` (default `[2 8]`). If you simply press enter when
prompted, the midpoint of this range is used. The exponential term uses this
value via `A2 = q/(kb*T*n2)`.


## Optimisation attempts

`config.optimization.max_attempts` sets how many times the final optimisation
stage is allowed to run.  Each retry increases the iteration limits and stops
early when the improvement between attempts falls below
`config.optimization.improvement_threshold`.

### Optimisation methods

Set `config.optimization.method` to `multistart`, `particleswarm` or `ga` to
select the global search strategy. `particleswarm` and `ga` require the
corresponding MATLAB toolboxes or Octave packages.

### Additional options

- `improvement_threshold` – minimal improvement (in percent) between retries
  before the optimisation stops early.
- `m_range`/`m_step` – when convergence fails the script sweeps the exponent `m`
  over this range to seek a better fit.


### Parallel computing

Set `config.parallel.use` to `true` in `loadConfig.m` to enable parallel
optimisation. When enabled, `runWithMultiStart` uses a parallel pool for the
MultiStart solver and the fallback random restarts. Specify
`config.parallel.poolSize` to control the size of the pool or leave it empty to
use the default. The Parallel Computing Toolbox (MATLAB) or the `parallel`
package (Octave) must be installed.
