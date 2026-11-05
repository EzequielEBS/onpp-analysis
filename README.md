# onpp-analysis

Repository for analysis of the ordered normalized power prior (ONPP): simulation studies and real-data applications comparing NPP (normalized power prior), ONPP (order-restricted NPP), and HLNNPP (hierarchical logistic-normal NPP, a coupled-but-unordered alternative) for borrowing strength from historical data in Bernoulli, Normal, linear-regression, and Poisson models.

## Repository structure

```
code/
  models/            Stan model files (.stan), one subfolder per data family
    bin/               Bernoulli models: eta_bin (NPP), gamma_bin (ONPP),
                        logitnorm_bin / logitnorm_simple_bin (HLNNPP), hier_beta_bin
    lm/                Linear-regression analogues of the above
    normal/            Normal (fixed-variance) NPP/ONPP models
    poi/               Poisson NPP/ONPP models (partial/unfinished)
  sim_exp/           Simulation-study scripts, one subfolder per data family
    bin/, lm/, normal/, poi/
      simulations_*.R / analysis_*.r      run + summarize the main Order x
                                           Discrepancy scenario grid
      analysis_seq_*.r                    sequential-borrowing variant
      prior_comparison_sample.R / _plot.R NPP vs ONPP vs HLNNPP comparison
      logitnorm_aux_*_plot.R              HLNNPP hyperparameter diagnostics
      get_summary_*.r                     MSE/WIS/coverage summary helpers
    cancellation_density_plots.R          illustrates the "cancellation effect"
                                           scenario shared across families
  poster/            Standalone figures built for poster/slide use (combined
                      cross-experiment comparisons, schematic illustrations)
  real_data/         NPP/ONPP/HLNNPP comparison on real clinical-trial datasets
                      (pediatric SLE trial, HOVON AML trial)
  aux_fun_sim.R      Shared helpers (Stan-fit wrappers, plotting theme,
                      gt-table styling) used across all of the above

results/
  samples/           Cached Stan fits (.qs2, not tracked -- regenerate by
                      running the corresponding sampling script)
  sim_data/          Cached simulated datasets (.qs2, tracked via Git LFS)
  figures/           Generated plots, mirroring code/'s per-family layout
  summary_tables/    Generated gt summary tables (.png)
```

Model files are compiled with `cmdstanr::cmdstan_model()` the first time a sampling script runs; the resulting `.exe` binaries are gitignored and rebuilt locally.

## Running the pipeline

All scripts are run from the project root (`Rscript code/...` or sourced
interactively with that as the working directory). For a given family (e.g.
`bin`), the usual order is:

1. `code/sim_exp/bin/simulations_bin.R` -- fit the main scenario grid, caching
   draws to `results/samples/bin/`.
2. `code/sim_exp/bin/analysis_bin.r` -- build the MSE/WIS/coverage tables and
   figures from that cache.
3. `code/sim_exp/bin/prior_comparison_sample.R` then `prior_comparison_plot.R`
   -- the NPP vs. ONPP vs. HLNNPP comparison (separate cache/script pair, since
   it fits an additional model).
4. `code/sim_exp/bin/logitnorm_aux_plot.R` -- HLNNPP hyperparameter diagnostics,
   reusing the cache from step 3.

The `lm` and `normal` families follow the same pattern; `poi` currently only
has `simulations_poi.R`/`analysis_poi.r`.

`code/real_data/real_data_sample.R` then `real_data_plot.R` run the same
NPP/ONPP/HLNNPP comparison on real datasets instead of simulated scenarios.
`code/poster/` scripts are standalone and only need the caches the relevant
`sim_exp`/`real_data` sampling scripts already produced.

## References

- Bracher, J., Ray, E. L., Gneiting, T., and Reich, N. G. (2021). Evaluating epidemic forecasts in an interval format. *PLoS Computational Biology*, 17(2):e1008618.
- Han, Z., Zhang, Q., Tiwari, R., and Bai, T. (2025). Bayesian borrowing with multiple heterogeneous historical studies using order restricted normalized power prior. *Statistics in Medicine*, 44(3-4):e10302.
- Ibrahim, J. G., Chen, M.-H., et al. (2000). Power prior distributions for regression models. *Statistical Science*, 15(1):46-60.
- Ibrahim, J. G., Chen, M.-H., Gwon, Y., and Chen, F. (2015). The power prior: theory and applications. *Statistics in Medicine*, 34(28):3724-3749.
- Neuenschwander, B., Branson, M., and Spiegelhalter, D. J. (2009). A note on the power prior. *Statistics in Medicine*, 28(28):3562-3566.
