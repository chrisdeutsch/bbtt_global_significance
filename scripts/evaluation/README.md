# Evaluation

These scripts are specific to the `X -> HH -> bbtautau` search. Use at
your own risk.


## Workflow

The scripts:

- `listFailedFits.py`
- `listMissingJobs.py`

list failed fits and missing batch jobs in csv format for retrying /
resubmission.

The script `evaluateRetries.py` evaluates retried toys usign
alternative optimizer settings. It takes an input directory containing
fit results in the form of tarballs,
e.g. `retry_idx6929_m300.tar.gz`. The tarball contains 7 csv files one
for each alternative setting. The logic for picking the 'good fit' is
outlined in the INT note.

The script `mergeRetriedToys.py` merges default setting toys with the
retried toys using alternative optimizer settings.

The main analysis is performed in:

- `evaluateLocalSignificance.py`: Make plots and csv files of the q0
  sampling distributions.
- `evaluateGlobalSignificanceAsymptotics.py`: Estimate the global
  significance using the asymptotic approximation for the local significances.
- `evaluateGlobalSignificanceToys.py`: Estimate the global
  significance using toy experiments for the local significances.






