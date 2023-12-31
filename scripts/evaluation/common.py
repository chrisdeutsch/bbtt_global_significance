import warnings

import numpy as np
import pandas as pd
from scipy.optimize import root_scalar, minimize
from scipy import stats

warnings.simplefilter(action='ignore', category=FutureWarning)
from statsmodels.distributions.empirical_distribution import ECDF


def load_toys(fn):
    df = pd.read_csv(fn, low_memory=False)

    # Convert to proper dtypes
    df = df.astype({
        "uncond_status": "int64",
        "cond_status": "int64",
        "uncond_covQual": "int64",
        "cond_covQual": "int64",
    })

    # Rename 'index' to 'toyindex' to avoid confusion with the index of
    # the dataframe
    if "index" in df.columns:
        df.rename(columns={"index": "toyindex"}, inplace=True)

    # Check that there are no duplicates
    assert not df.duplicated(["toyindex", "mass"]).any()

    # Add failed fit column
    df["failed_fit"] = (df["uncond_status"] != 0) | (df["cond_status"] != 0)

    # Transform q0 to one-sided discovery test statistic
    df.loc[df["muhat"] <= 0, "q0"] = 0.0

    # Add flag indiciating whether all fits for a given experiment
    # were successful
    df["good_toy"] = ~df.groupby("toyindex")["failed_fit"].transform(np.any)

    return df


class ToyPvalueCalculator(object):
    def __init__(self):
        self.q0_toys = {}
        self.q0_ecdf = {}

    def add_q0_distribution(self, mass, arr):
        q0 = np.array(arr).flatten()

        self.q0_toys[mass] = q0
        self.q0_ecdf[mass] = ECDF(q0)

    def get_pval(self, mass, q0):
        mass = np.array(mass)
        q0 = np.array(q0)
        assert mass.shape == q0.shape

        pval = np.zeros_like(q0)
        for m in np.unique(mass):
            pval[mass == m] = 1 - self.q0_ecdf[m](q0[mass == m])

        return pval

    def bootstrap(self, rng):
        bootstrapped = ToyPvalueCalculator()

        for mass in self.q0_toys:
            q0 = self.q0_toys[mass]
            q0_resampled = rng.choice(q0, size=len(q0))
            bootstrapped.add_q0_distribution(mass, q0_resampled)

        return bootstrapped


def binom_mle_interval(k, n, bracket=(0.01, 0.1)):
    llh = lambda p, k, n: k * np.log(p) + (n - k) * np.log(1 - p)
    delta_llh = lambda p: llh(p, k, n) - llh(k / n, k, n)

    root_lo = root_scalar(lambda p: -2 * delta_llh(p) - 1, bracket=[bracket[0], k / n])
    root_hi = root_scalar(lambda p: -2 * delta_llh(p) - 1, bracket=[k / n, bracket[1]])
    assert root_lo.converged and root_hi.converged

    return root_lo.root, root_hi.root


def fit_trial_factor(sig_max):
    # Constants removed
    log_pdf = lambda x, n: np.log(n) + (n - 1) * stats.norm.logcdf(x) - x**2 / 2
    nll = lambda n: -log_pdf(sig_max, n).sum()

    res = minimize(nll, x0=[15.], bounds=[[10., 21.]])
    if not res.success:
        print(res)

    assert res.success

    trial_factor, = res.x

    return trial_factor
