
# Rank-Ordered Logit from Scratch, with Unobserved Ranking Heterogeneity

This repository contains two companion documents that build a Rank-Ordered
Logit (ROL) model entirely from scratch in base R (no `apollo`, `mlogit`,
or other choice-modeling packages required for estimation) and use it to
illustrate a well-known limitation of the ROL model documented by
Fok, Paap, and van Dijk (2010, *Journal of Applied Econometrics*).

## What's in here

- **`ROL_from_scratch.Rmd`** lays out the ROL model
  mathematically (log-likelihood, gradient, Hessian, Newton–Raphson update,
  standard errors), then applies it to simulated data.
- **A companion `.R` script**: the same code as in the Rmd, meant to be run
  interactively/standalone (e.g. for the full 10,000-replication Monte Carlo
  study, which is impractical to re-run every time the Rmd is knitted).

## What the code does

1. **Simulates data** matching the design in Fok, Paap & van Dijk (2010):
   1,000 individuals rank 4 alternatives, but individuals differ in how many
   of their top items they can *actually* rank correctly, the rest of their
   ranking is filled in at random. This "ranking ability" is unobserved by
   the researcher.
2. **Estimates a standard ROL model** on this data: coded manually
   (log-likelihood, analytical gradient and Hessian, Newton–Raphson
   optimization, standard errors from the inverse Hessian) rather than via
   an existing package, so every step of the estimation is transparent.
3. **Runs a Monte Carlo study** (10,000 replications) to see how well the
   estimated parameters recover the true values on average.

## The main takeaway

Because the standard ROL model assumes every respondent can rank *all*
alternatives correctly, applying it to data where some respondents can't
actually do so is misspecified, and the Monte Carlo results show the
resulting parameter estimates are biased. This motivates the latent-class
extension proposed in Fok, Paap & van Dijk (2010), which lets ranking
ability itself be inferred from the data rather than assumed.

## Reproducing the results

The Monte Carlo loop (10,000 replications) is slow to re-run, so its output
is saved to disk (`mc_summary_ROL`) and simply loaded and printed when the
Rmd is knitted, rather than recomputed. To regenerate it from scratch, run
the (commented-out) Monte Carlo section in the `.R` script directly.


An important part of the code has been written using Claude Code. 
