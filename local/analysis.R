library(emerge)
library(data.table)
library(ggplot2)

# Ingredients
# Time since first case
# Rt series
# Overdispersion parameter
# Generation time
# Delay from infection to death
# DFR

Deaths = 88 # deaths observed by end of epidemic
Death_delay = pgamma(0:40, shape = 4.42, rate = 0.388)
Death_delay = Death_delay / max(Death_delay)
gen_shape = 2.71  # WHO Ebola Response Team 2014 (NEJM): SI mean 15.3, SD 9.3
gen_rate  = 0.177
Generation_delay = pgamma(0:80, shape = gen_shape, rate = gen_rate)
Generation_delay = Generation_delay / max(Generation_delay)

# Median R(t) from a doubling time via the discrete Euler-Lotka equation:
#   1 = R * sum_j g(j) exp(-r j),  r = ln(2) / T_d
# where g is the discrete generation-interval pmf derived from its CDF.
R_from_doubling = function(T_d, gen_cdf) {
    pmf = diff(gen_cdf)
    delays = seq_along(pmf) - 1
    r = log(2) / T_d
    1 / sum(pmf * exp(-r * delays))
}

doubling_time = 14  # days
R_mu = R_from_doubling(doubling_time, Generation_delay)
cfr_a = 55       # Beta prior on CFR: alpha
cfr_b = 169-55   # Beta prior on CFR: beta
# curve(dbeta(x, cfr_a, cfr_b), 0, 1)
K = 10000
seed = 0
i = 1
results = vector("list", K)
R_series = vector("list", K)
weights = numeric(K)
simlengths = integer(K)
D_sims = numeric(K)
total_infections = numeric(K)
set.seed(42)
system.time(
while (i <= K) {
    # TODO: draw simulation length from customizable prior
    simlength = sample(20:300, 1)

    # Generate R(t) sequence
    R = emerge:::random_Rt_m32(simlength, R_mu, 0.1, 10, seed)

    # Generate infection & death series
    x = emerge:::simulate_emergence(R = R, k = 0.1,
        generation_cdf = Generation_delay,
        outcome_cdfs = list(
            deaths_max = Death_delay
        ),
        seed = seed + 1)

    # Total simulated deaths by simulation end assuming 100% CFR
    D_sim = sum(x$deaths_max)

    # Only accept if enough potential deaths to ever match the observation
    if (D_sim >= Deaths) {
        # Beta-binomial marginal likelihood (CFR integrated out under beta(cfr_a, cfr_b))
        results[[i]] = x
        R_series[[i]] = R
        weights[i] = exp(lchoose(D_sim, Deaths) +
                         lbeta(Deaths + cfr_a, D_sim - Deaths + cfr_b) -
                         lbeta(cfr_a, cfr_b))
        simlengths[i] = simlength
        D_sims[i] = D_sim
        total_infections[i] = sum(x$infections)
        i = i + 1
    }
    seed = seed + 2
})

cat("Efficiency:", K/(seed/2), "\n")
# results2 = data.table::rbindlist(results, idcol = "run")
# results2[, t := (200 - .N + 1):200, by = run]
# results2[, alpha := weights[run] / max(weights)]
#
# results2[alpha>0.01, sum(infections), by = run][, hist(V1)]
# ggplot(results2[alpha > 0.01]) + geom_line(aes(x = t, y = infections, group = run, alpha = alpha))
# ggplot(results2[alpha > 0.01]) + geom_line(aes(x = t, y = deaths_max, group = run, alpha = alpha))

# Posterior samples: resample by weights, then draw CFR from its conjugate posterior
M = 10000
idx = sample.int(K, size = M, replace = TRUE, prob = weights)
posterior = data.table(
    cfr = rbeta(M, Deaths + cfr_a, D_sims[idx] - Deaths + cfr_b),
    simlength = simlengths[idx],
    total_infections = total_infections[idx]
)

# Posterior of R(t): subsample resampled indices, align trajectories to end at t=0
n_show = 300
R_dt = rbindlist(lapply(sample(idx, n_show), function(i) {
    R = R_series[[i]]
    data.table(plot_id = i, t = -(length(R) - 1):0, R = R)
}), idcol = "draw")
ggplot(R_dt) +
    geom_line(aes(x = t, y = R, group = draw), alpha = 0.08) +
    geom_hline(yintercept = R_mu, colour = "grey40", linetype = "dashed") +
    labs(x = "Days before present", y = "R(t)",
         title = sprintf("Posterior R(t) (dashed = prior median, T_d = %g d)",
                         doubling_time))

# Posterior of CFR (with Beta(cfr_a, cfr_b) prior shown for comparison)
ggplot(posterior, aes(x = cfr)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    stat_function(fun = function(c) dbeta(c, cfr_a, cfr_b),
                  colour = "grey40", linetype = "dashed") +
    labs(x = "CFR", y = "density",
         title = "Posterior of CFR (dashed = prior)")

# Posterior of time since emergence (uniform 100:200 prior shown for comparison)
ggplot(posterior, aes(x = simlength)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = 1, fill = "steelblue", alpha = 0.5) +
    geom_hline(yintercept = 1/101, colour = "grey40", linetype = "dashed") +
    labs(x = "Time since emergence (days)", y = "density",
         title = "Posterior of time since emergence (dashed = prior)")

# Posterior of total infections to date
ggplot(posterior, aes(x = total_infections)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    labs(x = "Total infections", y = "density",
         title = "Posterior of total infections to date")

