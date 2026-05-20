test_that("simulate_emergence returns a data.frame with expected shape", {
    n <- 50
    gen_cdf <- pgamma(0:40, shape = 2.71, rate = 0.177)
    gen_cdf <- gen_cdf / max(gen_cdf)
    death_cdf <- pgamma(0:40, shape = 4.42, rate = 0.388)
    death_cdf <- death_cdf / max(death_cdf)

    x <- emerge:::simulate_emergence(
        R = rep(1.5, n),
        k = 1,
        generation_cdf = gen_cdf,
        outcome_cdfs = list(deaths = death_cdf),
        seed = 1
    )

    expect_s3_class(x, "data.frame")
    expect_named(x, c("infections", "deaths"))
    expect_equal(nrow(x), n)
})
