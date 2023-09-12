# create a prior function for beta(2, 5). Note that the function compute log of priors and must be
# transformed
pr_fun <- build_prior(dist = "beta", p1 = 2, p2 = 5)$fun
round(exp(pr_fun(c(-1, 0, 0.1, 0.5, 0.9, 1, 2))), 3)

# create a prior function for a Uniform distribution with lower bound of -10, and upper bound of 10
pr_fun <- build_prior(dist = "uniform", p1 = -10, p2 = 10)$fun
round(exp(pr_fun(c(-11, -5, 0, 4, 10, 12))), 3)

# create a prior function for Gaussian distribution with mean of 1 and standard deviation of 2
pr_fun <- build_prior(dist = "normal", p1 = 1, p2 = 2)$fun
round(exp(pr_fun(c(-9, -5, -3, -1, 1, 3, 5, 7, 11))), 3)
