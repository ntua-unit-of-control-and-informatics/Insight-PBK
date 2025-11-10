using Distributions

μ = 1
σ = 10

log_m = log(μ) - 0.5 * log(1 + (σ/μ)^2)
sigma = sqrt(log(1 + (σ/μ)^2))

println("Lognormal parameters for mean=$μ and sd=$σ : log-mean = $log_m , log-sd = $sigma")

d = LogNormal(log_m, sigma)

samples = rand(d, 10000000)

mean(samples), std(samples)