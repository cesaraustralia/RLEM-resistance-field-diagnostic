using Turing

# We need a logistic function, which is provided by StatsFuns.
using LinearAlgebra
using StatsFuns: logistic
using StatsPlots
using CSV, DataFrames
using Random
d = CSV.read("./data/summary.csv", DataFrame)

@model function rlem(R, pop, time, n, nDead, nGeno, NDead)

    if R === missing
        R = Vector{Union{Missing,Float64}}(undef, max(pop...))
        muR ~ Beta(1, 1)
        sdR ~ Exponential(1)
        for i ∈ eachindex(R)
            R[i] ~ truncated(Normal(muR, sdR), 0, 1)
        end
    end


    # priors
    interceptRR ~ Normal(1, 10)
    # interceptRS ~ Normal(1, 10)
    interceptSS ~ Normal(1, 10)
    slope ~ Normal(0, 10)

    for i ∈ eachindex(time)
        p = R[pop[i]]
        q = 1 - p

        mRR = logistic(interceptRR + slope * time[i])
        # pRS = logistic(interceptRS + slope * time[i]) # assume comlete recessiveness
        mSS = logistic(interceptSS + slope * time[i]) # assume comlete recessiveness

        pDead = min(mRR * p^2 + mSS * 2 * p * q + mSS * q^2, 1)
        nDead[i] ~ Binomial(n[i], pDead)

        scale = mRR * p^2 + mSS * 2 * p * q + mSS * q^2
        NDead[i, :] ~ Multinomial(nGeno[i], [mRR * p^2; mSS * 2 * p * q; mSS * q^2] ./ scale)
    end
end

# test 
# x = sample(rlem(missing, [1; 1; 1], [1; 3; 8], [101; 101; 101], [94; 99; 100], [33; 38; 39], [0 13 20; 0 16 22; 0 16 23]), NUTS(), 1_000)

# real data
NDead = [d.nDeadRR d.nDeadRS d.nDeadSS]
Random.seed!(1234);
nsamp = 5_000
x = sample(rlem(missing, d.pop, d.time, d.n, d.nDead, d.nGeno, NDead), NUTS(), nsamp)
plot(x)


# save table of credible intervals
uniquepops = d[.!nonunique(d[:, [:popName]]), [:popName, :pop]]
function get_popName(parstring)
    reg = match(r"\d", parstring)
    if (typeof(reg) === Nothing)
        return parstring
    end
    index = parse(Int64, match(r"\d", parstring).match)
    return uniquepops.popName[index]
end

credint = DataFrame(hpd(x))
credint.parameters = String.(credint.parameters)
credint.parameters = map(x -> get_popName(x), credint.parameters)

CSV.write("./plots/credible_intervals.csv", credint)

# How many hours should we wait to count the alive mites?
t = [i for i = 1:0.1:8]
pRR = 100 .* logistic.(x[:interceptRR] * ones(1, length(t)) + x[:slope] * reshape(t, 1, length(t)));
pSS = 100 .* logistic.(x[:interceptSS] * ones(1, length(t)) + x[:slope] * reshape(t, 1, length(t)));
pRRmu = mean(pRR, dims=1)'
pSSmu = mean(pSS, dims=1)'
pRRsd = std(pRR, dims=1)'
pSSsd = std(pSS, dims=1)'


pSSmu[t.==1]
quantile(vec(pSS[:, t.==1]), [0.025, 0.975])

pRRmu[t.==1]
quantile(vec(pRR[:, t.==1]), [0.025, 0.975])


plot(t, pRRmu, grid=false, ribbon=(pRRsd, pRRsd), fillalpha=0.5, label="RR",
    ylabel="mortality (%)", xlabel="time (h)");
p = plot!(t, pSSmu, grid=false, ribbon=pSSsd, fillalpha=0.5, label="RS/SS")
savefig("./plots/time_vs_mortality.png")

[[i for i = 1:0.1:8] pRRmu pRRsd pSSmu pSSsd] |>
x -> DataFrame(x, [:time; :pRRmu; :pRRsd; :pSSmu; :pSSsd]) |>
     x -> CSV.write("./plots/fig4.csv", x)



# How many mites should we assay to be 90% confident of 1% resistance

@model function prob_correct_id(pRR, pSS; R=0.30)
    n_cont ~ Uniform(10, 1000)
    n = floor(Int, n_cont)
    R_pop ~ Bernoulli(0.5) # resistant population
    nRR ~ Binomial(n, R_pop * R^2)
    nSS = n - nRR
    survRR ~ Binomial(nRR, 1.0 - pRR)
    survSS ~ Binomial(nSS, 1.0 - pSS)
end

t = 1 # hours exposure
pRR = mean(logistic.(x[:interceptRR] + x[:slope] .* t), dims=1)[1];
pSS = mean(logistic.(x[:interceptSS] + x[:slope] .* t), dims=1)[1];

R = 0.30 # resistance allele genotype frequency
x2 = sample(prob_correct_id(pRR, pSS; R=R), Prior(), 100000)
df = DataFrame(x2)
df.n = round.(Int, df.n_cont);

surv_threshold = 3
df.correct_RR_id = df.R_pop .== 1 .&& ((df.survSS .+ df.survRR) .>= surv_threshold);
df.correct_SS_id = df.R_pop .== 0 .&& ((df.survSS .+ df.survRR) .< surv_threshold);
df.correct_id = df.correct_SS_id .|| df.correct_RR_id;

pdat = combine(groupby(df, :n),
    :correct_id => (x -> sum(x) / length(x)) => :accuracy,
    :correct_RR_id => (x -> sum(x) / length(x)) => :true_RR,
    :correct_SS_id => (x -> sum(x) / length(x)) => :true_SS,
)
plot(
    pdat.n,
    [pdat.accuracy pdat.true_RR pdat.true_SS],
    labels=["accuracy" "true R" "true S"],
    title="resistance threshold ≥ $surv_threshold surviving mites",
    xlabel="Number of mites assayed",
    ylabel="Rate"
)
savefig("./plots/accuracy_for_n$surv_threshold.png")

# make table of mean accuracy for different R and surv_threshold

@model function prob_correct_id(pRR, pSS; df)

    R_pop = Array{Int64}(undef, length(df.R))
    nRR = Array{Int64}(undef, length(df.R))
    nSS = Array{Int64}(undef, length(df.R))
    survRR = Array{Int64}(undef, length(df.R))
    survSS = Array{Int64}(undef, length(df.R))

    for i ∈ eachindex(df.R)
        R = df.R[i]
        n_mites = df.n_mites[i]

        R_pop[i] ~ Bernoulli(0.5) # sample resistant population
        nRR[i] ~ Binomial(n_mites, R_pop[i] * R^2) # sample mites from resistant population 
        nSS[i] = n_mites - nRR[i] # sample mites from susceptible population
        survRR[i] ~ Binomial(nRR[i], 1.0 - pRR) # survival of mites from resistant population
        survSS[i] ~ Binomial(nSS[i], 1.0 - pSS) # survival of mites from susceptible population        

    end

    correct_RR_id = R_pop .== 1 .&& ((survSS .+ survRR) .>= df.surv_threshold)
    correct_SS_id = R_pop .== 0 .&& ((survSS .+ survRR) .< df.surv_threshold)
    correct_id = correct_SS_id .|| correct_RR_id


    return correct_id
end


df = allcombinations(DataFrame,
    R=[0.3],
    surv_threshold=[3],
    n_mites=[10, 100, 300, 1000],
)

df = allcombinations(DataFrame,
    R=[0.01, 0.05, 0.10, 0.20, 0.30, 0.50, 0.80, 1.00],
    surv_threshold=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    n_mites=[1, 3, 10, 30, 100, 300, 1000]
)


mod1 = prob_correct_id(pRR, pSS; df)
mod1()
n_samples = 5000 # pretty slow for numbers higher than 1000
sims = [mod1()[:] for i in 1:n_samples]
meansims = []
for i in 1:length(sims[1])
    push!(meansims, mean([sims[j][i] for j in 1:n_samples]))
end
df.correct_id = meansims
println(df)
CSV.write("./plots/accuracy.csv", df)


# confidence added by each surviving mite
@model function rlem2(R, N, alive)
    p = R
    q = 1 - R
    alive ~ Binomial(N, 1 - (p^2 * pRR + 2 * p * q * pSS + q^2 * pSS))
end


rlem3 = rlem2(0, 100, 1) # initialise model

function get_prob(; R_i=0.30, N_i=100, alive=:0:10)
    [prob"alive = i | model = rlem3, R=R_i, N = N_i" for i = alive]
end
get_prob()

alive_vect = 0:10
N_vect = [30, 100, 300, 1000]
R_vect = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
labels = map(R_i -> string("R = ", R_i), R_vect)
for N_i in N_vect
    plot(
        alive_vect, [get_prob(; N_i=N_i, R_i=R_i) for R_i in R_vect],
        yl="Probability of observing N survivors | N=$N_i",
        xl="N survivors",
        legend=true,
        xticks=0:1:10,
        labels=reshape(labels, 1, length(labels)),
    )
    savefig("./plots/probality_of_n_survivors_for_N$(N_i)_pop.png")
end


N_vect = [30, 100, 300, 1000]
R_vect = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
labels = map(N_i -> string("N = ", N_i), N_vect)
for R_i in R_vect
    plot(
        alive_vect, [get_prob(; N_i=N_i, R_i=R_i) for N_i in N_vect],
        yl="Probability of observing N survivors | R=$R_i",
        xl="N survivors",
        legend=true,
        xticks=0:1:10,
        labels=reshape(labels, 1, length(labels)),
    )
    savefig("./plots/probality_of_n_survivors_for_R$(R_i)_pop.png")
end

[get_prob(; N_i=N_i, R_i=0) for N_i in N_vect] |>
x -> [collect(alive_vect) x...] |>
     x -> DataFrame(x, [:time; :n30; :n100; :n300; :n1000]) |>
          x -> CSV.write("./plots/fig5.csv", x)


