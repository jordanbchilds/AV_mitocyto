using Pkg;

using Statistics, Random, NamedArrays, DataStructures, CSV, DataFrames, GLM, StatsBase, Distributions, GLM;

include("./functions_bigHierarchy.jl")

full_data_prepped = CSV.read("../Data_prepped.csv", DataFrame) ;

mitochan = "VDAC"
channels = ["NDUFB8", "CYB", "MTCO1"]
nChan = length(channels)

ctrlIDs = ["C01", "C02", "C03", "C04"]

data = full_data_prepped[:,["ID", "patient_id", "VDAC", "NDUFB8", "CYB", "MTCO1"]]
rename!(data, [:fibreID, :sampleID, :VDAC, :NDUFB8, :MTCO1, :CYB])

data[!,[mitochan; channels]] = log.( data[!,[mitochan; channels]] )

data[!,:sbjType] = [id in ctrlIDs ? "control" : "patient" for id in data[:,:sampleID]]

data[:,:uniqueID] = lpad.(string.(data[:,:fibreID]), 5, "0") .*"_" .*data[:,:sampleID] ;

sbjIDs = unique(data[:,:sampleID])
ptsIDs = sbjIDs[ .!(in.(sbjIDs, Ref(ctrlIDs))) ]
nPts = length(ptsIDs);

if !ispath("Output_BH")
    mkdir("Output_BH")
end

nChains = 5

Threads.@threads for chan in channels
   # Threads.@threads for chain in 1:nChains
        output = gibbs_sampler(data, "VDAC", chan, warmup=10000000, iter=50000)
        mySaver(output, fileRoot="Output_BH/"*chan*"_chain"*lpad(chain, 2, "0")*"_" )
    # end
end

@time gibbs_sampler(data, "VDAC", "MTCO1", warmup=9000, iter=1000)