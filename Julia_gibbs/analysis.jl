using Pkg
Pkg.add("NBInclude")

using Statistics, Random, NamedArrays, DataStructures, CSV, DataFrames, GLM, StatsBase, Distributions, NBInclude;

include("./functions.jl")

full_data_prepped = CSV.read("../Data_prepped.csv", DataFrame);

mitochan = "VDAC"
channels = ["NDUFB8", "CYB", "MTCO1"]
nChan = length(channels)

raw_data = full_data_prepped[:,["ID", "patient_id", "VDAC", "NDUFB8", "CYB", "MTCO1"]] ;

data_lng = stack(raw_data, [:VDAC, :NDUFB8, :CYB, :MTCO1])

data_lng[:, :value] = log.(data_lng[:, :value])

rename!(data_lng, [:fibreID, :sampleID, :Channel, :Value])

ctrlID = ["C01", "C02", "C03", "C04"]
data_lng[!,:sbjType] = [xx in ctrlID ? "control" : "patient" for xx in data_lng[:,:sampleID]] ;

isControl(sampleID::String) = (sampleID in ctrlID)
isControl(sampleID::String3) = (sampleID in ctrlID)
isNotControl(sampleID::String) = !(sampleID in ctrlID)
isNotControl(sampleID::String3) = !(sampleID in ctrlID);

sbj = unique(data_lng[:,:sampleID])
pts = filter(isNotControl, sbj) ;

mkdir("Output") ;

for chan in channels
    Threads.@threads for pat in pts
        root = chan*"_"*pat
        dd = getData_mats(data_lng; mitochan="VDAC", chan=chan, pts=[pat])
        output::Dict{String, NamedArray} = gibbs_sampler(dd, warmup=100, iter=50)
        mySaver(output, fileRoot=string("Output/"*root*"__") )
        
        # Threads.@threads for i in 1:1
        #    output::Dict{String, NamedArray} = gibbs_sampler(dd, warmup=100, iter=50)
        #    mySaver(output, fileRoot=string("Output/"*root *"_chain_"*lpad(i, 2, "0")*"_") )
        # end
    end
end

Threads.nthreads()

nbexport("analysis.jl", "analysis.ipynb")