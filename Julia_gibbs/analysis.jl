using Pkg;

using Statistics, Random, NamedArrays, DataStructures, CSV, DataFrames, GLM, StatsBase, Distributions;

include("./functions.jl")

full_data_prepped = CSV.read("../Data_prepped.csv", DataFrame) ;

mitochan = "VDAC"
channels = ["NDUFB8", "CYB", "MTCO1"]
nChan = length(channels)

raw_data = full_data_prepped[:,["ID", "patient_id", "VDAC", "NDUFB8", "CYB", "MTCO1"]] ;

data_lng = stack(raw_data, [:VDAC, :NDUFB8, :CYB, :MTCO1])

data_lng[:, :value] = log.(data_lng[:, :value])

rename!(data_lng, [:fibreID, :sampleID, :Channel, :Value])
data_lng[!,:sampleID] = convert.(String, data_lng[!,:sampleID])

ctrlID = ["C01", "C02", "C03", "C04"]
data_lng[!,:sbjType] = [xx in ctrlID ? "control" : "patient" for xx in data_lng[:,:sampleID]] ;

sbjIDs = unique(data_lng[:,:sampleID])
ctrlIDs = filter(isControl, sbjIDs)
ptsIDs = filter(isNotControl, sbjIDs) ;

if isfile("Output")
    mkdir("Output")
end

for chan in channels
    Threads.@threads for pat in ptsIDs
        root = chan*"_"*pat
        dd = getData_mats(data_lng; mitochan="VDAC", chan=chan, pts=[pat], ctrlID=ctrlIDs)
        Threads.@threads for i in 1:5
           output::Dict{String, NamedArray} = gibbs_sampler(dd, warmup=400000, iter=20000)
           mySaver(output, fileRoot=string("Output/"*root *"_chain_"*lpad(i, 2, "0")*"_") )
        end
    end
end