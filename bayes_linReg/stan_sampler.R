# install.packages("devtools")
# library("devtools")
# install_github("jordanbchilds/analysis2Dmito", force=TRUE, quiet=TRUE)

library(analysis2Dmito)
library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")
library("rstan")

folder = "stan_sampler_gamma00001"

dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("Output", folder), showWarnings = FALSE)

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

raw_data = read.csv("../Data_prepped.csv", header=TRUE)

raw_data = raw_data[,c("ID", "patient_id", mitochan, channels)]
colnames(raw_data) = c("fibreID", "sampleID", mitochan, channels)

data_lng = pivot_longer(raw_data, cols=c(mitochan, channels), names_to="channel")
data_lng = as.data.frame(data_lng)

data = data_lng
data$value = log(data$value)

sbj = unique(data$sampleID)
ctrlID = grep("C", sbj, value=TRUE)
pts = sbj[!(sbj %in% ctrlID)]

data_list = list()
for (chan in channels) {
  for (i in seq_along(pts)) {
    data_list[[paste(chan, pts[i], sep = "__")]] = getData_mats(
      data,
      pts = pts[i],
      channels = c(mitochan, chan),
      ctrlID = ctrlID,
      getIndex = TRUE
    )
  }
}

# ncores = detectCores()-2
ncores = 8
cl  = makeCluster(ncores)
{
  output = parLapply(
    cl,
    data_list,
    stan_inference, 
    warmup=50000, 
    iter=52000,
    parameterVals = list(tau_def=0.0001))
}
stopCluster(cl)

for( rt in names(output) ){
  list_saver(output[[rt]], file.path("Output", folder, rt), rootSep="_")
}




