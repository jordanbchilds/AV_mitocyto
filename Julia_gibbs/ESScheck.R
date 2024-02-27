# --- convergence check
library("tidyr")
library("data.table")
library("coda")
library("mcmcse")
library("stringr")

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

raw_data = read.csv("../Data_prepped.csv", header=TRUE)

raw_data = raw_data[,c("ID", "patient_id", mitochan, channels)]
colnames(raw_data) = c("fibreID", "sampleID", mitochan, channels)

data_lng = pivot_longer(raw_data, cols=c(mitochan, channels), names_to="channel")
data_lng = as.data.frame(data_lng)

sbj = unique(data_lng$sampleID)
ctrlID = grep("C", sbj, value=TRUE)
pts = sbj[!(sbj %in% ctrlID)]

for( chan in channels ){
  for( pat in pts ){
    for( chain in 1:5){
      root = paste0(chan, "_", pat, "_chain", str_pad(ch, width=2, pad="0"))
      
      fn_post = file.path("Output", paste0(root, "_POST.csv"))
      
      post = as.data.frame( fread( fn_post ))
      print( root )
      print(effectiveSize( post ))
      print(multiESS(post))
    }
  }
}







