
library("data.table")
library("tidyr")


dat_raw = as.data.frame( fread("./Data_prepped.csv") )

mitochan = "VDAC"
channels = c("CYB", "MTCO1", "NDUFB8")

c("ID", "patient_id", mitochan, channels) %in% colnames(dat_raw)

dat_1 = dat_raw[,c("ID", "patient_id", mitochan, channels)] 

dat_2 = pivot_longer(dat_1, cols=c(mitochan, channels), names_to="channel")

dat_3 = dat_2
dat_3$cell_id = paste0(dat_3$patient_id, "_", dat_3$ID)
dat_3$patient_type = "patient"
dat_3$patient_type[ grep("C0", dat_3$patient_id) ] = "control"

dat_4 = dat_3[,c("value", "ID", "patient_id", "channel", "patient_type", "cell_id")]
colnames(dat_4) =  c("value", "id","patient_id", "channel", "patient_type", "cell_id")

head(dat_4)
write.table(dat_4, file="manual_classification/dat.txt", sep="\t")













