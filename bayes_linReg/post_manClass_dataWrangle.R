# --- COMPARISON TO MANUAL CLASSIFICATIONS
library("data.table")
library("analysis2Dmito")
library("tidyr")

mdat = as.data.frame( fread("../dat_with_class.txt") )

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")

sbjs = unique(mdat$patient_id)
ctrlIDs = grep("C0", sbjs, value=TRUE)
pts = grep("C0", sbjs, value=TRUE, invert=TRUE)

cols = c("blue","red")

for( chan in channels ){
  x = mdat[ mdat$channel==mitochan & mdat$patient_type=="control", "value"]
  y = mdat[ mdat$channel==chan & mdat$patient_type=="control", "value"]
  for( pat in pts ){
    xpat = mdat[ mdat$channel==mitochan & mdat$patient_id==pat, "value"]
    ypat = mdat[ mdat$channel==chan & mdat$patient_id==pat, "value"]
  
  classDown = mdat[ mdat$channel==paste0(chan, "_MCLASS_DEFICIENT") & mdat$patient_id==pat, "value"]
  classUp = mdat[ mdat$channel==paste0(chan, "_MCLASS_DEFICIENT") & mdat$patient_id==pat, "value"]
  
  class = classDown | classUp
  
  print( class != classUp)
  propdef = round(mean(class),3)
  plot(x,y, pch=20, cex=0.5, xlab=mitochan, ylab=chan, 
       xlim=range(c(x,xpat)), ylim=range(c(y,ypat)),
       col=alphaBlack(0.05), 
       main=bquote(.(pat)*"\n"*pi*"="*.(propdef)))
  points(xpat, ypat, col=cols[class+1], pch=20)
  }
}

tmp_wide = pivot_wider( mdat, 
                        id_cols=c("cell_id", "id", "patient_type", "patient_id"), 
                        names_from=channel, 
                        values_from=value )

tmp_wide$NDUFB8_value = tmp_wide$NDUFB8
tmp_wide$MTCO1_value = tmp_wide$MTCO1
tmp_wide$CYB_value = tmp_wide$CYB
tmp_wide$VDAC_value = tmp_wide$VDAC

tmp_wide$NDUFB8 = NULL
tmp_wide$MTCO1 = NULL
tmp_wide$CYB = NULL
tmp_wide$VDAC = NULL

tmp_wide$NDUFB8_diff = as.numeric( tmp_wide$NDUFB8_MCLASS_DEFICIENT | tmp_wide$NDUFB8_MCLASS_OVEREXP)
tmp_wide$MTCO1_diff = as.numeric( tmp_wide$MTCO1_MCLASS_DEFICIENT | tmp_wide$MTCO1_MCLASS_OVEREXP)
tmp_wide$CYB_diff = as.numeric( tmp_wide$CYB_MCLASS_DEFICIENT | tmp_wide$CYB_MCLASS_OVEREXP)

tmp_wide$NDUFB8_MCLASS_DEFICIENT = NULL
tmp_wide$NDUFB8_MCLASS_OVEREXP = NULL
tmp_wide$MTCO1_MCLASS_DEFICIENT = NULL
tmp_wide$MTCO1_MCLASS_OVEREXP = NULL
tmp_wide$CYB_MCLASS_DEFICIENT = NULL
tmp_wide$CYB_MCLASS_OVEREXP = NULL

tmp_wide$cell_id = NULL

tmp_wide_2 = pivot_longer(tmp_wide, cols=!c("id", "patient_type", "patient_id"), names_to=c("channel", ".value"), names_sep="_")

tmp_wide_2$cell_id = paste0(tmp_wide_2$channel, "_", tmp_wide_2$id)


colnames(tmp_wide_2) = c("id", "patient_type", "patient_id", "channel", "value", "jbcClass", "cell_id")

write.table(tmp_wide_2, file="../dat_with_class_prepped.txt", sep="\t", row.names=FALSE)

