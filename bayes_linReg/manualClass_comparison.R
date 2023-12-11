# --- COMPARISON TO MANUAL CLASSIFICATIONS
library("data.table")
library("analysis2Dmito")

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

tmp_wide_2 = pivot_longer( tmp_wide, )

tmp_wide$mclass = tmp_wide$NDUFB8_MCLASS_DEFICIENT | tmp_wide$NDUFB8_MCLASS_OVEREXP

newdat = mdat[mdat$channel %in% c(mitochan, chan), ]

newdat$V1 = NULL
newdat$mclass = 0






