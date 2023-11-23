library(readr)
library(data.table)

bluech=vector()
redch=vector()
bluech[1]=1
redch[1]=4
for(i in 1:20000){
  bluech[i+1]=bluech[i]+6
  redch[i+1]=redch[i]+6
}

DIR="D:/Maureen/STORM/Final analysis_R scripts/test diana/REcs"
IJname="IMG1"
partlist=list.files(DIR, pattern=IJname, full.names = TRUE)
partlist=grep(pattern=c("merge"), partlist,value=TRUE)
prot=grep(pattern=c("protocol"), partlist,value=FALSE)
partlist=partlist[-prot]

i=1
bluerecf=vector()
redrecf=vector()

frameno=vector()

for(d in 1:4){
  partno=paste0("part",d)
  
  part=partlist[grep(partno, partlist)]
  mydatVp1=read_csv(part)
  
  breci=mydatVp1$frame%in%bluech
  bluerec=mydatVp1[which(breci=="TRUE"),]
  
  rreci=mydatVp1$frame%in%redch
  redrec=mydatVp1[which(rreci=="TRUE"),]
  
  bluerecf=rbind(bluerecf,bluerec)
  redrecf=rbind(redrecf,redrec)
  
}

bluerecf$`uncertainty_z [nm]` =gsub("Inf", "NaN", bluerecf$`uncertainty_z [nm]`)
redrecf$`uncertainty_z [nm]`=gsub("Inf", "NaN", redrecf$`uncertainty_z [nm]`)

# fbname=paste0("405ch_",IJname)
# fbname=paste0(fbname,".csv")
# summpathblue=paste( DIR ,fbname, sep = "/")
# write_csv(bluerecf,summpathblue)
# 
# frname=paste0("Cy3ch_",IJname)
# frname=paste0(frname,".csv")
# summpathred=paste(DIR,frname, sep = "/")
# write_csv(redrecf,summpathred)

datVl=nrow(bluerecf)
datHl=nrow(redrecf)
chV=(1:datVl)
chH=(1:datHl)
chV[1:datVl]=1
chH[1:datHl]=2
mydatV=cbind(bluerecf,chV)
mydatH=cbind(redrecf,chH)
setnames(mydatV, "chV", "Channel")
setnames(mydatH, "chH", "Channel")
mydat=rbind(mydatV,mydatH)
sortname=paste0(IJname, "chsort.csv")
sortpath=paste(DIR, sortname, sep="/")
write_csv(mydat,sortpath)
# }

# names=c("id","frame", "x [nm]", "y [nm]", "z [nm]", "sigma1 [nm]", "sigma2 [nm]", "intensity [photon]", "offset [photon]", "bkgstd [photon]", "chi2", "uncertainty_xy [nm]", "uncertainty_z [nm]", "detections")
# setnames(bluerec,names)
# setnames(redrec,names)


