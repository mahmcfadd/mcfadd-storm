library(readr)
library(data.table)

DIR="D:/Maureen/STORM/Final analysis_R scripts/test diana/REcs"


file=list.files(path=DIR, pattern="chsort",full.names = TRUE)
filename=list.files(path=DIR, pattern="chsort",full.names = FALSE)
mydat=read_csv(file)

framechunk= max(mydat$frame)/40
resfile=list.files(path = DIR,pattern = "drift",full.names = TRUE)
driftxy=read.delim(resfile)
driftslicecorr=(driftxy$Slice*framechunk)
driftxy[,2]=driftxy[,2]*10
driftxy[,3]=driftxy[,3]*10
fstart=framechunk
mydatcorr=mydat
for(i in 1:39){
  fend=fstart+framechunk
  mydatcorr[which(mydat$frame>=fstart&mydat$frame<fend),3]=mydat[which(mydat$frame>=fstart&mydat$frame<fend),3]+driftxy[i,2]
  mydatcorr[which(mydat$frame>=fstart&mydat$frame<fend),4]=mydat[which(mydat$frame>=fstart&mydat$frame<fend),4]+driftxy[i,3]
  fend=fend+framechunk
  fstart=fstart+framechunk
}
mydatcorr$`uncertainty_z [nm]` =gsub("Inf", "NaN", mydatcorr$`uncertainty_z [nm]`)

mydatV=mydatcorr[which(mydatcorr$Channel==1),]
mydatH=mydatcorr[which(mydatcorr$Channel==2),]

colpath=gsub(pattern = "chsort.csv", replacement = ".csv",filename)
bpath=paste0(DIR,"/","405ch_",colpath)
rpath=paste0(DIR,"/","Cy3ch_",colpath)
write_csv(mydatV,bpath)
write_csv(mydatH,rpath)
