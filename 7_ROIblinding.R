#Libraries
library(svDialogs)
library(readr)
library(data.table)
library(caroline)

#Method
ROIdirs=vector()
Randnames=vector()
export=vector()
Rd=vector()
Rtot=vector()
  
folderdir="D:/Maureen/STORM/2colorVamp2Hmr/Vamp2-Homer1 DZcal/SD414DIV22"
blinddir=paste(folderdir,"Blinded",sep = "/")
dirrename=dir.exists(blinddir)

while(dirrename==TRUE){
  blinddir=paste0(blinddir,"New")
  dirrename=dir.exists(blinddir)
}

dir.create(blinddir)

fno=dlgInput(message="Indicate number of ROI folders to blind", default="6")$res
fno=as.numeric(fno)

for(i in 1:fno){
  ROIdir=dlgDir(default = folderdir, "Select ROI directory")$res
  ROIdirs=c(ROIdirs,ROIdir)
  Randname=paste0("Blind",i)
  Randnames=c(Randnames,Randname)
  export=c(export,"NaN")
  Rd=c(Rd,"Nan")
  Rtot=c(Rtot,"Nan")
  
}

rand=runif(fno,min = 1, max=fno)
ord=order(rand)
ROIdirsx=ROIdirs[ord]
Blindres=cbind(Randnames,ROIdirsx,export,Rd,Rtot)
brespath=paste(blinddir,"Blindtable.xls",sep="/")
write.delim(Blindres,brespath, sep ="\t")
print(ROIdirs)