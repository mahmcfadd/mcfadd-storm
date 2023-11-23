a=Sys.time()
patt="WIN"
tiflistdir="Y:/Maureen/2color Vamp2Hmr/161015 SD446 DIV23 VAMP2405HmrCy3/Extracted 2"
tiflistdir=gsub("/","\\", tiflistdir, fixed=TRUE)
stackdir=list.files(tiflistdir, pattern=patt, full.names = TRUE)
stkname=list.files(tiflistdir,pattern=patt, full.names = FALSE)
stkl=length(stackdir)
for(i in 1:stkl){
  tifodir=stackdir[i]
  tifofiles=list.files(tifodir,pattern = ".tif", all.files = FALSE, full.names = TRUE, include.dirs = FALSE)
  rejid1=grep("c1", tifofiles)
  rejid2=grep("c5",tifofiles)
  rejid=c(rejid1,rejid2)
  tifofiles=tifofiles[-rejid]
  tifochunk=length(tifofiles)/4
  tifolist=gsub("/","\\", tifofiles, fixed=TRUE)
  tifolistfi=tifolist[1:tifochunk]
  tifolistsec=tifolist[(tifochunk+1):(2*tifochunk)]
  tifolistthird=tifolist[(2*tifochunk+1):(3*tifochunk)]
  tifolistfour=tifolist[(3*tifochunk+1):length(tifofiles)]
  tifpath=paste(tiflistdir,stkname[i],sep = "/")
  inamepfirst=paste0(tifpath,"part1.txt")
  inamepsec=paste0(tifpath,"part2.txt")
  inamepthird=paste0(tifpath,"part3.txt")
  inamepfour=paste0(tifpath,"part4.txt")
  write(tifolistfi,file=inamepfirst, sep = "\n")
  write(tifolistsec,file=inamepsec, sep = "\n")
  write(tifolistthird,file=inamepthird, sep = "\n")
  write(tifolistfour,file=inamepfour, sep = "\n")
}
b=Sys.time()
print(c(a,b))