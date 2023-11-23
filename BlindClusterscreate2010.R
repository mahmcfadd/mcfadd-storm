##Libraries
library(dbscan)
library(rgl)
library(alphashape3d)
library(FNN)
library(geometry)
library(data.table)
library(heplots)
library(cluster)
library(readr)
library(svDialogs)
library(caroline)
library(stats)
#Variables
set.seed(23)
iname="ST001"
mpts=10
epsi=500
steepcut=0.0001
epcut=50
alpha=30
rad=125
totroi=0
skipp=0


#Functions
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Method
#ROIdir=dlgDir(default= "C:/Users/Maureen/Documents/STORM", "Select ROI directory")

Btabledir="D:/Maureen/STORM/2colorBsnHmr/BHSD414RECsp1/BlindedNew"
Btablefile=paste(Btabledir,"Blindtable.xls",sep = "/")
Btable=read.delim(Btablefile)
Bno=nrow(Btable)

Clusterdone=0
clusttype=dlgList(c("Vamp2","Bassoon"),title = "Select the type of staining to cluster:")$res

blindn=as.numeric(dlgInput(message = "Enter a folder number to analyse", default="1")$res)
if(clusttype=="Vamp2") vollim=10000000 else vollim=2000000
if(is.na(Btable[blindn,4])==FALSE) {if(Btable[blindn,4]==Btable[blindn,5]) Clusterdone=1}

bfolder=paste(Btabledir,Btable[blindn,1],sep = "/")
bstat=dir.exists(bfolder)
if(bstat==FALSE)dir.create(bfolder)
fstart=Btable[blindn,4]
if(is.na(Btable[blindn,4])==TRUE) fstart=1  else fstart=fstart+1

ROIdir=as.character(Btable[blindn,2])
FileVamp=list.files(ROIdir,pattern = "405_TS3D", all.files = FALSE, full.names = TRUE, include.dirs = FALSE)
FileCy3=list.files(ROIdir,pattern = "Cy3_TS3D", all.files = FALSE, full.names = TRUE, include.dirs = FALSE)
NameVamp=list.files(ROIdir,pattern = "405_TS3D", all.files = FALSE, full.names = FALSE, include.dirs = FALSE)
NameHomer=list.files(ROIdir,pattern = "Cy3_TS3D", all.files = FALSE, full.names = FALSE, include.dirs = FALSE)
modnameV=gsub(iname, "file", NameVamp,fixed = TRUE)
modnameV=gsub("405_TS3D", "file",modnameV,fixed = TRUE)
modnameV=gsub(".csv", ".xls",modnameV,fixed = TRUE)
modnameH=gsub(iname, "file", NameHomer,fixed = TRUE)
modnameH=gsub("Cy3_TS3D", "file",modnameH,fixed = TRUE)
modnameH=gsub(".csv", ".xls",modnameH,fixed = TRUE)
flength=length(modnameV)
bfoldername=basename(bfolder)
Clusterdir=paste0(bfoldername,"Clusters")
Clusterdir=paste(bfolder,Clusterdir,sep = "/")
dir.create(Clusterdir)
Btable[blindn,3]=Clusterdir
Btable[blindn,5]=flength

for(filen in fstart:flength) { ##start file loop
  if(Clusterdone==1) break
  # filen=1
  Resclusters=vector()
  channelinfo=vector()
  mydatVo=read.csv(FileVamp[filen])
  # mydatV=read.csv("Y:/Maureen/STORMVamp2Hmr Old Quant/2color Vamp2 Hmr/Treated/160524 SD426 VampHmr/Veh1/ST002ROIs/ST002ROI20405_TS3D.csv")
  nameh=grep(modnameV[filen], modnameH, value=FALSE)
  mydatHo=read.csv(FileCy3[nameh])
  # mydatH=read.csv("Y:/Maureen/STORMVamp2Hmr Old Quant/2color Vamp2 Hmr/Treated/160524 SD426 VampHmr/Veh1/ST002ROIs/ST002ROI20Cy3_TS3D.csv")
  mydatV=mydatVo[,2:4]
  mydatH=mydatHo[,2:4]
  open3d()
  rgl.bg(color="black")
  rgl.points(mydatV[,1:3],col="Red", alpha=0.3)
  rgl.points(mydatH[,1:3],col="Green", alpha=0.3)
  
  ##Vamp2 Cluster determination
  op=optics(mydatV, epsi, minPts = mpts)
  
  ###Randomized distribution
  Vn=nrow(mydatV)
  Vx=mydatV[,1]
  Vy=mydatV[,2]
  Vz=mydatV[,3]
  randx=runif(Vn, min = min(Vx), max= max(Vx) )
  randy=runif(Vn, min = min(Vy), max= max(Vy) )
  randz=runif(Vn, min = min(Vz), max= max(Vz) )
  randdist=cbind(randx,randy,randz)
  
  oprand=optics(randdist, epsi, minPts = mpts)
  randreachmed=median(oprand$coredist[which(oprand$coredist!=Inf)])
  if(is.na(randreachmed)==TRUE) next
  randreachmad=mad(oprand$coredist[which(oprand$coredist!=Inf)])
  
  ##Bassoon preliminary cluster
  simfact=5
  refine=1
  
  while(refine==1){
    maincut=randreachmed-(simfact*randreachmad)
    resrand=extractDBSCAN(oprand,maincut)
    resrandclust=cbind(randdist,resrand$cluster)
    res=extractDBSCAN(op, maincut)
    concatV=cbind(mydatV, res$cluster, res$order,res$reachdist,res$coredist)
    setnames(concatV, "res$coredist", "Coredist")
    setnames(concatV, "res$reachdist","ReachDist")
    setnames(concatV, "res$order", "Order")
    setnames(concatV, "res$cluster", "Cluster")
    maxCV=max(concatV$Cluster)
    concatVtemp=concatV[which(concatV$Cluster!=0),]
    centVm=vector()
    clids=vector()
    concatVtemp=as.matrix(concatVtemp)
    cvpoints=1
    
    for(i in 1:maxCV){
      nr=length(concatVtemp[which(concatVtemp[,4]==i),4])
      if(nr>50){
        aV=ashape3d(concatVtemp[which(concatVtemp[,4]==i),1:3],alpha=1)
        calpha=2*mean(aV$triang[,6])
        aV=ashape3d(concatVtemp[which(concatVtemp[,4]==i),1:3],alpha=calpha)
        VolV=volume_ashape3d(aV)
        #concatVtemp=concatVtemp[which(concatVtemp[,4]!=i),]
        if(VolV>vollim){
          open3d()
          rgl.bg(color="black")
          rgl.points(concatVtemp[which(concatVtemp[,4]==i),1:3], col="coral",alpha=0.3)
          rgl.points(mydatV[,1:3],col="Gray", alpha=0.2)
          rgl.points(mydatH[,1:3],col="White", alpha=0.2)
          title3d(main = i, col="White")
          centv=c(mean(concatV[which(concatV$Cluster==i),1]),mean(concatV[which(concatV$Cluster==i),2]),mean(concatV[which(concatV$Cluster==i),3]), i, VolV)
          centVm=rbind(centVm,centv)
          clids[cvpoints]=i
          cvpoints=cvpoints+1
        }
      }
    }
    
    Vclustchoice=dlgList(c("Unrefine",clids,"Refine","None"), multiple=FALSE, title="Select Bassoon cluster")
    Mainclust=Vclustchoice$res
    if(Mainclust=="Unrefine") simfact=simfact-0.4
    if(Mainclust!="Refine"& Mainclust!="Unrefine") refine=0 
    if(Mainclust=="Refine") simfact=simfact+0.2
    dlist=rgl.dev.list()
    dlistl=length(dlist)
    if(dlistl>0){
      for(i in 1:dlistl){
        rgl.close()
      }
    }
  }
 
  if(Mainclust!="None"){
    #if(Mainclust=="None") next()
    Mainclust=as.numeric(Vclustchoice$res)
    mainVclust=concatV[which(concatV$Cluster==Mainclust),]
    
    ##Homer1 main clusters
    op=optics(mydatH, epsi, minPts = mpts)
    
    ###Randomized distribution for Homer1 main cluster
    Hn=nrow(mydatH)
    Hx=mydatH[,1]
    Hy=mydatH[,2]
    Hz=mydatH[,3]
    randx=runif(Hn, min = min(Hx), max= max(Hx) )
    randy=runif(Hn, min = min(Hy), max= max(Hy) )
    randz=runif(Hn, min = min(Hz), max= max(Hz) )
    randdist=cbind(randx,randy,randz)
    oprand=optics(randdist, epsi, minPts = mpts)
    
    randreachmed=median(oprand$coredist[which(oprand$coredist!=Inf)])
    randreachmad=mad(oprand$coredist[which(oprand$coredist!=Inf)])
    
    ###Homer1 preliminary main cluster
    simfact=5
    refine=1
    
    while(refine==1){
      maincut=randreachmed-(simfact*randreachmad)
      resrand=extractDBSCAN(oprand,maincut)
      resrandclust=cbind(randdist,resrand$cluster)
      res=extractDBSCAN(op, maincut)
      concatH=cbind(mydatH, res$cluster, res$order,res$reachdist,res$coredist)
      setnames(concatH, "res$coredist", "Coredist")
      setnames(concatH, "res$reachdist","ReachDist")
      setnames(concatH, "res$order", "Order")
      setnames(concatH, "res$cluster", "Cluster")
      maxCH=max(concatH$Cluster)
      concatHtemp=concatH[which(concatH$Cluster!=0),]
      concatHtemp=as.matrix(concatHtemp)
      centHm=vector()
      clids=vector()
      cvpoints=1
      
      for(i in 1:maxCH){
        nr=length(concatHtemp[which(concatHtemp[,4]==i),4])
        if(nr>50){
          aH=ashape3d(concatHtemp[which(concatHtemp[,4]==i),1:3],alpha=1)
          calpha=2*mean(aH$triang[,6])
          aH=ashape3d(concatHtemp[which(concatHtemp[,4]==i),1:3],alpha=calpha)
          VolH=volume_ashape3d(aH)
          #concatHtemp=concatHtemp[which(concatHtemp[,4]!=i),]
          if(VolH>2000000){
            open3d()
            rgl.bg(color="black")
            rgl.points(concatHtemp[which(concatHtemp[,4]==i),1:3], col="chartreuse4",alpha=0.3)
            rgl.points(mydatH[,1:3],col="Gray", alpha=0.2)
            rgl.points(mainVclust[,1:3],col="red", alpha=0.2)
            title3d(main = i, col="White")
            clids[cvpoints]=i
            cvpoints=cvpoints+1
          }
        }
      }
      Hclustchoice=dlgList(c("Unrefine",clids,"Refine","None"), multiple=FALSE, title="Select Homer1 cluster")
      Mainclust=Hclustchoice$res
      #if(Mainclust=="None") next()
      if(Mainclust=="Unrefine") simfact=simfact-0.4
      if(Mainclust!="Refine"& Mainclust!="Unrefine") refine=0 
      if(Mainclust=="Refine") simfact=simfact+0.2
      dlist=rgl.dev.list()
      dlistl=length(dlist)
      if(dlistl>0){
        for(i in 1:dlistl){
          rgl.close()
        }
      }
    }
    dlist=rgl.dev.list()
    dlistl=length(dlist)
    if(dlistl>0){
      for(i in 1:dlistl){
        rgl.close()
      }
    }
    if(Mainclust!="None"){
      Mainclust=as.numeric(Hclustchoice$res)
      mainHclust=concatH[which(concatH$Cluster==Mainclust),]
      mainHclustx=mainHclust
      mainVclustx=mainVclust
      #Homer1 vesicle cluster determination
      # mainHclust=as.matrix(mainHclust)
      mainHclust=mainHclust[,1:3]
      mainVclust=mainVclust[,1:3]
      bluech=nrow(mainVclust)
      redch=nrow(mainHclust)
      channelinfo[1:bluech]=1
      channelinfo[(bluech+1):(bluech+redch)]=2
      Resclusters=rbind(mainVclust,mainHclust)
      Resclusters=cbind(Resclusters,channelinfo)
      Respath=paste0("Clustertable",filen,bfoldername,".xls")
      Respath=paste(Clusterdir,Respath, sep = "/")
      write.delim(Resclusters,Respath, sep ="\t")
    }
    
  }
  Btable[blindn,4]=filen
  write.delim(Btable,Btablefile,sep ="\t")
  # filen=filen+1
} ##end file loop
# ROIsumm=as.data.frame(ROIsumm)
# names=c("ROI","Bassoon Volume", "Bsn Loc no.","Bsn nc no.", "median Bsn nc diameter", "median Bsn ncloc.no./total","Bsn extra nc nnd","median Bsn extrannd/ncnnd", "Homer1 Volume", "Homer1 cluster no.", "median Hmr nc diameter", "median Hmr ncloc.no./total", "Hmr extra nc nnd","median Hmr extrannd/ncnnd", "nanocolumn no."," no. SVs in bbox", "% Bsn locs in bbox","% Bsn nc locs in bbox", "Bsn-Hmr centroid dist.")
# setnames(ROIsumm,names)
# summpath=paste(ROIdir,iname,sep = "/")
# summpath=paste(summpath,".xls",sep = "q")
# write.delim(ROIsumm,summpath, sep ="\t")
print(filen-1)
