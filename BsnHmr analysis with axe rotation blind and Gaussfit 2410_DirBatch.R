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
library(bda)
library(xlsx)

#Variables
set.seed(23)
iname="ST003"
mpts=10
epsi=500
steepcut=0.0001
epcut=50
alpha=30
rad=124
totroi=0
skipp=0


#Functions
sqr=function(x){
  y=x*x
  return(y)
}
Rsqrgauss=function(x,mn,param){#param are parameters of gaussian function (mean, sd, and variance)
  valx=x[,2]-mn
  valfity=(param[3]*exp(-1/2*(x[,1]-param[1])^2/param[2]^2))-mn
  y=sum(sqr(valfity))/sum(sqr(valx))
  return(y)
}


Bdir="D:/Maureen/STORM/2colorBsnHmr/Blinded directories.xlsx"
BlindDirs=as.vector(read.xlsx(Bdir,1,header = FALSE))
blength=length(BlindDirs[,1])
Btabledir=as.character(BlindDirs[,])

pb=winProgressBar(title=paste0("Blind directory ",1,"of",blength),label =paste0("0% of 0"), min=0, max=100,initial=0)


for(bd in 1:blength) {
  #Method
  Btablefile=paste(Btabledir[bd],"Blindtable.xls",sep = "/")
  Btable=read.delim(Btablefile)
  Bno=nrow(Btable)
  
  exportdir=paste(Btabledir[bd],"Analysis", sep="/")
  
  dirrename=dir.exists(exportdir)
  
  while(dirrename==TRUE){
    exportdir=paste0(exportdir,"New")
    dirrename=dir.exists(exportdir)
  }
  dir.create(exportdir)
  
  Clusterdir=paste(exportdir,"ClusterInfo", sep="/")
  dir.create(Clusterdir)
  
  for(blindfold in 1:Bno){ #start ROI folder loop
    
    Objsumm=vector()
    Anasumm=vector()
    
    ROIdir=as.character(Btable[blindfold,1])
    ROIdir=paste0(Btabledir[bd],"/",ROIdir,"/",ROIdir,"Clusters")
    Clusterfiles=list.files(ROIdir, all.files = FALSE, full.names = TRUE, include.dirs = FALSE)
    nofull=length(grep("full.xls",Clusterfiles))
    if(nofull!=0) Clusterfiles=Clusterfiles[-grep("full.xls",Clusterfiles)]
    flength=length(Clusterfiles)
    
    expname=basename(as.character(Btable[blindfold,2]))
    clusterexp=paste0(Clusterdir,"/", expname,"Clusters")
    dir.create(clusterexp)
    
    for(filen in 1:flength) { ##start file loop
      # filen=1
      mydat=read.delim(Clusterfiles[filen],sep="\t")
      mydatV=as.matrix(mydat[which(mydat[,4]==1),1:3])
      mydatH=as.matrix(mydat[which(mydat[,4]==2),1:3])
      
      #Cross-talk substraction
      Hdist=knn.dist(mydatH, k=10)[,10]
      Bdist=knn.dist(mydatV,k=10)[,10]
      
      Hcut=median(Hdist)+mad(Hdist)
      Vcut=median(Bdist)+mad(Bdist)
      
      BtH=knnx.dist(mydatH[,1:3],mydatV[,1:3],k=10)[,10]
      HtB=knnx.dist(mydatV[,1:3],mydatH[,1:3],k=10)[,10]
      
      Hxt=HtB/Hdist
      Bxt=BtH/Bdist
      
      mydatH=mydatH[which(Hxt>1),]
      mydatV=mydatV[which(Bxt>1),]
      
      Hdb=dbscan(mydatH,Hcut,minPts = 10)
      maxref=max(Hdb$cluster)
      num=0
      for(i in 1:maxref){
        no=length(which(Hdb$cluster==i))
        stat=no-num
        if(stat>0){
          num=no
          clust=i
        }
      }
      mydatH=mydatH[which(Hdb$cluster==clust),]
      
      Vdb=dbscan(mydatV,Vcut,minPts = 10)
      maxref=max(Vdb$cluster)
      num=0
      for(i in 1:maxref){
        no=length(which(Vdb$cluster==i))
        stat=no-num
        if(stat>0){
          num=no
          clust=i
        }
      }
      mydatV=mydatV[which(Vdb$cluster==clust),]
      
      #Vamp2 main cluster info
      centv=c(mean(mydatV[,1]),mean(mydatV[,2]),mean(mydatV[,3]))
      aV=ashape3d(mydatV[,1:3],alpha=1)
      calpha=2*mean(aV$triang[,6])
      aV=ashape3d(mydatV[,1:3],alpha=calpha)
      VolmainV=volume_ashape3d(aV)
      
      #Bsn nc cluster determination
      Vn=nrow(mydatV)
      Vx=mydatV[,1]
      Vy=mydatV[,2]
      Vz=mydatV[,3]
      randx=runif(Vn*20, min = min(Vx), max= max(Vx) )
      randy=runif(Vn*20, min = min(Vy), max= max(Vy) )
      randz=runif(Vn*20, min = min(Vz), max= max(Vz) )
      randdist=cbind(randx,randy,randz)
      inaV=inashape3d(aV, indexAlpha = 1, randdist[,1:3])
      randdist=cbind(randdist,inaV)
      randdist=randdist[which(randdist[,4]==TRUE),]
      if(length(randdist[,1])<Vn) next 
      randdist=randdist[1:Vn,]
      oprand=optics(randdist, epsi, minPts = mpts)
      randreachmed=median(oprand$coredist[which(oprand$coredist!=Inf)])
      randreachmad=mad(oprand$coredist[which(oprand$coredist!=Inf)])
      vcutlow=randreachmed-(3*randreachmad)
      opvrandlow=extractDBSCAN(oprand, vcutlow)
      opv=optics(mydatV[,1:3], epsi, minPts = mpts)
      resvlow=extractDBSCAN(opv, vcutlow)
      mydatV=cbind(mydatV, resvlow$cluster, resvlow$order,resvlow$reachdist,resvlow$coredist)
      
      centroidsV=vector()
      ckeepv=vector()
      locnov=vector()
      cnndv=vector()
      diav=vector()
      vno=max(mydatV[,4])
      
      for(i in 1:vno){
        nclust=length(mydatV[which(mydatV[,4]==i),2])
        if(nclust>16) {
          diax=c(max(mydatV[which(mydatV[,4]==i),1])-min(mydatV[which(mydatV[,4]==i),1]))
          diay=c(max(mydatV[which(mydatV[,4]==i),2])-min(mydatV[which(mydatV[,4]==i),2]))
          cnnd=knn.dist(mydatV[which(mydatV[,4]==i),1:3],k=1)
          cnndv=c(cnndv,median(cnnd))
          locnov=c(locnov,nclust)
          cV=c(i,mean(mydatV[which(mydatV[,4]==i),1]),mean(mydatV[which(mydatV[,4]==i),2]),mean(mydatV[which(mydatV[,4]==i),3]))
          diav=c(diav,median(diax,diay))
          centroidsV=rbind(centroidsV,cV)
          ckeepv=c(ckeepv,i)
        } 
      }
      
      c2keepv=mydatV[,4]%in%ckeepv
      ConcatVcout=mydatV[-which(c2keepv==TRUE),1:3]
      outcno=nrow(ConcatVcout)
      outcnnd=knnx.dist(mydatV[,1:3],ConcatVcout,k=2)
      
      outcnnd=median(outcnnd[,2])
      locrat=locnov/Vn
      cnndrat=outcnnd/cnndv
      
      vclustno=length(ckeepv)
      
      
      
      #Homer1 main Cluster info
      centH=c(mean(mydatH[,1]),mean(mydatH[,2]),mean(mydatH[,3]))
      aH=ashape3d(mydatH[,1:3],alpha=1)
      calpha=2*mean(aH$triang[,6])
      aH=ashape3d(mydatH[,1:3],alpha=calpha)
      VolH=volume_ashape3d(aH)
      
      #Homer1 vesicle cluster determination
      
      Hn=nrow(mydatH)
      Hx=mydatH[,1]
      Hy=mydatH[,2]
      Hz=mydatH[,3]
      randx=runif(Hn*20, min = min(Hx), max= max(Hx) )
      randy=runif(Hn*20, min = min(Hy), max= max(Hy) )
      randz=runif(Hn*20, min = min(Hz), max= max(Hz) )
      randdist=cbind(randx,randy,randz)
      inaH=inashape3d(aH, indexAlpha = 1, randdist[,1:3])
      randdist=cbind(randdist,inaH)
      randdist=randdist[which(randdist[,4]==TRUE),]
      randdist=randdist[1:Hn,]
      oprand=optics(randdist, epsi, minPts = mpts)
      randreachmed=median(oprand$reachdist[-which(oprand$reachdist==Inf)])
      randreachmad=mad(oprand$reachdist[-which(oprand$reachdist==Inf)])
      Hcutlow=randreachmed-(3*randreachmad)
      opH=optics(mydatH[,1:3], epsi, minPts = mpts)
      resHlow=extractDBSCAN(opH, Hcutlow)
      mydatH=cbind(mydatH[,1:3], resHlow$cluster, resHlow$order,resHlow$reachdist,resHlow$coredist)
      Halphap=aH$triang[which(aH$triang[,9]==2),]
      Halphapt=t(Halphap[,1:3])
      
      centroidsVv=vector()
      locnoh=vector()
      cnndh=vector()
      ckeep=vector()
      diaH=vector()
      
      hclustno=unique(mydatH[,4])
      hcno=max(mydatH[,4])
      
      clusthno=0
      if(hcno>0){
        for(i in 1:hcno){
          nclusth=length(mydatH[which(mydatH[,4]==i),1])
          if(nclusth>16){
            diax=c(max(c(mydatH[which(mydatH[,4]==i),1]))-min(c(mydatH[which(mydatH[,4]==i),1])))
            diay=c(max(c(mydatH[which(mydatH[,4]==i),2]))-min(c(mydatH[which(mydatH[,4]==i),2])))
            diaH=c(diaH,median(diax,diay))
            cnnd=knn.dist(mydatH[which(mydatH[,4]==i),1:3],k=1)
            cnndh=c(cnndh,median(cnnd))
            locnoh=c(locnoh,nclusth)
            chinfo=c(paste0("Homer",i),median(mydatH[which(mydatH[,4]==i),1]),median(mydatH[which(mydatH[,4]==i),2]),median(mydatH[which(mydatH[,4]==i),3]))
            centroidsVv=rbind(centroidsVv,chinfo)
            clusthno=clusthno+1
            ckeep=c(ckeep,i)
          }
        }
        
        c2keep=mydatH[,4]%in%ckeep
        ConcatHcout=mydatH[-which(c2keep==TRUE),1:3]
        outcnoh=nrow(ConcatHcout)
        outcnndh=knnx.dist(mydatH[,1:3],ConcatHcout,k=2)
        outcnndh=median(outcnndh[,2])
        
        locnohrat=locnoh/Hn
        cnndhrat=outcnndh/cnndh
      }  
      
      #ROI Bsn and Homer1 nanocluster info
      bno=length(ckeepv)
      hcno=length(ckeep)
      if(hcno>0){
        hna=vector()
        hna[hcno]=NA
        centroidsVv=cbind(centroidsVv,locnoh,diaH,locnohrat,cnndh,cnndhrat,hna,hna)
      } else {
        centroidsVv=c("Homer", centH,NA,NA,NA,NA,NA,NA,NA)
      }
      if(bno>0){
        bna=vector()
        bna[bno]=NA
        centroidsV=cbind(centroidsV,locnov,diav,locrat,cnndv,cnndrat,bna,bna)
      } else {
        centroidsV=c("Bassoon", centH,NA,NA,NA,NA,NA,NA,NA)
      }
      
     
      
      centroidsVr=rbind(centroidsV,centroidsVv)
      centroidsVr=as.data.table(centroidsVr)
      names=c("Cluster","x [nm]", "y [nm]","z [nm]", "Loc. no.", "NC diameter [nm]","Loc NC/out ratio","NC nnd [nm]", "NC/out nnd ratio", "SV nnd [nm]", "SV dist. to Hmr node [nm]")
      setnames(centroidsVr,names)
      centroidsVtable=paste0(clusterexp,"/","Info",basename(Clusterfiles[filen]))
      write.delim(centroidsVr,centroidsVtable, sep ="/t")
      
      #Homer1 ellipse
      Hellihull=ellipsoidhull(mydatH[,1:3],tol = 1000)
      Haxes=ellipse3d.axes(Hellihull$cov, centre = Hellihull$loc)
      rgl.close()
      
      #translation to origin
      mydatVx=rbind(cbind(mydatV[,1]-Hellihull$loc[1],mydatV[,2]-Hellihull$loc[2],mydatV[,3]-Hellihull$loc[3]),c(centv[1]-Hellihull$loc[1],centv[2]-Hellihull$loc[2],centv[3]-Hellihull$loc[3]))
      mydatHx=cbind(mydatH[,1]-Hellihull$loc[1],mydatH[,2]-Hellihull$loc[2],mydatH[,3]-Hellihull$loc[3])
      Haxesx=cbind(Haxes[,1]-Hellihull$loc[1],Haxes[,2]-Hellihull$loc[2],Haxes[,3]-Hellihull$loc[3])
      
      ##rotation around z axis (xycoordinat rotation)
      vAC=as.vector(0-Haxesx[6,1:2])
      vCO=as.vector(c(1,0))
      if(Haxesx[6,2]<0){
        angle=2*pi-acos(dot(vAC,vCO)/(knn.dist(rbind(c(0,0),vAC),k=1)[1]*knn.dist(rbind(c(0,0),vCO),k=1)[1]))
      } else {
        angle=acos(dot(vAC,vCO)/(knn.dist(rbind(c(0,0),vAC),k=1)[1]*knn.dist(rbind(c(0,0),vCO),k=1)[1]))
      }
      
      mydatVx1=mydatVx
      mydatVx1[,1]=mydatVx[,1]*cos(angle)-mydatVx[,2]*sin(angle)
      mydatVx1[,2]=mydatVx[,2]*cos(angle)+mydatVx[,1]*sin(angle)
      
      mydatHx1=mydatHx
      mydatHx1[,1]=mydatHx[,1]*cos(angle)-mydatHx[,2]*sin(angle)
      mydatHx1[,2]=mydatHx[,2]*cos(angle)+mydatHx[,1]*sin(angle)
      
      Haxesx1=Haxesx
      Haxesx1[,1]=Haxesx[,1]*cos(angle)-Haxesx[,2]*sin(angle)
      Haxesx1[,2]=Haxesx[,2]*cos(angle)+Haxesx[,1]*sin(angle)
      
      mydatVx=mydatVx1
      mydatHx=mydatHx1
      Haxesx=Haxesx1
      
      ##rotation around y axis (xz coordinate rotation)
      vAC=as.vector(0-c(Haxesx[6,1],Haxesx[6,3]))
      vCO=as.vector(c(1,0))
      if(Haxesx[6,3]<0){
        angle=2*pi-acos(dot(vAC,vCO)/(knn.dist(rbind(c(0,0),vAC),k=1)[1]*knn.dist(rbind(c(0,0),vCO),k=1)[1]))
      } else {
        angle=acos(dot(vAC,vCO)/(knn.dist(rbind(c(0,0),vAC),k=1)[1]*knn.dist(rbind(c(0,0),vCO),k=1)[1]))
      }
      mydatVx1=mydatVx
      mydatVx1[,1]=(mydatVx[,1]*cos(angle)-mydatVx[,3]*sin(angle))
      mydatVx1[,3]=(mydatVx[,3]*cos(angle)+mydatVx[,1]*sin(angle))
      
      mydatHx1=mydatHx
      mydatHx1[,1]=mydatHx[,1]*cos(angle)-mydatHx[,3]*sin(angle)
      mydatHx1[,3]=mydatHx[,3]*cos(angle)+mydatHx[,1]*sin(angle)
      
      Haxesx1=Haxesx
      Haxesx1[,1]=Haxesx[,1]*cos(angle)-Haxesx[,3]*sin(angle)
      Haxesx1[,3]=Haxesx[,3]*cos(angle)+Haxesx[,1]*sin(angle)
      
      mydatVx=mydatVx1
      mydatHx=mydatHx1
      Haxesx=Haxesx1
      
      ##rotation around x axis (yz coordinate rotation)
      vAC=as.vector(0-c(Haxesx[2,2],Haxesx[2,3]))
      vCO=as.vector(c(1,0))
      
      if(Haxesx[2,3]<0){
        angle=2*pi-acos(dot(vAC,vCO)/(knn.dist(rbind(c(0,0),vAC),k=1)[1]*knn.dist(rbind(c(0,0),vCO),k=1)[1]))
      } else {
        angle=acos(dot(vAC,vCO)/(knn.dist(rbind(c(0,0),vAC),k=1)[1]*knn.dist(rbind(c(0,0),vCO),k=1)[1]))
      }
      
      mydatVx1=mydatVx
      mydatVx1[,2]=(mydatVx[,2]*cos(angle)-mydatVx[,3]*sin(angle))
      mydatVx1[,3]=(mydatVx[,3]*cos(angle)+mydatVx[,2]*sin(angle))
      
      mydatHx1=mydatHx
      mydatHx1[,2]=mydatHx[,2]*cos(angle)-mydatHx[,3]*sin(angle)
      mydatHx1[,3]=mydatHx[,3]*cos(angle)+mydatHx[,2]*sin(angle)
      
      Haxesx1=Haxesx
      Haxesx1[,2]=Haxesx[,2]*cos(angle)-Haxesx[,3]*sin(angle)
      Haxesx1[,3]=Haxesx[,3]*cos(angle)+Haxesx[,2]*sin(angle)
      
      mydatVx=mydatVx1
      mydatHx=mydatHx1
      Haxesx=Haxesx1
      
      xval=abs(mydatVx[length(mydatVx[,1]),1])
      centvrot=mydatVx[length(mydatVx[,1]),]
      mydatVx=mydatVx[1:length(mydatVx[,1])-1,]
      
      #loc chunk for gaussian analysis
      datchunkV=mydatVx[which(-100<mydatVx[,2]&mydatVx[,2]<100),1]
      if(length(datchunkV)<2)next()
      datchunkH=mydatHx[which(-100<mydatHx[,2]&mydatHx[,2]<100),1]
      if(length(datchunkH)<2)next()
      chunkinfo=rbind(cbind(datchunkV,1),cbind(datchunkH,2))
      
      #Gaussian fitting
      Hbin=binning(datchunkH[order(datchunkH)],bw=5)
      Vbin=binning(datchunkV[order(datchunkV)],bw=5)
      
      #hmr peak fit
      x <- Hbin$mids
      f <- function(par)
      {
        m <- par[1]
        sd <- par[2]
        k <- par[3]
        rhat <- k * exp(-0.5 * ((x - m)/sd)^2)
        sum((Hbin$counts - rhat)^2)
      }
      
      resfith=optim(c(mean(datchunkH), 20, 10), f, method="BFGS", control=list(reltol=1e-9))
      v=resfith$par
      Rhmr=Rsqrgauss(cbind(Hbin$mids,Hbin$counts),mean(Hbin$counts),v)
      if(Rhmr>1)Rhmr=2-Rhmr
      
      if(Rhmr<0.75) next()
      
      #bsn peak fit
      x <- Vbin$mids
      f <- function(par)
      {
        m <- par[1]
        sd <- par[2]
        k <- par[3]
        rhat <- k * exp(-0.5 * ((x - m)/sd)^2)
        sum((Vbin$counts - rhat)^2)
      }
      resfitb=optim(c(mean(datchunkV), 20, 10), f, method="BFGS", control=list(reltol=1e-9))
      v=resfitb$par
      Rbsn=Rsqrgauss(cbind(Vbin$mids,Vbin$counts),mean(Vbin$counts),v)
      if(Rbsn>1)Rbsn=2-Rbsn
      
      clusteriexp=paste0(clusterexp,"/", "gausslocs",filen,".xls")
      write_delim(as.data.frame(chunkinfo),clusteriexp,delim = "\t")
      
      ptop=abs(resfitb$par[1]-resfith$par[1])
      plot(x=c(Vbin$mids,Hbin$mids),y=c(Vbin$counts,Hbin$counts),main = ptop)
      newhlen=max(mydatHx[,2])-min(mydatHx[,2])
      newhdepth=4.71*resfith$par[2]
      newhwidth=max(mydatHx[,3])-min(mydatHx[,3])
      
      newblen=max(mydatVx[,2])-min(mydatVx[,2])
      newbdepth=4.71*resfitb$par[2]
      newbwidth=max(mydatVx[,3])-min(mydatVx[,3])
      
      #Homer1 alpha translation along ellipse axes
      
      #Homer1 bbox projection
      if(centvrot[1]>0) rad=124 else rad=-124
      EBsn=rbind(c(rad,max(mydatHx[,2]),0),c(rad,min(mydatHx[,2]),0),c(rad,0,max(mydatHx[,3])),c(rad,0,min(mydatHx[,3])),c(rad+(newhdepth/2)+20,0,0),c(rad-(newhdepth/2)-20,0,0))
      plengthax=dist(EBsn[1:2,])*0.05
      pwidthax=dist(EBsn[3:4,])*0.05
      pdepthax=dist(EBsn[5:6,])*0.05
      EBsn=rbind(c(rad,max(mydatHx[,2])+plengthax,0),c(rad,min(mydatHx[,2])-plengthax,0),c(rad,0,max(mydatHx[,3])+pwidthax),c(rad,0,min(mydatHx[,3])-pwidthax),c(rad+(newhdepth/2)+20,0,0),c(rad-(newhdepth/2)-20,0,0))
      
      #Bounding box
      centaxes=c(rad, 0, 0)
      AxeX=EBsn[1:2,]
      projtop=EBsn[3:6,]
      Boxpoints=vector()
      Boxproj=vector()
      interpoints=vector()
      
      #AxeX point projection to Axes plane 
      for(i in 1:2) {
        Norm=c((centaxes[1]-projtop[i,1]),(centaxes[2]-projtop[i,2]),(centaxes[3]-projtop[i,3]))
        dplane=Norm[1]*projtop[i,1] + Norm[2]*projtop[i,2] + Norm[3]*projtop[i,3]
        plane=c(Norm, dplane)
        t=as.numeric((dplane-(AxeX[1,1]*Norm[1])-(AxeX[1,2]*Norm[2])-(AxeX[1,3]*Norm[3]))/((Norm[1]^2)+(Norm[2]^2)+(Norm[3]^2)))
        bpoint1=c((AxeX[1,1]+t*Norm[1]),(AxeX[1,2]+t*Norm[2]),(AxeX[1,3]+t*Norm[3]))
        bpoint2=c((AxeX[2,1]+t*Norm[1]),(AxeX[2,2]+t*Norm[2]),(AxeX[2,3]+t*Norm[3]))
        Boxproj=rbind(Boxproj,bpoint1,bpoint2)
      }
      Boxpoints=Boxproj
      
      #Bbox area
      lengthax=as.numeric(dist(EBsn[1:2,]))/1000
      widthax=as.numeric(dist(EBsn[3:4,]))/1000
      depthax=as.numeric(dist(EBsn[5:6,]))/1000
      hmrarea=lengthax*widthax
      bboxvol=lengthax*widthax*depthax
      
      
      for(i in 3:4){
        bpoints=Boxproj
        Norm=c((centaxes[1]-projtop[i,1]),(centaxes[2]-projtop[i,2]),(centaxes[3]-projtop[i,3]))
        dplane=Norm[1]*projtop[i,1] + Norm[2]*projtop[i,2] + Norm[3]*projtop[i,3]
        plane=c(Norm, dplane)
        t=as.numeric((dplane-(AxeX[1,1]*Norm[1])-(AxeX[1,2]*Norm[2])-(AxeX[1,3]*Norm[3]))/((Norm[1]^2)+(Norm[2]^2)+(Norm[3]^2)))
        trans=c(t*Norm[1],t*Norm[2],t*Norm[3])
        bpoints[,1]=Boxproj[,1]+trans[1]
        bpoints[,2]=Boxproj[,2]+trans[2]   
        bpoints[,3]=Boxproj[,3]+trans[3] 
        interpoints=rbind(interpoints,bpoints)
        Boxpoints=rbind(Boxpoints,bpoints)
      }
      # Boxpoints=rbind(Boxpoints,EBsn)
      boxa=ashape3d(Boxpoints,alpha=10000, pert = TRUE)
      btri=boxa$triang
      btrit=t(btri[which(btri[,9]==2),1:3])
      
      #Test of Bsn locs within bbox
      cinBa=inashape3d(boxa,indexAlpha = 1, mydatVx[,1:3])
      cinBano=(length(cinBa)-length(cinBa[which(cinBa==FALSE)]))*100/length(cinBa)
      
      open3d()
      rgl.points(mydatVx[,1:3],col="red",alpha=0.4,size=5)
      rgl.points(mydatHx[,1:3],col="chartreuse4",alpha=0.4,size=5)
      rgl.triangles(Boxpoints[btrit,],col="black",alpha=0.3)
      segments3d(Haxesx)
      
      #Summary of Homer and Vamp2 object metrics
      Objdat=c(basename(Clusterfiles[filen]),Rhmr,Rbsn,cinBano,ptop,knn.dist(rbind(centH,centv),k=1)[1],xval, VolmainV, Vn, VolH,Hn,lengthax, widthax, depthax,bboxvol,newhlen,newhwidth,newhdepth,newblen,newbwidth,newbdepth)
      if(bno==0)   Objdat=c(Objdat,0, NA,outcnnd, NA,NA) else Objdat=c(Objdat,bno, median(diav),median(locrat),outcnnd, median(cnndrat)) 
      if(hcno==0)  Objdat=c(Objdat,0, NA,outcnndh, NA,NA) else Objdat=c(Objdat,hcno, median(diaH),median(locnohrat),outcnndh, median(cnndhrat))
      
      Objsumm=rbind(Objsumm,Objdat)
    
      ijstat=trunc(filen*100/flength)
      setWinProgressBar(pb,ijstat, title=paste0("Blind directory ",bd,"of",blength),label = paste0(ijstat,"% of ",blindfold,"/",Bno))
      
    } ##end file loop
    
    Objsumm=as.data.table(Objsumm)
    names=c("ROI","Hmr fit R2","Bsn fit R2","Bsn locs in bbox", "peak to peak distance","centroid distances","centroid xval dist", "Bassoon Volume", "Bassoon Loc no.","Homer1 Volume","Hmr Loc. no.","Homer bbox length (nm)", "Homer bbox width (nm)","Homer Bbox depth (nm)", "Homer bbox volume (µm3)","newhlen","newhwidth","newhdepth","newblen","newbwidth","newbdepth","Bsn nc no.","median Bnc diameter","median Bncloc.no./total","extra Bnc nnd","median Bnc extrannd/ncnnd", "Homer1 cluster no.", "median Hmr nc diameter", "median Hmr ncloc.no./total", "Hmr extra nc nnd","median Hmr extrannd/ncnnd")
    setnames(Objsumm,names)
    objexp=paste0(exportdir,"/","Objects",basename(as.character(Btable[blindfold,2])),".xls")
    write.delim(Objsumm,objexp, sep ="\t")
  }#end blind folder loop
}#end Blind directory loop

cpb=close(pb)
print("Done")
