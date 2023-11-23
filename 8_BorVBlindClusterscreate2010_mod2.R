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
mpts=100 #initial value of the minimum points paramater in dbscan
epsi=500 #initial value of the radius paramater in dbscan
# steepcut=0.0001
epcut=50
alpha=30
rad=125
totroi=0
skipp=0

xy_edge=200
z_edge=50
COM_dlim=300

#Functions
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

COM_apply=function(xyz_data,cluster_list){
  max_clust=max(cluster_list)
  COM_table=matrix(ncol = 4, nrow = max_clust)
  for(i in 1:max_clust){
    COM_table[i,]=c(median(xyz_data[which(cluster_list==i),1]),median(xyz_data[which(cluster_list==i),2]),median(xyz_data[which(cluster_list==i),3]),i)
  }
  return(COM_table)
}

rgl.clusters=function(dataset, cluster_vector){
  max_clust=max(cluster_vector)
  rain=rainbow(max_clust)
  for(i in 1:max_clust){
    rgl.points(dataset[which(cluster_vector==i),1:3], color=rain[i], alpha=0.3)
  }
  
}

dist_3d=function(dataset,query){
  dist_vec=matrix(nrow=nrow(query),ncol=nrow(dataset))
  for(i in 1:nrow(query)){
    dist_vec[i,]=apply(dataset,MARGIN = 1,function(x){((x[1]-query[i,1])^2+(x[2]-query[i,2])^2+(x[3]-query[i,3])^2)^(1/2)})
    dist_vec[i,]=sort(dist_vec[i,])
  }
  return(dist_vec)
}


#Method
#ROIdir=dlgDir(default= "C:/Users/Maureen/Documents/STORM", "Select ROI directory")

Btabledir="E:/MM-Lenkei/STORM-Recs/2 color Bassoon-Homer1/SD414/Recs_Final/BH_SD414_Veh2_Recs/Blinded"
Btablefile=paste(Btabledir,"Blindtable.xls",sep = "/")
Btable=read.delim(Btablefile)
Bno=nrow(Btable)

Clusterdone=0
clusttype=select.list(c("Vamp2","Bassoon"),title = "Select the type of staining to cluster:",graphics = TRUE)

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
  mydatVo=mydatVo[which(mydatVo$uncertainty_xy..nm.<50 & mydatVo$uncertainty_z..nm.< 200),]
  nameh=grep(modnameV[filen], modnameH, value=FALSE)
  mydatHo=read.csv(FileCy3[nameh])
  mydatHo=mydatHo[which(mydatHo$uncertainty_xy..nm.<50 & mydatHo$uncertainty_z..nm.< 200),]
  mydatV=mydatVo[which(mydatVo$uncertainty_xy..nm.<50 & mydatVo$uncertainty_z..nm.< 200),2:4]
  mydatH=mydatHo[which(mydatHo$uncertainty_xy..nm.<50 & mydatHo$uncertainty_z..nm.< 200),2:4]
  
  open3d()
  rgl.bg(color="black")
  rgl.points(mydatV[,1:3],col="Red", alpha=0.3)
  rgl.points(mydatH[,1:3],col="Green", alpha=0.3)
  
  ##Vamp2 Cluster determination
  op=optics(mydatV, epsi, minPts = mpts)
  
  ###Uniform distances_405
  Vn=nrow(mydatV)
  
  VolumeV=(max(mydatV[,1])-min(mydatV[,1]))*(max(mydatV[,2])-min(mydatV[,2]))*(max(mydatV[,3])-min(mydatV[,3]))
  ptVol_V=VolumeV/Vn
  Uni_V=2*(3*ptVol_V/(4*3.1415926))^(1/3)
  
  ###Uniform distance_cy3
  Hn=nrow(mydatH)
  
  VolumeH=(max(mydatH[,1])-min(mydatH[,1]))*(max(mydatH[,2])-min(mydatH[,2]))*(max(mydatH[,3])-min(mydatH[,3]))
  ptVol_H=VolumeH/Hn
  Uni_H=2*(3*ptVol_H/(4*3.1415926))^(1/3)
  
  
  ##405dye preliminary cluster
  simfact=4
  refine=1
  
  while(refine==1){
    maincut=Uni_V-Uni_V^(1/2)*simfact
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
    
    # for(i in 1:maxCV){
    #   nr=length(concatVtemp[which(concatVtemp[,4]==i),4])
    #   if(nr>50){
    #     aV=ashape3d(concatVtemp[which(concatVtemp[,4]==i),1:3],alpha=1)
    #     calpha=2*mean(aV$triang[,6])
    #     aV=ashape3d(concatVtemp[which(concatVtemp[,4]==i),1:3],alpha=calpha)
    #     VolV=volume_ashape3d(aV)
    #     #concatVtemp=concatVtemp[which(concatVtemp[,4]!=i),]
    #     if(VolV>vollim){
    #       open3d()
    #       rgl.bg(color="black")
    #       rgl.points(concatVtemp[which(concatVtemp[,4]==i),1:3], col="coral",alpha=0.3)
    #       rgl.points(mydatV[,1:3],col="Gray", alpha=0.2)
    #       rgl.points(mydatH[,1:3],col="White", alpha=0.2)
    #       title3d(main = i, col="White")
    #       centv=c(mean(concatV[which(concatV$Cluster==i),1]),mean(concatV[which(concatV$Cluster==i),2]),mean(concatV[which(concatV$Cluster==i),3]), i, VolV)
    #       centVm=rbind(centVm,centv)
    #       clids[cvpoints]=i
    #       cvpoints=cvpoints+1
    #     }
    #   }
    # }
    open3d()
    rgl.bg(color="black")
    #rgl.points(concatVtemp[,1:3], col="darkorchid2",alpha=0.3)
    rgl.clusters(concatVtemp[,1:3], concatVtemp[,4])
    rgl.points(mydatV[,1:3],col="Gray", alpha=0.3)
    rgl.points(mydatH[,1:3],col="White", alpha=0.3)
    
    Vclustchoice=select.list(c("Unrefine","***Perfect***","Refine","None"), multiple=FALSE, title="Select 405 cluster",graphics = TRUE)
    Mainclust=Vclustchoice
    if(Mainclust=="Unrefine") simfact=simfact-1
    if(Mainclust!="Refine"& Mainclust!="Unrefine") refine=0 
    if(Mainclust=="Refine") simfact=simfact+0.5
    dlist=rgl.dev.list()
    dlistl=length(dlist)
    if(dlistl>0){
      for(i in 1:dlistl){
        rgl.close()
      }
    }
  }
 
  if(Mainclust!="None"){
    # Mainclust=as.numeric(Vclustchoice)
    # mainVclust=concatV[which(concatV$Cluster==Mainclust),]
    
    ##Homer1 main clusters
    op=optics(mydatH, epsi, minPts = mpts)
    
    ##Cy3dye preliminary cluster
    simfact=5
    refine=1
    
    while(refine==1){
      maincut=Uni_H-Uni_H^(1/2)*simfact
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
      # 
      # for(i in 1:maxCH){
      #   nr=length(concatHtemp[which(concatHtemp[,4]==i),4])
      #   if(nr>50){
      #     aH=ashape3d(concatHtemp[which(concatHtemp[,4]==i),1:3],alpha=1)
      #     calpha=2*mean(aH$triang[,6])
      #     aH=ashape3d(concatHtemp[which(concatHtemp[,4]==i),1:3],alpha=calpha)
      #     VolH=volume_ashape3d(aH)
      #     #concatHtemp=concatHtemp[which(concatHtemp[,4]!=i),]
      #     if(VolH>2000000){
      #       open3d()
      #       rgl.bg(color="black")
      #       rgl.points(concatHtemp[which(concatHtemp[,4]==i),1:3], col="chartreuse4",alpha=0.3)
      #       rgl.points(mydatH[,1:3],col="Gray", alpha=0.2)
      #       rgl.points(mainVclust[,1:3],col="red", alpha=0.2)
      #       title3d(main = i, col="White")
      #       clids[cvpoints]=i
      #       cvpoints=cvpoints+1
      #     }
      #   }
      # }
      open3d()
      rgl.bg(color="black")
      rgl.clusters(concatHtemp[,1:3], concatHtemp[,4])
      rgl.points(mydatH[,1:3],col="Gray", alpha=0.2)
      rgl.points(mydatV[,1:3],col="White", alpha=0.6)
      
      
      Hclustchoice=select.list(c("Unrefine","***Perfect***","Refine","None"), multiple=FALSE, title="Select Homer1 cluster",graphics = TRUE)
      Mainclust=Hclustchoice
      #if(Mainclust=="None") next()
      if(Mainclust=="Unrefine") simfact=simfact-1
      if(Mainclust!="Refine"& Mainclust!="Unrefine") refine=0 
      if(Mainclust=="Refine") simfact=simfact+0.5
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
      # Mainclust=as.numeric(Hclustchoice)
      # mainHclust=concatH[which(concatH$Cluster==Mainclust),]
      COM_V=COM_apply(concatV[,1:3],concatV[,4])
      COM_H=COM_apply(concatH[,1:3],concatH[,4])
      Vollimits=c(min(c(min(mydatV[,1]),min(mydatH[,1]))),max(c(max(mydatV[,1]),max(mydatH[,1]))),
                  min(c(min(mydatV[,2]),min(mydatH[,2]))),max(c(max(mydatV[,2]),max(mydatH[,2]))),
                  min(c(min(mydatV[,3]),min(mydatH[,3]))),max(c(max(mydatV[,3]),max(mydatH[,3])))
                  )
      COM_V=COM_V[which(COM_V[,1]>(Vollimits[1]+xy_edge) & COM_V[,1]<(Vollimits[2]-xy_edge) &
                        COM_V[,2]>(Vollimits[3]+xy_edge) & COM_V[,2]<(Vollimits[4]-xy_edge) &  
                        COM_V[,3]>(Vollimits[5]+z_edge) & COM_V[,3]<(Vollimits[6]-z_edge)),]
      
      COM_H=COM_H[which(COM_H[,1]>(Vollimits[1]+xy_edge) & COM_H[,1]<(Vollimits[2]-xy_edge) &
                          COM_H[,2]>(Vollimits[3]+xy_edge) & COM_H[,2]<(Vollimits[4]-xy_edge) &  
                          COM_H[,3]>(Vollimits[5]+z_edge) & COM_H[,3]<(Vollimits[6]-z_edge)),]
    
      if(length(COM_H)==0 | length(COM_V)==0 ) next
      
      if(is.matrix(COM_V)==FALSE) COM_V=cbind(COM_V[1],COM_V[2],COM_V[3],COM_V[4])
      
      if(length(COM_H)>length(COM_V)){
        if((length(COM_H)/4)>2){
          COM_knn=kNN(COM_H, 2, query=COM_V)
         
          Vclusters=which(COM_knn$dist[,1]<COM_dlim)
          Hclusters=COM_knn$id[Vclusters,1]
        }else{
          COM_knn=dist_3d(COM_H[,1:3],COM_V[,1:3])
          COM_knn[which(COM_knn[,2]>COM_dlim),2]='NA'
          Hclusters=which(COM_knn[,1]<COM_dlim)
          Vclusters=COM_knn$id[Hclusters,1]
        }
      }else if((length(COM_V)/4)>2){
        COM_knn=kNN(COM_V, 2, query=COM_H)
        Hclusters=which(COM_knn$dist[,1]<COM_dlim)
        Vclusters=COM_knn$id[Hclusters,1]
      }else{
        COM_knn=dist_3d(COM_H[,1:3],COM_V[,1:3])
      }
      }
       

# else if(length(COM_H[,1])<=1 & length(COM_V[,1])>2){
#        COM_knn=kNN(COM_V,2, query=COM_H)
#        Hclusters=which(COM_knn$dist[,1]<300)
#        Vclusters=COM_knn$id[Hclusters,1]
#      } else if (length(COM_H[,1])==1 & length(COM_V[,1])== 1){
#        COM_knn=((COM_H[1]-COM_V[1])^2+(COM_H[2]-COM_V[2])^2+(COM_H[3]-COM_V[3])^2)^(1/2)
#      }
       
     
      
      open3d()
      rgl.bg(color = "black")
      rgl.points(concatHtemp[which(concatHtemp[,4]%in%Hclusters),1:3],color="green")
      rgl.points(concatVtemp[which(concatVtemp[,4]%in%Vclusters),1:3],color="red")
      rgl.points(mydatV[,1:3],color="gray", alpha=0.2)
      # 
      # #Homer1 vesicle cluster determination
      # # mainHclust=as.matrix(mainHclust) 
      # mainHclust=mainHclust[,1:3]
      # mainVclust=mainVclust[,1:3]
      # bluech=nrow(mainVclust)
      # redch=nrow(mainHclust)
      # channelinfo[1:bluech]=1
      # channelinfo[(bluech+1):(bluech+redch)]=2
      # Resclusters=rbind(mainVclust,mainHclust)
      # Resclusters=cbind(Resclusters,channelinfo)
      # Respath=paste0("Clustertable",filen,bfoldername,".xls")
      # Respath=paste(Clusterdir,Respath, sep = "/")
      # write.delim(Resclusters,Respath, sep ="\t")
    }
    
  }
  Btable[blindn,4]=filen
  write.delim(Btable,Btablefile,sep ="\t")
  # filen=filen+1
} ##end file loop
 

print(filen-1)
