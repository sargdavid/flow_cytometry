library(flowCore)
library(flowAI)
library(raster)
library(pracma)
library(MASS)

#flowset
myfiles<- list.files(path="/Users/md1518/Downloads/flow_cytometry/HEUvsUE/HEU/CONTROL/",pattern=".FCS$")
HEU_control<- flowCore::read.flowSet(myfiles,path="/Users/md1518/Downloads/flow_cytometry/HEUvsUE/HEU/CONTROL/")
spillover(HEU_control[[1]])
HEU_control_comp<- compensate(HEU_control,spillover(HEU_control[[1]])$SPILL)
trans2_HEU_control<-estimateLogicle(HEU_control_comp[[1]],colnames(HEU_control_comp[,3:10]))
HEU_control_comp_trans<- transform(HEU_control_comp,trans2_HEU_control)

HEU_control<- list()
for(i in 1:20){
  HEU_control[[i]] <- HEU_control_comp_trans[[i]]@exprs
}

#the column name does not work
#0 is control, 1 is stimulated
for(i in 1:20){
  HEU_control[[i]]<- cbind(HEU_control[[i]],0)
  names(HEU_control[[i]][,12])<- "Stimulation"
}
nu<- c(36,39,19,13,26,3,10,43,42,29,44,41,28,23,22,16,12,38,20,17)
for(i in 1:20){
  HEU_control[[i]]<- cbind(HEU_control[[i]],nu[i])
  names(HEU_control[[i]][,13])<- "Patient"
}
nuu<- c(22,29,36,43,50,57,71,106,113,162,169,176,190,218,225,232,239,260,274,281)
for(i in 1:20){
  HEU_control[[i]]<- cbind(HEU_control[[i]],nuu[i])
  names(HEU_control[[i]][,14])<- "File"
}

myfiles2<- list.files(path="/Users/md1518/Downloads/flow_cytometry/HEUvsUE/HEU/LPS/",pattern=".FCS$")
HEU_LPS<- flowCore::read.flowSet(myfiles2,path="/Users/md1518/Downloads/flow_cytometry/HEUvsUE/HEU/LPS/")
HEU_LPS_comp<- compensate(HEU_LPS,spillover(HEU_LPS[[1]])$SPILL)
trans2_HEU_LPS<-estimateLogicle(HEU_LPS_comp[[1]],colnames(HEU_LPS_comp[,3:10]))
HEU_LPS_comp_trans<- transform(HEU_LPS_comp,trans2_HEU_LPS)

HEU_LPS<- list()
for (i in 1:20){
  HEU_LPS[[i]] <- HEU_LPS_comp_trans[[i]]@exprs
}
#0 is control, 1 is stimulated
for(i in 1:20){
  HEU_LPS[[i]]<- cbind(HEU_LPS[[i]],1)
  names(HEU_LPS[[i]][,12])<- "Stimulation"
}
nu<- c(36,39,19,13,26,3,10,43,42,29,44,41,28,23,22,16,12,38,20,17)
for(i in 1:20){
  HEU_LPS[[i]]<- cbind(HEU_LPS[[i]],nu[i])
  names(HEU_LPS[[i]][,13])<- "Patient"
}
nuu<- c(24,31,38,45,52,59,73,108,115,164,171,178,192,220,227,234,241,262,276,283)
for(i in 1:20){
  HEU_LPS[[i]]<- cbind(HEU_LPS[[i]],nuu[i])
  names(HEU_LPS[[i]][,14])<- "File"
}

myfiles3<- list.files(path="/Users/md1518/Downloads/flow_cytometry/HEUvsUE/UE/CONTROL/",pattern=".FCS$")
UE_control<- flowCore::read.flowSet(myfiles3,path="/Users/md1518/Downloads/flow_cytometry/HEUvsUE/UE/CONTROL/")
UE_control_comp<- compensate(UE_control,spillover(UE_control[[1]])$SPILL)
trans2_UE_control<-estimateLogicle(UE_control_comp[[1]],colnames(UE_control_comp[,3:10]))
UE_control_comp_trans<- transform(UE_control_comp,trans2_UE_control)

UE_control<- list()
for (i in 1:24){
  UE_control[[i]]<- UE_control_comp_trans[[i]]@exprs
}
#0 is control, 1 is stimulated
for(i in 1:20){
  UE_control[[i]]<- cbind(UE_control[[i]],0)
  names(UE_control[[i]][,12])<- "Stimulation"
}
nu<- c(21,8,37,7,6,30,27,4,31,25,11,1,9,2,35,18,34,32,14,15,24,5,40,33)
for(i in 1:20){
  UE_control[[i]]<- cbind(UE_control[[i]],nu[i])
  names(UE_control[[i]][,13])<- "Patient"
}
nuu<- c(1,8,15,64,78,85,92,99,120,127,134,141,148,155,183,197,204,211,246,253,267,288,295,302)
for(i in 1:20){
  UE_control[[i]]<- cbind(UE_control[[i]],nuu[i])
  names(UE_control[[i]][,14])<- "File"
}

myfiles4<- list.files(path="/Users/md1518/Downloads/flow_cytometry/HEUvsUE/UE/LPS/",pattern=".FCS$")
UE_LPS<- flowCore::read.flowSet(myfiles4,path="/Users/md1518/Downloads/flow_cytometry/HEUvsUE/UE/LPS/")
UE_LPS_comp<- compensate(UE_LPS,spillover(UE_LPS[[1]])$SPILL)
trans2_UE_LPS<-estimateLogicle(UE_LPS_comp[[1]],colnames(UE_LPS_comp[,3:10]))
UE_LPS_comp_trans<- transform(UE_LPS_comp,trans2_UE_LPS)

UE_LPS_nontransformed<- list()
for (i in 1:24){
  UE_LPS_nontransformed[[i]] <- UE_LPS_comp[[i]]@exprs
}

UE_LPS<- list()
for (i in 1:24){
  UE_LPS[[i]] <- UE_LPS_comp_trans[[i]]@exprs
}
#0 is control, 1 is stimulated
for(i in 1:20){
  UE_LPS[[i]]<- cbind(UE_LPS[[i]],1)
  names(UE_LPS[[i]][,12])<- "Stimulation"
}
nu<- c(21,8,37,7,6,30,27,4,31,25,11,1,9,2,35,18,34,32,14,15,24,5,40,33)
for(i in 1:20){
  UE_LPS[[i]]<- cbind(UE_LPS[[i]],nu[i])
  names(UE_LPS[[i]][,13])<- "Patient"
}
nuu<- c(3,10,17,66,80,87,94,101,122,129,136,143,150,157,185,199,206,213,248,255,269,290,297,304)
for(i in 1:20){
  UE_LPS[[i]]<- cbind(UE_LPS[[i]],nuu[i])
  names(UE_LPS[[i]][,14])<- "File"
}

f.plot <- function(x,y,nr=20,nc=20, scale="raw",alpha=0.1) {
  zx = c(1:nr,rep(1,nc),1+trunc( nr*(x- min(x))/(max(x)-min(x)) ))
  zx[zx>nr] = nr 
  zy = c(rep(1,nr),1:nc,1+trunc( nc*(y- min(y))/(max(y)-min(y)) ))
  zy[zy>nc] = nc
  z = table(zx,zy); z[,1]=z[,1]-1; z[1,]=z[1,]-1; 
  if (scale=="l") z= log(1+z) 
  if (scale=="d") z= log(1+log(1+z))
  z[max(z)*alpha>z]=0
  image(z=t(z),x=seq(length=nr+1,from=min(x),to=max(x)),
        y= seq(length=nc+1,from=min(y),to=max(y)),
        xlab="",ylab="", col=topo.colors(100))
  list(z = t(z), x = seq(length=nr,from=min(x),to=max(x)), y = seq(length=nc,from=min(y),to=max(y)),zz=cbind(zx,zy)[-(1:(nr+nc)),])
}

search.ind = function(x, idx, idy, threshold = 2){
  neighbor_x = x[max((ceil(idx)-threshold),1):min((ceil(idx)+threshold),nrow(x)),
                 max((ceil(idy)-threshold),1):min((ceil(idy)+threshold),ncol(x))]
  ind = which(neighbor_x == max(neighbor_x),arr.ind = TRUE)
  idx_f = max((ceil(idx)-threshold),1)+ind[1,1]-1
  idy_f = max((ceil(idy)-threshold),1)+ind[1,2]-1
  c(idx_f, idy_f)
}

find_area <- function(x, idx, idy, t1 = 35, t2 = 0.85){
  gr = cbind(rep(1:nrow(x),ncol(x)),rep(1:ncol(x),rep(nrow(x),ncol(x))))
  d = as.matrix(dist(gr))
  id = idx+(idy-1)*nrow(x)
  gr[d[id,]<t1 & x>t2*x[id],2:1]
}

assign = rep(0,nrow(dat))
for (i in 1:nrow(dat)) {
  c_values = NULL
  for (j in 1:5) {
    c = -0.5*(t(c(f.density$x[dat[i,1]],f.density$y[dat[i,2]])-c(f.density$x[idx_f[j,1]],f.density$y[idx_f[j,2]]))%*%solve(cov[[j]])%*%(c(f.density$x[dat[i,1]],f.density$y[dat[i,2]])-c(f.density$x[idx_f[j,1]],f.density$y[idx_f[j,2]])))+log(x[idx_f][j])
    c_values = c(c_values,c)
  }
  assign[i] = which.max(c_values)
}

landmark_finder100<- function(xx){
  x0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,2])
  y0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,1])
  f.density = f.plot(x0,y0,nr=100,nc=100,scale="l")
  msize <- 100
  x <- matrix(as.numeric(f.density$z),nrow = 100,ncol = 100)
  r <- raster(x)
  extent(r) <- extent(c(0, msize, 0, msize) + 0.1)
  f <- function(X) if(max(X, na.rm=TRUE) > quantile(f.density$z,0.85)){max(X, na.rm=TRUE)}else{-1}
  ww <- matrix(1, nrow=9, ncol=9) ## Weight matrix for cells in moving window
  localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)
  r2 <- r==localmax
  maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
  maxXY_rot = cbind(100-maxXY[,2],maxXY[,1])
  idx_f = t(as.matrix(apply(as.matrix(maxXY_rot), 1, function(t){search.ind(x, t[1],t[2])})))
  idx_f = idx_f[idx_f[,1] <90 & idx_f[,2] <90,]
  idx_f = idx_f[x[idx_f] > quantile(x,0.9),]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  for (i in (nrow(idx_f)-2):1){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i+1,][1] & idx_f[i,][2]==idx_f[i+1,][2]){
      idx_f= idx_f[c(1:i,(i+2):nrow(idx_f)),]
    }
  }
  for (i in nrow(idx_f):3){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i-1,][1] & idx_f[i,][2]==idx_f[i-1,][2]){
      idx_f= idx_f[c(1:(i-2),i:nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[1,1]> (idx_f[2,1]-11) & x[idx_f][1]< x[idx_f][2]){
    idx_f= idx_f[2:nrow(idx_f),]
  }
  for (i in nrow(idx_f):2){
    #if (idx_f[i,][1]< (idx_f[1,1]+10) & idx_f[i,][2]> idx_f[1,2]){
    if (nrow(idx_f)>3 & idx_f[i,][1]< (idx_f[1,1]+11) & x[idx_f][i]<=x[idx_f][1]){
      idx_f= idx_f[c(1:(i-1),(i+1):nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[nrow(idx_f),1]< (idx_f[(nrow(idx_f)-1),1]+11) & x[idx_f][nrow(idx_f)]< x[idx_f][(nrow(idx_f)-1)]){
    idx_f= idx_f[1:(nrow(idx_f)-1),]
  }
  for (i in (nrow(idx_f)-1):1){
    #if (idx_f[i,][1]> (idx_f[nrow(idx_f),1]-11) & idx_f[i,][1]> (idx_f[(nrow(idx_f)-2),1]) & x[idx_f][i]<=x[idx_f][nrow(idx_f)]){
    if (nrow(idx_f)>3 & idx_f[i,][1]> (idx_f[nrow(idx_f),1]-11) & x[idx_f][i]<=x[idx_f][nrow(idx_f)]){
      idx_f= idx_f[-i,]
    }
  }
  if (nrow(idx_f)==4 & idx_f[2,][1]> idx_f[3,][1]-7){
    idx_f= idx_f[-3,]
  }
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  idx_f = idx_f[sort(x[idx_f],index.return = TRUE,decreasing = T)$ix[1:3],]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(cbind(f.density$x[idx_f[,1]],f.density$y[idx_f[,2]]))
  return(idx_f)
}

#number of pixels in dat threshold:
dat_finder<- function(xx){
  x0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,2])
  y0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,1])
  f.density = f.plot(x0,y0,nr=100,nc=100,scale="l")
  msize <- 100
  x <- matrix(as.numeric(f.density$z),nrow = 100,ncol = 100)
  r <- raster(x)
  extent(r) <- extent(c(0, msize, 0, msize) + 0.1)
  f <- function(X) if(max(X, na.rm=TRUE) > quantile(f.density$z,0.85)){max(X, na.rm=TRUE)}else{-1}
  ww <- matrix(1, nrow=9, ncol=9) ## Weight matrix for cells in moving window
  localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)
  r2 <- r==localmax
  maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
  maxXY_rot = cbind(100-maxXY[,2],maxXY[,1])
  idx_f = t(as.matrix(apply(as.matrix(maxXY_rot), 1, function(t){search.ind(x, t[1],t[2])})))
  idx_f = idx_f[idx_f[,1] <90 & idx_f[,2] <90,]
  idx_f = idx_f[x[idx_f] > quantile(x,0.9),]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  for (i in (nrow(idx_f)-2):1){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i+1,][1] & idx_f[i,][2]==idx_f[i+1,][2]){
      idx_f= idx_f[c(1:i,(i+2):nrow(idx_f)),]
    }
  }
  for (i in nrow(idx_f):3){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i-1,][1] & idx_f[i,][2]==idx_f[i-1,][2]){
      idx_f= idx_f[c(1:(i-2),i:nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[1,1]> (idx_f[2,1]-11) & x[idx_f][1]< x[idx_f][2]){
    idx_f= idx_f[2:nrow(idx_f),]
  }
  for (i in nrow(idx_f):2){
    if (nrow(idx_f)>3 & idx_f[i,][1]< (idx_f[1,1]+11) & x[idx_f][i]<=x[idx_f][1]){
      idx_f= idx_f[c(1:(i-1),(i+1):nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[nrow(idx_f),1]< (idx_f[(nrow(idx_f)-1),1]+11) & x[idx_f][nrow(idx_f)]< x[idx_f][(nrow(idx_f)-1)]){
    idx_f= idx_f[1:(nrow(idx_f)-1),]
  }
  for (i in (nrow(idx_f)-1):1){
    
    if (nrow(idx_f)>3 & idx_f[i,][1]> (idx_f[nrow(idx_f),1]-11) & x[idx_f][i]<=x[idx_f][nrow(idx_f)]){
      idx_f= idx_f[-i,]
    }
  }
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  idx_f = idx_f[sort(x[idx_f],index.return = TRUE,decreasing = T)$ix[1:3],]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.65)
  for (i in nrow(dat):1){
    if (dat[i,][1]>90 | dat[i,][2]>90){
      dat= dat[-i,]
    }
  }
  return(length(f.density$x[dat[,2]]))
}

qda100<- function(xx){
  x0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,2])
  y0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,1])
  f.density = f.plot(x0,y0,nr=100,nc=100,scale="l")
  msize <- 100
  x <- matrix(as.numeric(f.density$z),nrow = 100,ncol = 100)
  r <- raster(x)
  extent(r) <- extent(c(0, msize, 0, msize) + 0.1)
  f <- function(X) if(max(X, na.rm=TRUE) > quantile(f.density$z,0.85)){max(X, na.rm=TRUE)}else{-1}
  ww <- matrix(1, nrow=9, ncol=9) ## Weight matrix for cells in moving window
  localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)
  r2 <- r==localmax
  maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
  maxXY_rot = cbind(100-maxXY[,2],maxXY[,1])
  idx_f = t(as.matrix(apply(as.matrix(maxXY_rot), 1, function(t){search.ind(x, t[1],t[2])})))
  idx_f = idx_f[idx_f[,1] <90 & idx_f[,2] <90,]
  idx_f = idx_f[x[idx_f] > quantile(x,0.9),]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  for (i in (nrow(idx_f)-2):1){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i+1,][1] & idx_f[i,][2]==idx_f[i+1,][2]){
      idx_f= idx_f[c(1:i,(i+2):nrow(idx_f)),]
    }
  }
  for (i in nrow(idx_f):3){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i-1,][1] & idx_f[i,][2]==idx_f[i-1,][2]){
      idx_f= idx_f[c(1:(i-2),i:nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[1,1]> (idx_f[2,1]-11) & x[idx_f][1]< x[idx_f][2]){
    idx_f= idx_f[2:nrow(idx_f),]
  }
  for (i in nrow(idx_f):2){
    #if (idx_f[i,][1]< (idx_f[1,1]+10) & idx_f[i,][2]> idx_f[1,2]){
    if (nrow(idx_f)>3 & idx_f[i,][1]< (idx_f[1,1]+11) & x[idx_f][i]<=x[idx_f][1]){
      idx_f= idx_f[c(1:(i-1),(i+1):nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[nrow(idx_f),1]< (idx_f[(nrow(idx_f)-1),1]+11) & x[idx_f][nrow(idx_f)]< x[idx_f][(nrow(idx_f)-1)]){
    idx_f= idx_f[1:(nrow(idx_f)-1),]
  }
  for (i in (nrow(idx_f)-1):1){
    #if (idx_f[i,][1]> (idx_f[nrow(idx_f),1]-11) & idx_f[i,][1]> (idx_f[(nrow(idx_f)-2),1]) & x[idx_f][i]<=x[idx_f][nrow(idx_f)]){
    if (nrow(idx_f)>3 & idx_f[i,][1]> (idx_f[nrow(idx_f),1]-11) & x[idx_f][i]<=x[idx_f][nrow(idx_f)]){
      idx_f= idx_f[-i,]
    }
  }
  # if (nrow(idx_f)==4 & idx_f[2,][1]> idx_f[3,][1]-7){
  #   idx_f= idx_f[-3,]
  #}
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  idx_f = idx_f[sort(x[idx_f],index.return = TRUE,decreasing = T)$ix[1:3],]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  #image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  #points(cbind(f.density$x[idx_f[,1]],f.density$y[idx_f[,2]]))
  dat1 = find_area(x,idx_f[1,1],idx_f[1,2],t1 = 8, t2 = 0.8)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat1[,2]], y = f.density$y[dat1[,1]])
  dat2 = find_area(x,idx_f[2,1],idx_f[2,2],t1 = 11, t2 = 0.9)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat2[,2]], y = f.density$y[dat2[,1]])
  dat3 = find_area(x,idx_f[3,1],idx_f[3,2],t1 = 10, t2 = 0.95)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat3[,2]], y = f.density$y[dat3[,1]])
  dat123<- rbind(cbind(dat1,0),cbind(dat2,1),cbind(dat3,2))
  dat123=data.frame(dat123)
  qda2<- qda(X3~.,data=dat123)
  dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.65)
  for (i in nrow(dat):1){
    if (dat[i,][1]>90 | dat[i,][2]>90){
      dat= dat[-i,]
    }
  }
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]], y = f.density$y[dat[,1]])
  dat<-as.data.frame(dat)
  names(dat)=names(dat123)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  #plot(dat,col=as.numeric(pop123))
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]], y = f.density$y[dat[,1]],col=as.numeric(pop123)+5)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]][pop123==1], y = f.density$y[dat[,1]][pop123==1],col="black")
}

qda100_<- function(xx){
  #xx<-HEU_control[[12]]
  x0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,2])
  y0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,1])
  f.density = f.plot(x0,y0,nr=100,nc=100,scale="l")
  msize <- 100
  x <- matrix(as.numeric(f.density$z),nrow = 100,ncol = 100)
  r <- raster(x)
  extent(r) <- extent(c(0, msize, 0, msize) + 0.1)
  f <- function(X) if(max(X, na.rm=TRUE) > quantile(f.density$z,0.85)){max(X, na.rm=TRUE)}else{-1}
  ww <- matrix(1, nrow=9, ncol=9) ## Weight matrix for cells in moving window
  localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)
  r2 <- r==localmax
  maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
  maxXY_rot = cbind(100-maxXY[,2],maxXY[,1])
  idx_f = t(as.matrix(apply(as.matrix(maxXY_rot), 1, function(t){search.ind(x, t[1],t[2])})))
  idx_f = idx_f[idx_f[,1] <90 & idx_f[,2] <90,]
  idx_f = idx_f[x[idx_f] > quantile(x,0.9),]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  for (i in (nrow(idx_f)-2):1){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i+1,][1] & idx_f[i,][2]==idx_f[i+1,][2]){
      idx_f= idx_f[c(1:i,(i+2):nrow(idx_f)),]
    }
  }
  for (i in nrow(idx_f):3){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i-1,][1] & idx_f[i,][2]==idx_f[i-1,][2]){
      idx_f= idx_f[c(1:(i-2),i:nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[1,1]> (idx_f[2,1]-11) & x[idx_f][1]< x[idx_f][2]){
    idx_f= idx_f[2:nrow(idx_f),]
  }
  for (i in nrow(idx_f):2){
    #if (idx_f[i,][1]< (idx_f[1,1]+10) & idx_f[i,][2]> idx_f[1,2]){
    if (nrow(idx_f)>3 & idx_f[i,][1]< (idx_f[1,1]+11) & x[idx_f][i]<=x[idx_f][1]){
      idx_f= idx_f[c(1:(i-1),(i+1):nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[nrow(idx_f),1]< (idx_f[(nrow(idx_f)-1),1]+11) & x[idx_f][nrow(idx_f)]< x[idx_f][(nrow(idx_f)-1)]){
    idx_f= idx_f[1:(nrow(idx_f)-1),]
  }
  for (i in (nrow(idx_f)-1):1){
    #if (idx_f[i,][1]> (idx_f[nrow(idx_f),1]-11) & idx_f[i,][1]> (idx_f[(nrow(idx_f)-2),1]) & x[idx_f][i]<=x[idx_f][nrow(idx_f)]){
    if (nrow(idx_f)>3 & idx_f[i,][1]> (idx_f[nrow(idx_f),1]-11) & x[idx_f][i]<=x[idx_f][nrow(idx_f)]){
      idx_f= idx_f[-i,]
    }
  }
  # if (nrow(idx_f)==4 & idx_f[2,][1]> idx_f[3,][1]-7){
  #   idx_f= idx_f[-3,]
  #}
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  idx_f = idx_f[sort(x[idx_f],index.return = TRUE,decreasing = T)$ix[1:3],]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  #image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  #points(cbind(f.density$x[idx_f[,1]],f.density$y[idx_f[,2]]))
  dat1 = find_area(x,idx_f[1,1],idx_f[1,2],t1 = 8, t2 = 0.8)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat1[,2]], y = f.density$y[dat1[,1]])
  dat2 = find_area(x,idx_f[2,1],idx_f[2,2],t1 = 11, t2 = 0.70)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat2[,2]], y = f.density$y[dat2[,1]])
  dat3 = find_area(x,idx_f[3,1],idx_f[3,2],t1 = 10, t2 = 0.89)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat3[,2]], y = f.density$y[dat3[,1]])
  dat12<- rbind(cbind(dat1,0),cbind(dat2,1),cbind(dat3,0))
  dat12=data.frame(dat12)
  qda2<- qda(X3~.,data=dat12)
  dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.75)
  for (i in nrow(dat):1){
    if (dat[i,][1]>90 | dat[i,][2]>90){
      dat= dat[-i,]
    }
  }
  if (length(f.density$x[dat[,2]])<1000){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.73)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  if (length(f.density$x[dat[,2]])>1700){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.78)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]], y = f.density$y[dat[,1]])
  dat<-as.data.frame(dat)
  names(dat)=names(dat123)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  #plot(dat,col=as.numeric(pop123))
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]], y = f.density$y[dat[,1]],col=as.numeric(pop123)+5)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]][pop123==1], y = f.density$y[dat[,1]][pop123==1],col="deeppink2",pch=16)
}

lymph_finder_<- function(xx){
  x0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,2])
  y0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,1])
  f.density = f.plot(x0,y0,nr=100,nc=100,scale="l")
  msize <- 100
  x <- matrix(as.numeric(f.density$z),nrow = 100,ncol = 100)
  r <- raster(x)
  extent(r) <- extent(c(0, msize, 0, msize) + 0.1)
  f <- function(X) if(max(X, na.rm=TRUE) > quantile(f.density$z,0.85)){max(X, na.rm=TRUE)}else{-1}
  ww <- matrix(1, nrow=9, ncol=9) ## Weight matrix for cells in moving window
  localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)
  r2 <- r==localmax
  maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
  maxXY_rot = cbind(100-maxXY[,2],maxXY[,1])
  idx_f = t(as.matrix(apply(as.matrix(maxXY_rot), 1, function(t){search.ind(x, t[1],t[2])})))
  idx_f = idx_f[idx_f[,1] <90 & idx_f[,2] <90,]
  idx_f = idx_f[x[idx_f] > quantile(x,0.9),]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  for (i in (nrow(idx_f)-2):1){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i+1,][1] & idx_f[i,][2]==idx_f[i+1,][2]){
      idx_f= idx_f[c(1:i,(i+2):nrow(idx_f)),]
    }
  }
  for (i in nrow(idx_f):3){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i-1,][1] & idx_f[i,][2]==idx_f[i-1,][2]){
      idx_f= idx_f[c(1:(i-2),i:nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[1,1]> (idx_f[2,1]-11) & x[idx_f][1]< x[idx_f][2]){
    idx_f= idx_f[2:nrow(idx_f),]
  }
  for (i in nrow(idx_f):2){
    if (nrow(idx_f)>3 & idx_f[i,][1]< (idx_f[1,1]+11) & x[idx_f][i]<=x[idx_f][1]){
      idx_f= idx_f[c(1:(i-1),(i+1):nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[nrow(idx_f),1]< (idx_f[(nrow(idx_f)-1),1]+11) & x[idx_f][nrow(idx_f)]< x[idx_f][(nrow(idx_f)-1)]){
    idx_f= idx_f[1:(nrow(idx_f)-1),]
  }
  for (i in (nrow(idx_f)-1):1){
    if (nrow(idx_f)>3 & idx_f[i,][1]> (idx_f[nrow(idx_f),1]-11) & x[idx_f][i]<=x[idx_f][nrow(idx_f)]){
      idx_f= idx_f[-i,]
    }
  }
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  idx_f = idx_f[sort(x[idx_f],index.return = TRUE,decreasing = T)$ix[1:3],]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  dat1 = find_area(x,idx_f[1,1],idx_f[1,2],t1 = 8, t2 = 0.8)
  dat2 = find_area(x,idx_f[2,1],idx_f[2,2],t1 = 11, t2 = 0.70)
  dat3 = find_area(x,idx_f[3,1],idx_f[3,2],t1 = 10, t2 = 0.89)
  dat12<- rbind(cbind(dat1,0),cbind(dat2,1),cbind(dat3,0))
  dat12=data.frame(dat12)
  qda2<- qda(X3~.,data=dat12)
  dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.75)
  for (i in nrow(dat):1){
    if (dat[i,][1]>90 | dat[i,][2]>90){
      dat= dat[-i,]
    }
  }
  if (length(f.density$x[dat[,2]])<1000){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.73)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  if (length(f.density$x[dat[,2]])<1000){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.71)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  if (length(f.density$x[dat[,2]])>1300){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.78)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  if (length(f.density$x[dat[,2]])>1300){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.80)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  if (length(f.density$x[dat[,2]])>1300){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.82)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  if (length(f.density$x[dat[,2]])>1300){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.84)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  if (length(f.density$x[dat[,2]])>1300){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.86)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]], y = f.density$y[dat[,1]],col=as.numeric(pop123)+5)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]][pop123==1], y = f.density$y[dat[,1]][pop123==1],col="deeppink2",pch=16)
  return(length(f.density$x[dat[,2]][pop123==1]))
}

lymph_finder<- function(xx){
  x0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,2])
  y0<- sqrt((-min(min(xx[,2]),min(xx[,1])))+3+xx[,1])
  f.density = f.plot(x0,y0,nr=100,nc=100,scale="l")
  msize <- 100
  x <- matrix(as.numeric(f.density$z),nrow = 100,ncol = 100)
  r <- raster(x)
  extent(r) <- extent(c(0, msize, 0, msize) + 0.1)
  f <- function(X) if(max(X, na.rm=TRUE) > quantile(f.density$z,0.85)){max(X, na.rm=TRUE)}else{-1}
  ww <- matrix(1, nrow=9, ncol=9) ## Weight matrix for cells in moving window
  localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)
  r2 <- r==localmax
  maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
  maxXY_rot = cbind(100-maxXY[,2],maxXY[,1])
  idx_f = t(as.matrix(apply(as.matrix(maxXY_rot), 1, function(t){search.ind(x, t[1],t[2])})))
  idx_f = idx_f[idx_f[,1] <90 & idx_f[,2] <90,]
  idx_f = idx_f[x[idx_f] > quantile(x,0.9),]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  for (i in (nrow(idx_f)-2):1){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i+1,][1] & idx_f[i,][2]==idx_f[i+1,][2]){
      idx_f= idx_f[c(1:i,(i+2):nrow(idx_f)),]
    }
  }
  for (i in nrow(idx_f):3){
    if (nrow(idx_f)>3 & idx_f[i,][1]==idx_f[i-1,][1] & idx_f[i,][2]==idx_f[i-1,][2]){
      idx_f= idx_f[c(1:(i-2),i:nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[1,1]> (idx_f[2,1]-11) & x[idx_f][1]< x[idx_f][2]){
    idx_f= idx_f[2:nrow(idx_f),]
  }
  for (i in nrow(idx_f):2){
    if (nrow(idx_f)>3 & idx_f[i,][1]< (idx_f[1,1]+11) & x[idx_f][i]<=x[idx_f][1]){
      idx_f= idx_f[c(1:(i-1),(i+1):nrow(idx_f)),]
    }
  }
  if (nrow(idx_f)>3 & idx_f[nrow(idx_f),1]< (idx_f[(nrow(idx_f)-1),1]+11) & x[idx_f][nrow(idx_f)]< x[idx_f][(nrow(idx_f)-1)]){
    idx_f= idx_f[1:(nrow(idx_f)-1),]
  }
  for (i in (nrow(idx_f)-1):1){
    if (nrow(idx_f)>3 & idx_f[i,][1]> (idx_f[nrow(idx_f),1]-11) & x[idx_f][i]<=x[idx_f][nrow(idx_f)]){
      idx_f= idx_f[-i,]
    }
  }
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  idx_f = idx_f[sort(x[idx_f],index.return = TRUE,decreasing = T)$ix[1:3],]
  idx_f = idx_f[order(idx_f[,1],decreasing=FALSE),]
  dat1 = find_area(x,idx_f[1,1],idx_f[1,2],t1 = 10, t2 = 0.74)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat1[,2]], y = f.density$y[dat1[,1]])
  dat2 = find_area(x,idx_f[2,1],idx_f[2,2],t1 = 12, t2 = 0.71)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat2[,2]], y = f.density$y[dat2[,1]])
  dat3 = find_area(x,idx_f[3,1],idx_f[3,2],t1 = 12, t2 = 0.83)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat3[,2]], y = f.density$y[dat3[,1]])
  dat12<- rbind(cbind(dat1,0),cbind(dat2,1),cbind(dat3,0))
  dat12=data.frame(dat12)
  qda2<- qda(X3~.,data=dat12)
  dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.75)
  for (i in nrow(dat):1){
    if (dat[i,][1]>90 | dat[i,][2]>90){
      dat= dat[-i,]
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])<220 & length(f.density$x[dat[,2]])<1100){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.74)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])<220 & length(f.density$x[dat[,2]])<1100){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.73)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])<220 & length(f.density$x[dat[,2]])<1100){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.72)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])<220 & length(f.density$x[dat[,2]])<1100){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.71)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.77)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.79)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.81)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.83)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.85)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.86)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.87)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.88)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.89)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.90)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.91)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.92)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.93)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  if (length(f.density$x[dat[,2]][pop123==1])>300 | length(f.density$x[dat[,2]])>1150){
    dat<- find_area(x,idx_f[3,1],idx_f[3,2],t1 = 80, t2 = 0.94)
    for (i in nrow(dat):1){
      if (dat[i,][1]>90 | dat[i,][2]>90){
        dat= dat[-i,]
      }
    }
  }
  dat<-as.data.frame(dat)
  names(dat)=names(dat12)[1:2]
  pop123<- predict(qda2,newdata=dat)$class
  #for (i in length(f.density$x[dat[,2]][pop123==1]):2){
  #  if (((sort(f.density$x[dat[,2]][pop123==1]))[i]-(sort(f.density$x[dat[,2]][pop123==1]))[i-1])>9 | ((sort(f.density$y[dat[,2]][pop123==1]))[i]-(sort(f.density$y[dat[,2]][pop123==1]))[i-1])>9){
  #    f.density$x[dat[,2]][pop123==1]= f.density$x[dat[,2]][pop123==1][sort(f.density$x[dat[,2]][pop123==1])[-i]]
  #  }
  #}
  #lymph_x= f.density$x[dat[,2]][pop123==1]
  #lymph_y= f.density$y[dat[,2]][pop123==1]
  #if (((sort(f.density$x[dat[,2]][pop123==1]))[length(f.density$x[dat[,2]][pop123==1])]-(sort(f.density$x[dat[,2]][pop123==1]))[length(f.density$x[dat[,2]][pop123==1])-1])>9){
  #  lymph_x= lymph_x[-which.max(f.density$x[dat[,2]][pop123==1])]
  #   lymph_y= lymph_y[-which.max(f.density$x[dat[,2]][pop123==1])]
  # }
  # if (((sort(f.density$y[dat[,2]][pop123==1]))[length(f.density$x[dat[,2]][pop123==1])]-(sort(f.density$y[dat[,2]][pop123==1]))[length(f.density$x[dat[,2]][pop123==1])-1])>9){
  #  lymph_y= lymph_y[-which.max(f.density$y[dat[,2]][pop123==1])]
  #   lymph_x= lymph_x[-which.max(f.density$y[dat[,2]][pop123==1])]
  #  }
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]], y = f.density$y[dat[,1]],col=as.numeric(pop123)+5)
  image(z = f.density$z, x = f.density$x, y = f.density$y,col = topo.colors(100))
  points(x = f.density$x[dat[,2]][pop123==1], y = f.density$y[dat[,1]][pop123==1],col="deeppink2",pch=16)
  return(length(f.density$x[dat[,2]][pop123==1]))
}

#going back to data from pixels
nr=100
vx<- c(rep(1:nr,rep(nr,nr)))
vy<- c(rep(1:nr,nr))
cla<- pop123
lab<- as.character(dat[,1]*1000+dat[,2])
lab1<- as.character(vx*1000+vy)
names(cla)<- lab
nc=nr
zx = c(1:nr,rep(1,nc),1+trunc( nr*(x0- min(x0))/(max(x0)-min(x0)) ))
zx[zx>nr] = nr 
zy = c(rep(1,nr),1:nc,1+trunc( nc*(y0- min(y0))/(max(y0)-min(y0)) ))
zy[zy>nc] = nc
zw<- as.character(zx*1000+zy)
zw<- zw[-(1:200)]
#dat[cla==1,]
lab_lympho<-lab[cla==1]
lab_total<-match(lab1,lab_lympho)
lab_total[is.na(lab_total)]=0
#table(lab_total)
names(lab_total)=lab1
#lab_total[zw]
i<- lab_total[zw]!=0
f.density1<- f.plot(x0[i],y0[i],nr=100,nc=100,scale="l")
plot(y0,x0)
plot(y0[i],x0[i])
xx[i,]
#zz<- cla[zw]
#x0[zz==1]
#y0[zz==1]
#image(z = f.density1$z, x = f.density1$x, y = f.density1$y,col = topo.colors(100))
#zz<- lab_total[zw]
#x0[zz]
#y0[zz]
plot(y0,x0,pch=20,col="gray")
points(y0[i],x0[i],col="blue",pch=20)

