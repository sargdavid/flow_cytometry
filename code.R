rm(list = ls())
require(flowCore)
require(data.table)
require(ggplot2)
require(datanugget)
require(DNAMR)
require(Hmisc)
require(flowAI)
require(ggcyto)

tmp2 <- read.FCS(filename = "Data/HEU/0001.FCS")
HEUvsUE1_comp<- compensate(tmp2,spillover(tmp2)$SPILL)
HEUvsUE1_comp_clean<- flow_auto_qc(HEUvsUE1_comp)
trans1<- estimateLogicle(HEUvsUE1_comp_clean,colnames(HEUvsUE1_comp_clean[,3:10]))
HEUvsUE1_comp_clean_trans<- transform(HEUvsUE1_comp_clean,trans1)
HEUvsUE1<-HEUvsUE1_comp_clean_trans@exprs
#tmp8 <- read.FCS(filename = "Data/HEU/0002.FCS")
#HEUvsUE2<-tmp8@exprs
#tmp9 <- read.FCS(filename = "Data/HEU/0003.FCS")
#HEUvsUE3<-tmp9@exprs
tmp10 <- read.FCS(filename = "Data/HEU/0008.FCS")
HEUvsUE8_comp<- compensate(tmp10,spillover(tmp10)$SPILL)
HEUvsUE8_comp_clean<- flow_auto_qc(HEUvsUE8_comp)
trans8<- estimateLogicle(HEUvsUE8_comp_clean,colnames(HEUvsUE8_comp_clean[,3:10]))
HEUvsUE8_comp_clean_trans<- transform(HEUvsUE8_comp_clean,trans8)
HEUvsUE8<-HEUvsUE8_comp_clean_trans@exprs
tmp11 <- read.FCS(filename = "Data/HEU/0022.FCS")
HEUvsUE22_comp<- compensate(tmp11,spillover(tmp11)$SPILL)
HEUvsUE22_comp_clean<- flow_auto_qc(HEUvsUE22_comp)
trans22<- estimateLogicle(HEUvsUE22_comp_clean,colnames(HEUvsUE22_comp_clean[,3:10]))
HEUvsUE22_comp_clean_trans<- transform(HEUvsUE22_comp_clean,trans22)
HEUvsUE22<-HEUvsUE22_comp_clean_trans@exprs
tmp12 <- read.FCS(filename = "Data/HEU/0029.FCS")
HEUvsUE29_comp<- compensate(tmp12,spillover(tmp12)$SPILL)
HEUvsUE29_comp_clean<- flow_auto_qc(HEUvsUE29_comp)
trans29<- estimateLogicle(HEUvsUE29_comp_clean,colnames(HEUvsUE29_comp_clean[,3:10]))
HEUvsUE29_comp_clean_trans<- transform(HEUvsUE29_comp_clean,trans29)
HEUvsUE29<-HEUvsUE29_comp_clean_trans@exprs

autoplot(HEUvsUE1_comp)
autoplot(HEUvsUE1_comp_clean)
autoplot(HEUvsUE1_comp_clean_trans)

HEUvsUE1<- as.data.frame(HEUvsUE1)
HEUvsUE8<- as.data.frame(HEUvsUE8)
HEUvsUE22<- as.data.frame(HEUvsUE22)
HEUvsUE29<- as.data.frame(HEUvsUE29)

HEUvsUE<- rbind(cbind(HEUvsUE1,1),cbind(HEUvsUE8,2),cbind(HEUvsUE22,3),cbind(HEUvsUE29,4))
group<- HEUvsUE[,12]
HEUvsUE_SC<- HEUvsUE[,1:2]
HEUvsUE_SC<- as.data.frame(HEUvsUE_SC)
HEUvsUE_dyes<- HEUvsUE[,3:10]
HEUvsUE_dyes<- as.data.frame(HEUvsUE_dyes)
HEUvsUE<- HEUvsUE[,1:10]
HEUvsUE<- as.data.frame(HEUvsUE)

myfiles<- list.files(path="C:/Users/Mahandg/Downloads/flow_cytometry/Data/HEU/",pattern=".FCS$")
tmp<- flowCore::read.flowSet(myfiles,path="C:/Users/Mahandg/Downloads/flow_cytometry/Data/HEU/")


hist.data.frame(HEUvsUE)

hist(HEUvsUE$`FSC-A`,nclass=100)
hist(HEUvsUE$`SSC-A`,nclass=100)
hist(HEUvsUE$`FITC-A`,nclass=20000,xlim=c(-100,1100))
hist(HEUvsUE$`PE-A`,nclass=5000,xlim=c(-200,5000))
hist(HEUvsUE$`PerCP-Cy5-5-A`,nclass=7000,xlim=c(-200,5000))
hist(HEUvsUE$`PE-Cy7-A`,nclass=2000,xlim=c(-1000,12000))
hist(HEUvsUE$`APC-A`,nclass=7000,xlim=c(-200,4000))
hist(HEUvsUE$`APC-Cy7-A`,nclass=10000,xlim=c(-200,3000))
hist(HEUvsUE$`Pacific Blue-A`,nclass=10000,xlim=c(-200,2000))
hist(HEUvsUE$`Alex 700-A`,nclass=10000,xlim=c(-200,3000))

tr1<-transgap(HEUvsUE$`FSC-A`)
hist(tr1[[1]],nclass=100)
tr2<-transgap(HEUvsUE$`SSC-A`)
hist(tr2[[1]],nclass=100)
tr3<-trans2t(HEUvsUE$`FITC-A`)
hist(tr3[[1]],nclass=100)
tr4<-trans2t(HEUvsUE$`PE-A`)
hist(tr4[[1]],nclass=100)
tr5<-trans2t(HEUvsUE$`PerCP-Cy5-5-A`)
hist(tr5[[1]],nclass=100)
tr6<-trans2t(HEUvsUE$`PE-Cy7-A`)
hist(tr6[[1]],nclass=100)
tr7<-trans2t(HEUvsUE$`APC-A`)
hist(tr7[[1]],nclass=100)
tr8<-trans2t(HEUvsUE$`APC-Cy7-A`)
hist(tr8[[1]],nclass=100)
tr9<-trans2t(HEUvsUE$`Pacific Blue-A`)
hist(tr9[[1]],nclass=100)
tr10<-trans2t(HEUvsUE$`Alex 700-A`)
hist(tr10[[1]],nclass=100)


#col1
j= HEUvsUE[,1]>0
HEUvsUE[j,1] = log(1+0.3*HEUvsUE[j,1])
HEUvsUE[!j,1] = -log(1-0.5*HEUvsUE[!j,1])
hist(HEUvsUE[,1],100,col=7,xlim=c(5,15))
skewness(HEUvsUE$`FSC-A`)
hist(tr1[[1]],nclass=100)
skewness(tr1$x)
#col2
j= HEUvsUE[,2]>0
HEUvsUE[j,2] = log(1+0.005*HEUvsUE[j,2])
HEUvsUE[!j,2] = -log(1-0.007*HEUvsUE[!j,2])
hist(HEUvsUE[,2],100,col=7,xlim=c(0,10))
skewness(HEUvsUE$`SSC-A`)
hist(tr2[[1]],nclass=100)
skewness(tr2$x)
#protein1
j= HEUvsUE[,3]>0
HEUvsUE[j,3] = log(1+0.01*HEUvsUE[j,3])
HEUvsUE[!j,3] = -log(1-0.08*HEUvsUE[!j,3])
hist(HEUvsUE[,3],100,col=7,xlim=c(-3,4))
skewness(HEUvsUE$`FITC-A`)
hist(tr3[[1]],nclass=100)
skewness(tr3$x)
#protein2
j= HEUvsUE[,4]>0
HEUvsUE[j,4] = log(1+0.05*HEUvsUE[j,4])
HEUvsUE[!j,4] = -log(1-0.04*HEUvsUE[!j,4])
hist(HEUvsUE[,4],100,col=7)
skewness(HEUvsUE$`PE-A`)
hist(tr4[[1]],nclass=100)
skewness(tr4$x)
#protein3
j= HEUvsUE[,5]>0
HEUvsUE[j,5] = log(1+0.005*HEUvsUE[j,5])
#HEUvsUE[j,5] = 1/(HEUvsUE[j,5])
HEUvsUE[!j,5] = -log(1-0.02*HEUvsUE[!j,5])
hist(HEUvsUE[,5],100,col=7,xlim=c(-1,5))
skewness(HEUvsUE$`PerCP-Cy5-5-A`)
hist(tr5[[1]],nclass=100)
skewness(tr5$x)
#protein4
j= HEUvsUE[,6]>0
#HEUvsUE[j,6] = log(1+0.02*HEUvsUE[j,6])
HEUvsUE[j,6] = 1/(HEUvsUE[j,6])
HEUvsUE[!j,6] = -log(1-0.009*HEUvsUE[!j,6])
hist(HEUvsUE[,6],100,col=7,xlim=c(-0.5,1))
skewness(HEUvsUE$`PE-Cy7-A`)
hist(tr6[[1]],nclass=100)
skewness(tr6$x)
#protein5
j= HEUvsUE[,7]>0
HEUvsUE[j,7] = log(1+0.02*HEUvsUE[j,7])
#HEUvsUE[j,7] = 1/(HEUvsUE[j,7])
HEUvsUE[!j,7] = -log(1-0.03*HEUvsUE[!j,7])
hist(HEUvsUE[,7],100,col=7,xlim=c(-1,1))
skewness(HEUvsUE$`APC-A`)
hist(tr7[[1]],nclass=100)
skewness(tr7$x)
#protein6
j= HEUvsUE[,8]>0
HEUvsUE[j,8] = log(1+0.008*HEUvsUE[j,8])
#HEUvsUE[j,8] = 1/(HEUvsUE[j,8])
HEUvsUE[!j,8] = -log(1-0.025*HEUvsUE[!j,8])
hist(HEUvsUE[,8],100,col=7)
skewness(HEUvsUE$`APC-Cy7-A`)
hist(tr8[[1]],nclass=100)
skewness(tr8$x)
#protein7
j= HEUvsUE[,9]>0
HEUvsUE[j,9] = log(1+0.01*HEUvsUE[j,9])
#HEUvsUE[j,9] = 1/(HEUvsUE[j,9])
HEUvsUE[!j,9] = -log(1-0.1*HEUvsUE[!j,9])
hist(HEUvsUE[,9],100,col=7)
skewness(HEUvsUE$`Pacific Blue-A`)
hist(tr9[[1]],nclass=100)
skewness(tr9$x)
#protein8
j= HEUvsUE[,10]>0
#HEUvsUE[j,10] = log(1+0.01*HEUvsUE[j,10])
HEUvsUE[j,10] = 1/(HEUvsUE[j,10])
HEUvsUE[!j,10] = -log(1-0.025*HEUvsUE[!j,10])
hist(HEUvsUE[,10],100,col=7,xlim=c(-1,1))
skewness(HEUvsUE$`Alex 700-A`)
hist(tr10[[1]],nclass=100)
skewness(tr10$x)


par(mfrow=c(2,2),mar=c(3,3,1,1))
hist(HEUvsUE[group==1,4],100,col=7,xlim=c(-3,8))
hist(HEUvsUE[group==2,4],100,col=7,xlim=c(-3,8))
hist(HEUvsUE[group==3,4],100,col=7,xlim=c(-3,8))


par(mfrow=c(3,3),mar=c(5,3,1,1)) 
for(i in 1:9) hist(x[,i],100,main=nx[i])
## Doing Transformations
xx = x
d = c(0.015,0.01,0.1,0.3,0.1,0.1,0.1,0.02,0.05)
e = c(0.015,0.01,0.1,0.3,0.1,0.1,0.1,0.02,0.05)
for(i in 1:9) { 
  j= xx[,i]>0
  xx[j,i] = log(1+d[i]*xx[j,i])
  xx[!j,i] = -log(1-e[i]*xx[!j,i])
}
par(mfrow=c(3,3),mar=c(5,3,0,0)); 
for(i in 1:9) hist(xx[,i],100,col=7,main=nx[i])



my.DN = create.DN(x = HEUvsUE_dyes,
                  RS.num = 5000,
                  DN.num1 = 1389801,
                  DN.num2 = 2000,
                  dist.metric = "euclidean")

my.DN.refined = refine.DN(x = HEUvsUE_dyes,
                          DN= my.DN)

datanugg=my.DN.refined[1]
datanuggAss=my.DN.refined[2]
datanugg=data.frame(datanugg)
datanuggAss=data.frame(datanuggAss)
unref_datanugg = my.DN[1]
unref_datanuggAss = my.DN[2]
unref_datanugg=data.frame(unref_datanugg)
unref_datanuggAss=data.frame(unref_datanuggAss)
DN.info.all = unref_datanugg
DN.info.refined = datanugg
colnames(DN.info.refined)[10] = "Weight"
colnames(DN.info.all)[10] = "Weight"


# find quantiles of the data nugget weights for assigning color gradient
#DN.all.weight.q = quantile(as.numeric(DN.info.all[, "Weight"]),
#                           probs = seq(0,1-.001,.001))

#DN.refined.weight.q = quantile(as.numeric(DN.info.refined[, "Weight"]),
#                               probs = seq(0,1-.001,.001))

# initialize vectors for holding numbers for assigning color to each data nugget
#DN.all.plot.weights = rep(0, nrow(DN.info.all))
#DN.refined.plot.weights = rep(0, nrow(DN.info.refined))

# cycle through the quantiles
#for (i in 1:1000){
  
  # assign the current quantile index to all data nuggets with weights above the
  # current quantile
#  DN.all.plot.weights[which(DN.info.all[, "Weight"] > DN.all.weight.q[i])] = i 
  
#  DN.refined.plot.weights[which(DN.info.refined[, "Weight"] > 
#                                  DN.refined.weight.q[i])] = i 
  
#}

wmean = function (x, y,w) {
  ny=max(y)
  res = NULL
  for(i in 1:ny)
    res = cbind(res,   apply(x[y==i,,drop=F],2,weighted.mean,w=w[y==i]))
  res
}

#initialize
init= function(x,w=rep(1,nrow(x)),k=3, u = NULL){
  if(is.null(u)){
    n = nrow(x)
    ff = function(k,n)  sample(c(1:k,sample(k,n-k,rep=T)))
    y = ff(k,n)
    while(min(table(y))<2)  y = ff(k,n)
    u = wmean(x,y,w)}
  z = NULL
  for(i in 1:k)
    z =cbind(z,apply(w*(t(t(x)-u[,i]))^2,1,sum))
  zz=apply(z,1,which.min)
  while(length(unique(zz)) <k) zz = init(x,w,k,u)
  zz
}



#Weighted within-cluster Sum of Squares
wss4 = function(x,y,w = rep(1, length(y)),groupSum = F){
  ss0 = lm(x~factor(y),weights = w)$resid^2*w
  if(!groupSum){
    sum(ss0)
  }else{
    output = sapply(1:max(y), function(t){sum(ss0[y == t,])})
    names(output) = 1:max(y)
    sum.output = sum(output)
    return(list(WWCSS = output,
                TotalWWCSS = sum.output))
  }
}


## for column clusters , calculates cluster means 
#cluster means using wmean
# f.cmat = function (x, gr,w=rep(1,nrow(x)))
# {
#   x <- x[, sort.list(gr)]
#   w = w[sort.list(gr)]
#   tk <- c(table(gr))
#   p <- length(tk)
#   ttk <- cbind(rep(tk, p), 0)
#   ttk[1 + (p + 1) * (0:(p - 1)), 2] <- 1
#   z <- array(rep(ttk[, 2], ttk[, 1]), c(length(gr), p))
#   t(t(x %*% z)/tk)
# }

#Wkmeans for one initialization
Wkmeans_single = function(x,y,u = NULL, w = rep(1, length(y)),k=3,K=30) {
  if(missing(y)) y = init(x,w,k,u) 
  n = length(y)
  #if(missing(K)) K=2^(40/log10(n))
  k = max(y)
  tb= table(y)
  yy=y
  for(l in 1:K) {
    message(paste("Iteration ", 
                  l, sep = ""))
    uu0 = uu00 = wss4(x,yy,w) 
    for(i in 1:n) { for (j in 1:k){
      y2 = yy;
      if(y2[i]!=j & tb[y2[i]]>1) {
        y2[i]=j
        uu1= wss4(x, y2,w)
        if(uu1< uu0) {uu0= uu1; yy = y2; tb= table(yy)  }
      }
    }
    }
    if(uu0>=uu00) break;
  }
  list(uu0,yy)
}

# Function of general Wkmeans with multiple initializations
Wkmeans <- function (dataset, k, cl.centers = NULL, obs.weights, num.init = 1, 
                     max.iterations = 10, seed = 291102) 
{
  set.seed(seed)
  k = floor(k)
  num.init = floor(num.init)
  max.iterations = floor(max.iterations)
  dataset = as.matrix(dataset)
  obs.weights = as.vector(obs.weights)
  
  best.clusters = NULL
  for (U in 1:num.init) {
    message(paste("Initialization ", 
                  U, ":", sep = ""))
    result = Wkmeans_single(x = dataset, u = cl.centers, w = obs.weights,k=k,K=max.iterations)
    clusterAss = result[[2]]
    wss = result[[1]]
    print(wss)
    if (is.null(best.clusters) == TRUE) {
      best.clusters = clusterAss
      best.WWCSS.results = wss
    }
    else if (best.WWCSS.results > wss) {
      best.clusters = clusterAss
      best.WWCSS.results = wss
    }
  }
  
  best.cl.centers = t(wmean(x = dataset,y = best.clusters,w = obs.weights))
  rownames(best.cl.centers) = 1:k
  best.WWCSS.results = wss4(x = dataset, y =  best.clusters,w = obs.weights,groupSum = TRUE)
  output = list(best.clusters, best.cl.centers, best.WWCSS.results)
  names(output) = c("Cluster Assignments", "Cluster Centers","Weighted WCSS")
  return(output)
}

#predict clusters for Wkmeans
Wkmeans.predict = function(x,cl,newx) {
  k = max(cl)
  n = nrow(newx)
  clp=rep(0,k)
  clnew = rep(0,n)
  cl1 = cl
  n0 = length(cl)+1
  for(i in 1:n) {
    xx = rbind(x,newx[i,])
    for(j in 1:k) {
      cl1[n0] = j
      clp[j] = wss4(xx,cl1)$TotalWWCSS
    }
    clnew[i] = which.min(clp)
  }
  clnew
}

# align_clusters = function(x) {
#   u = cbind(1:2,2:1,c(1,3),c(3,1),2:3,3:2)
#   u = rbind(c(3,3,2,2,1,1),u)
#   m = sum(x)
#   res = NULL; for(i in 1:6) res[i] =m - sum(diag(x[u[,i],]))
#   i = which.min(res)
#   list(res[i],x[u[,i],])
# }

DN.clus_2 = Wkmeans(dataset = DN.info.refined[, 2:9],
                    k = 2,
                    obs.weights = DN.info.refined[, "Weight"],
                    num.init = 1,
                    max.iterations = 50)

DN.clus_3 = Wkmeans(dataset = DN.info.refined[, 2:9],
                    k = 3,
                    obs.weights = DN.info.refined[, "Weight"],
                    num.init = 1,
                    max.iterations = 50)

DN.clus_4 = Wkmeans(dataset = DN.info.refined[, 2:9],
                    k = 4,
                    obs.weights = DN.info.refined[, "Weight"],
                    num.init = 1,
                    max.iterations = 50)

DN.clus_5 = Wkmeans(dataset = DN.info.refined[, 2:9],
                    k = 5,
                    obs.weights = DN.info.refined[, "Weight"],
                    num.init = 1,
                    max.iterations = 50)

DN.clus_6 = Wkmeans(dataset = DN.info.refined[, 2:9],
                    k = 6,
                    obs.weights = DN.info.refined[, "Weight"],
                    num.init = 5,
                    max.iterations = 50)

DN.clus_7 = Wkmeans(dataset = DN.info.refined[, 2:9],
                    k = 7,
                    obs.weights = DN.info.refined[, "Weight"],
                    num.init = 5,
                    max.iterations = 50)

DN.clus_8 = Wkmeans(dataset = DN.info.refined[, 2:9],
                    k = 8,
                    obs.weights = DN.info.refined[, "Weight"],
                    num.init = 5,
                    max.iterations = 50)

DN.clus_9 = Wkmeans(dataset = DN.info.refined[, 2:9],
                    k = 9,
                    obs.weights = DN.info.refined[, "Weight"],
                    num.init = 5,
                    max.iterations = 50)

DN.clus_10 = Wkmeans(dataset = DN.info.refined[, 2:9],
                     k = 10,
                     obs.weights = DN.info.refined[, "Weight"],
                     num.init = 10,
                     max.iterations = 50)

DN.clus_11 = Wkmeans(dataset = DN.info.refined[, 2:9],
                     k = 11,
                     obs.weights = DN.info.refined[, "Weight"],
                     num.init = 10,
                     max.iterations = 50)

DN.clus_12 = Wkmeans(dataset = DN.info.refined[, 2:9],
                     k = 12,
                     obs.weights = DN.info.refined[, "Weight"],
                     num.init = 10,
                     max.iterations = 50)

DN.clus_13 = Wkmeans(dataset = DN.info.refined[, 2:9],
                     k = 13,
                     obs.weights = DN.info.refined[, "Weight"],
                     num.init = 10,
                     max.iterations = 50)

DN.clus_14 = Wkmeans(dataset = DN.info.refined[, 2:9],
                     k = 14,
                     obs.weights = DN.info.refined[, "Weight"],
                     num.init = 10,
                     max.iterations = 50)

DN.clus_15 = Wkmeans(dataset = DN.info.refined[, 2:9],
                     k = 15,
                     obs.weights = DN.info.refined[, "Weight"],
                     num.init = 10,
                     max.iterations = 50)

wss = c(DN.clus_2$`Weighted WCSS`$TotalWWCSS,
        DN.clus_3$`Weighted WCSS`$TotalWWCSS,
        DN.clus_4$`Weighted WCSS`$TotalWWCSS,
        DN.clus_5$`Weighted WCSS`$TotalWWCSS,
        DN.clus_6$`Weighted WCSS`$TotalWWCSS,
        DN.clus_7$`Weighted WCSS`$TotalWWCSS,
        DN.clus_8$`Weighted WCSS`$TotalWWCSS,
        DN.clus_9$`Weighted WCSS`$TotalWWCSS,
        DN.clus_10$`Weighted WCSS`$TotalWWCSS,
        DN.clus_11$`Weighted WCSS`$TotalWWCSS,
        DN.clus_12$`Weighted WCSS`$TotalWWCSS,
        DN.clus_13$`Weighted WCSS`$TotalWWCSS,
        DN.clus_14$`Weighted WCSS`$TotalWWCSS,
        DN.clus_15$`Weighted WCSS`$TotalWWCSS)
plot(2:15,wss,ylab = "Total Weighted WSS",type = "l")
points(2:15,wss)

krange = 2:15
#nd: number of differences
nd = 2
dd = if(nd == 1) wss[-1] else -diff(wss)
plot(krange[c(-1,-length(krange))],-diff(dd)/dd[-length(dd)]*100)
krange[c(-1,-length(krange))][sort.list(diff(dd)/dd[-length(dd)]*100)[1:3]]

wss = c(DN.clus_7$`Weighted WCSS`$TotalWWCSS,
        DN.clus_8$`Weighted WCSS`$TotalWWCSS,
        DN.clus_9$`Weighted WCSS`$TotalWWCSS,
        DN.clus_10$`Weighted WCSS`$TotalWWCSS,
        DN.clus_11$`Weighted WCSS`$TotalWWCSS,
        DN.clus_12$`Weighted WCSS`$TotalWWCSS,
        DN.clus_13$`Weighted WCSS`$TotalWWCSS,
        DN.clus_14$`Weighted WCSS`$TotalWWCSS,
        DN.clus_15$`Weighted WCSS`$TotalWWCSS)
plot(7:15,wss,ylab = "Total Weighted WSS",type = "l")
points(7:15,wss)

krange = 7:15
#nd: number of differences
nd = 2
dd = if(nd == 1) wss[-1] else -diff(wss)
plot(krange[c(-1,-length(krange))],-diff(dd)/dd[-length(dd)]*100,xlab = "k", ylab = "Relative 2nd difference of TWWSS")
krange[c(-1,-length(krange))][sort.list(diff(dd)/dd[-length(dd)]*100)[1:3]]

png(file="C:/Users/Mahandg/Downloads/flow_cytometry/wkmeans.png",width=1000,height = 400)
par(mfrow = c(1,2))
plot(7:15,wss,xlab = "k",ylab = "Total Weighted WSS",type = "l")
points(7:15,wss)
plot(krange[c(-1,-length(krange))],-diff(dd)/dd[-length(dd)]*100,xlab = "k", ylab = "Relative 2nd difference of TWWSS")
dev.off()

install.packages("RColorBrewer")
library(RColorBrewer)
# append the clusters to the refined data nuggets
DN.info.refined[, "Cluster"] = factor(DN.clus_12$`Cluster Assignments`)

# create separate datasets for each protein measurement that forms the data 
# nugget center
DN.info.refined1 = DN.info.refined[, c("Data.Nuggets.Center1", "Cluster")]
DN.info.refined2 = DN.info.refined[, c("Data.Nuggets.Center2", "Cluster")]
DN.info.refined3 = DN.info.refined[, c("Data.Nuggets.Center3", "Cluster")]
DN.info.refined4 = DN.info.refined[, c("Data.Nuggets.Center4", "Cluster")]
DN.info.refined5 = DN.info.refined[, c("Data.Nuggets.Center5", "Cluster")]
DN.info.refined6 = DN.info.refined[, c("Data.Nuggets.Center6", "Cluster")]
DN.info.refined7 = DN.info.refined[, c("Data.Nuggets.Center7", "Cluster")]
DN.info.refined8 = DN.info.refined[, c("Data.Nuggets.Center8", "Cluster")]

# give the separate datasets a common first column name
colnames(DN.info.refined1)[1] = "Level"
colnames(DN.info.refined2)[1] = "Level"
colnames(DN.info.refined3)[1] = "Level"
colnames(DN.info.refined4)[1] = "Level"
colnames(DN.info.refined5)[1] = "Level"
colnames(DN.info.refined6)[1] = "Level"
colnames(DN.info.refined7)[1] = "Level"
colnames(DN.info.refined8)[1] = "Level"

# create the masked protein names
protein.names = paste("Protein ",
                      LETTERS[1:8],
                      sep = "")

# create data for box plot
for.box.plot = 
  cbind.data.frame(data.frame(Protein = rep(protein.names, 
                                            each = nrow(DN.info.refined))),
                   rbind.data.frame(DN.info.refined1,
                                    DN.info.refined2,
                                    DN.info.refined3,
                                    DN.info.refined4,
                                    DN.info.refined5,
                                    DN.info.refined6,
                                    DN.info.refined7,
                                    DN.info.refined8))
BP.colors = hcl(h = seq(5, 
                        355, 
                        length = 13), 
                l = 65, 
                c = 100)[1:12]

png(file="C:/Users/Mahandg/Downloads/flow_cytometry/boxplot.png",width=1500,height = 900)
# create the box plot with clusters
ggplot(for.box.plot) +
  geom_boxplot(aes(x = Cluster,
                   y = Level,
                   fill = Cluster),
               width = 0.5) +
  facet_wrap(~Protein,
             nrow = 2) +
  scale_y_continuous("Level of Expression") +
  theme_gray(base_size = 26) +
  scale_fill_brewer(palette= "Paired") +
  guides(fill = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        element_blank())
dev.off()



wss = DN.clus_5$`Weighted WCSS`$TotalWWCSS
install.packages("RColorBrewer")
library(RColorBrewer)
# append the clusters to the refined data nuggets
DN.info.refined[, "Cluster"] = factor(DN.clus_5$`Cluster Assignments`)
DN.info.refined1 = DN.info.refined[, c("Data.Nuggets.Center1", "Cluster")]
DN.info.refined2 = DN.info.refined[, c("Data.Nuggets.Center2", "Cluster")]
DN.info.refined3 = DN.info.refined[, c("Data.Nuggets.Center3", "Cluster")]
DN.info.refined4 = DN.info.refined[, c("Data.Nuggets.Center4", "Cluster")]
DN.info.refined5 = DN.info.refined[, c("Data.Nuggets.Center5", "Cluster")]
DN.info.refined6 = DN.info.refined[, c("Data.Nuggets.Center6", "Cluster")]
DN.info.refined7 = DN.info.refined[, c("Data.Nuggets.Center7", "Cluster")]
DN.info.refined8 = DN.info.refined[, c("Data.Nuggets.Center8", "Cluster")]
colnames(DN.info.refined1)[1] = "Level"
colnames(DN.info.refined2)[1] = "Level"
colnames(DN.info.refined3)[1] = "Level"
colnames(DN.info.refined4)[1] = "Level"
colnames(DN.info.refined5)[1] = "Level"
colnames(DN.info.refined6)[1] = "Level"
colnames(DN.info.refined7)[1] = "Level"
colnames(DN.info.refined8)[1] = "Level"
# create the masked protein names
protein.names = paste("Protein ",
                      LETTERS[1:8],
                      sep = "")

# create data for box plot
for.box.plot = 
  cbind.data.frame(data.frame(Protein = rep(protein.names, 
                                            each = nrow(DN.info.refined))),
                   rbind.data.frame(DN.info.refined1,
                                    DN.info.refined2,
                                    DN.info.refined3,
                                    DN.info.refined4,
                                    DN.info.refined5,
                                    DN.info.refined6,
                                    DN.info.refined7,
                                    DN.info.refined8))
BP.colors = hcl(h = seq(5, 
                        355, 
                        length = 13), 
                l = 65, 
                c = 100)[1:12]

png(file="C:/Users/Mahandg/Downloads/flow_cytometry/boxplot.png",width=1500,height = 900)
# create the box plot with clusters
ggplot(for.box.plot) +
  geom_boxplot(aes(x = Cluster,
                   y = Level,
                   fill = Cluster),
               width = 0.5) +
  facet_wrap(~Protein,
             nrow = 2) +
  scale_y_continuous("Level of Expression") +
  theme_gray(base_size = 26) +
  scale_fill_brewer(palette= "Paired") +
  guides(fill = FALSE) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        element_blank())
dev.off()
