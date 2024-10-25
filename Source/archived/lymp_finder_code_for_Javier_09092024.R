# Functions

## f.plot

2D density plots
```{r}
f.plot <- function(x,
                   y,
                   nr = 20,
                   nc = 20,
                   scale = "raw",
                   alpha = 0.1) {
  zx <- c(1:nr,
          rep(1,
              nc),
          1 + trunc(nr*(x - min(x))/(max(x) - min(x))))
  zx[zx > nr] <- nr 
  zy <- c(rep(1,
              nr),
          1:nc,
          1 + trunc(nc*(y - min(y))/(max(y) - min(y)) ))
  zy[zy > nc] <- nc
  z <- table(zx,
             zy)
  z[,1] <- z[,1] - 1
  z[1,] <- z[1,] - 1 
  if (scale == "l") z <- log(1+z) 
  if (scale == "d") z <- log(1+log(1+z))
  z[max(z)*alpha > z] <- 0
  list(z = t(z), 
       x = seq(length = nr,
               from = min(x),
               to = max(x)),
       y = seq(length = nc,
               from = min(y),
               to = max(y)),
       zz = cbind(zx,
                  zy)[-(1:(nr + nc)),],
       zx = zx,
       zy = zy)
}
```

## search.ind

Finding landmarks
```{r}
search.ind <- function(x,
                       idx,
                       idy,
                       threshold = 2) {
  neighbor_x <- x[max((ceil(idx) - threshold), 1):min((ceil(idx) + threshold),
                                                      nrow(x)),
                  max((ceil(idy) - threshold), 1):min((ceil(idy) + threshold),
                                                      ncol(x))]
  ind <- which(neighbor_x == max(neighbor_x),
               arr.ind = TRUE)
  idx_f <- max((ceil(idx) - threshold),1) + ind[1,1] - 1
  idy_f <- max((ceil(idy) - threshold),1) + ind[1,2] - 1
  c(idx_f,
    idy_f)
}
```

## find_area

Identifying clusters around each landmark (i.e. cell subpopulations in FSC vs SSC plot)
```{r}
find_area <- function(x,
                      idx,
                      idy,
                      t1 = 35, # distance threshold from a landmark
                      t2 = 0.85) { # density threshold
  gr <- cbind(rep(1:nrow(x),
                  ncol(x)),
              rep(1:ncol(x),
                  rep(nrow(x),
                      ncol(x))))
  d <- as.matrix(dist(gr))
  id <- idx+(idy-1)*nrow(x)
  gr[d[id,] < t1 & x > t2*x[id],2:1]
}
```

## lympho_finder_convex

NOTE: make this generic to find any subpopulation in FSC vs. SSC plot
```{r}
lympho_finder_convex<- function(fsc_a, 
                                ssc_a) { 
  x0 <- sqrt((-min(min(ssc_a),
                   min(fsc_a))) + 3 + ssc_a)
  y0 <- sqrt((-min(min(ssc_a),
                   min(fsc_a))) + 3 + fsc_a)
  f.density <- f.plot(x0,
                      y0,
                      nr = 100,
                      nc=100,
                      scale="l")
  msize <- 100
  x <- matrix(as.numeric(f.density$z),
              nrow = 100,
              ncol = 100)
  r <- raster(x)
  extent(r) <- extent(c(0,
                        msize,
                        0,
                        msize) + 0.1)
  f <- function(X) if(max(X, na.rm = TRUE) > quantile(f.density$z,0.85)){max(X, na.rm = TRUE)} else {-1}
  ww <- matrix(1, 
               nrow = 9,
               ncol = 9) # weight matrix for cells in moving window
  localmax <- focal(r, 
                    fun = f,
                    w = ww,
                    pad = TRUE,
                    padValue = NA)
  r2 <- r==localmax
  maxXY <- xyFromCell(r2,
                      Which(r2 == 1, 
                            cells = TRUE))
  maxXY_rot <- cbind(100 - maxXY[, 2],
                     maxXY[, 1])
  idx_f <- t(as.matrix(apply(as.matrix(maxXY_rot),
                             1,
                             function(t) {search.ind(x, t[1],t[2])})))
  idx_f <- idx_f[idx_f[,1] <90 & 
                   idx_f[,2] <90,]
  idx_f <- idx_f[x[idx_f] > quantile(x,0.90),]
  idx_f <- unique(idx_f)
  for (i in (nrow(idx_f)-1):1) {
    if (isTRUE(nrow(idx_f) > 3 &
               (idx_f[i,][1] >= (idx_f[(i+1),][1]-8)) &
               (idx_f[i,][2]>=(idx_f[(i+1),][2]-12)) &
               (idx_f[i,][2]<=(idx_f[(i+1),][2]+12)) &
               (x[idx_f][i]<x[idx_f][(i+1)]))){
      idx_f= idx_f[-i,]
    } else if (isTRUE(nrow(idx_f) > 3 &
                      (idx_f[i,][1]>=(idx_f[(i+1),][1]-8)) &
                      (idx_f[i,][2]>=(idx_f[(i+1),][2]-12)) &
                      (idx_f[i,][2]<=(idx_f[(i+1),][2]+12)) &
                      (x[idx_f][i]>=x[idx_f][(i+1)]))) {
      idx_f= idx_f[-(i+1),]
    }
  }
  for (i in (nrow(idx_f) - 2):1){
    if (isTRUE((nrow(idx_f) > 3 & (idx_f[i,][1]>=(idx_f[(i+2),][1]-8)) &
                (idx_f[i,][2] >= (idx_f[(i+2),][2]-12)) &
                (idx_f[i,][2] <= (idx_f[(i+2),][2]+12)) &
                (x[idx_f][i] < x[idx_f][(i+2)])))) {
      idx_f= idx_f[-i,]
    }
    else if (isTRUE(nrow(idx_f) > 3 &
                    (idx_f[i,][1] >= (idx_f[(i+2),][1]-8)) &
                    (idx_f[i,][2] >= (idx_f[(i+2),][2]-12)) &
                    (idx_f[i,][2] <= (idx_f[(i+2),][2]+12)) &
                    (x[idx_f][i] >= x[idx_f][(i+2)]))) {
      idx_f= idx_f[-(i+2),]
    }
  }
  for (i in (nrow(idx_f) - 3):1) {
    if (nrow(idx_f) > 3 &
        (idx_f[i,][1] >= (idx_f[(i+3),][1]-8)) &
        (idx_f[i,][2] >= (idx_f[(i+3),][2]-12)) &
        (idx_f[i,][2] <= (idx_f[(i+3),][2]+12)) &
        (x[idx_f][i] <=x [idx_f][(i+3)])){
      idx_f= idx_f[-i,]
    }
    else if (nrow(idx_f) > 3 &
             (idx_f[i,][1] >= (idx_f[(i+3),][1]-8)) &
             (idx_f[i,][2] >= (idx_f[(i+3),][2]-12)) &
             (idx_f[i,][2] <= (idx_f[(i+3),][2]+12)) &
             (x[idx_f][i] > x[idx_f][(i+3)])) {
      idx_f= idx_f[-(i+3),]
    }
  }
  idx_f <- idx_f[sort(x[idx_f],
                      index.return = TRUE,
                      decreasing = T)$ix,]
  for (i in nrow(idx_f):2) {
    if (nrow(idx_f) > 3 &
        idx_f[i,][1] > idx_f[1,][1] &
        idx_f[i,][1] < (idx_f[1,][1]+15) &
        (idx_f[i,][2] <= (idx_f[1,][2]+20))) {
      idx_f= idx_f[-i,]
    }
  }
  idx_f <- idx_f[sort(x[idx_f],
                      index.return = TRUE,
                      decreasing = T)$ix[1:3],]
  idx_f <- idx_f[order(idx_f[,1],
                       decreasing=FALSE),]
  idx_f <- idx_f[sort(x[idx_f],
                      index.return = TRUE,
                      decreasing = T)$ix,]
  dat1 <- find_area(x,
                    idx_f[1,1],
                    idx_f[1,2],
                    t1 = 8,
                    t2 = 0.75)
  dat2 <- find_area(x,
                    idx_f[2,1],
                    idx_f[2,2],
                    t1 = 6,
                    t2 = 0.80)
  dat3 <- find_area(x,
                    idx_f[3,1],
                    idx_f[3,2],
                    t1 = 6,
                    t2 = 0.85)
  dat12 <- rbind(cbind(dat1,1),
                 cbind(dat2,0),
                 cbind(dat3,0))
  dat12 <- data.frame(dat12)
  qda2 <- qda(X3~.,
              data=dat12)
  dat_conv <- find_area(x,
                        idx_f[1,1],
                        idx_f[1,2],
                        t1 = 50,
                        t2 = 0.4)
  datt <- cbind(f.density$x[dat_conv[,2]],
                f.density$y[dat_conv[,1]])
  dat <- find_area(x,idx_f[1,1],
                   idx_f[1,2],
                   t1 = 80,
                   t2 = 0.61)
  for (i in nrow(dat):1){
    if (dat[i,][1] > 90 | 
        dat[i,][2] > 90){
      dat= dat[-i,]
    }
  }
  dat <- as.data.frame(dat)
  names(dat) <- names(dat12)[1:2]
  pop123 <- predict(qda2,
                    newdata = dat)$class
  v <- dat[pop123 == 1,]
  m <- as.matrix(dist(v))
  mm <- t(apply(m,2,function(x) c(sum(x ==1),sum(x==sqrt(2)))))
  j3 <- mm[,1] == 0 |
    (mm[,1] == 1 &
       mm[,2] < 2)
  v1 <- v[!j3,]
  m <- as.matrix(dist(v1))
  mm <- t(apply(m,2,function(x) c(sum(x ==1),sum(x==sqrt(2)))))
  j3 <-  mm[,1]==0 | ( mm[,1]==1 & mm[,2]<2)
  v2 <- v1[!j3,]
  chnew <- convhulln(v2)
  jj <- inhulln(chnew, f.density$zz)
  xx <- as.matrix(data.table(fsc_a = fsc_a,
                             ssc_a = ssc_a))
  chnewnew <- convhulln(xx[jj, ])
  return(inhulln(chnewnew,
                 xx))
}
```

## Separate lymphosites

```{r,fig.width=5,fig.height=5}
row_keep <- lapply(X = dt_exp,
                   FUN = function(a) {
                     out <- lympho_finder_convex(fsc_a = a[, 1],
                                                 ssc_a = a[, 2])
                     print(table(out))
                     return(out)
                   })
```