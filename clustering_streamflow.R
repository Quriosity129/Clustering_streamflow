
require(lubridate)

#### AUXILIARY FUNCTIONS ####

compDoY <- function(yy,mm,dd){
  if ((yy %% 4)==0){
    prevDays <- c(0,31,60,91,121,152,182,213,244,274,305,335)
  } else{
    prevDays <- c(0,31,59,90,120,151,181,212,243,273,304,334)
  }
  return(prevDays[mm]+dd)
}

myDecluster <- function(v,x,r){
  
  # v: data series (e.g. daily precipitation)
  # x: binary time series of extreme event occurrences
  # r: run length for declustering
  y <- x
  w <- rle(x)
  # Set runs of zeros with length<r to 1
  final_pos <- cumsum(w$lengths)
  for (uu in 1:(length(w$lengths)-1)){
    if ((w$lengths[uu]<r) & (w$values[uu]==0)){
      y[(final_pos[uu]-w$lengths[uu]+1):final_pos[uu]] <- 1
    }
  }
  if (x[length(x)]==1){
    uu <- w$lengths[length(w$lengths)]
    if ((w$lengths[uu]<r) & (w$values[uu]==0)){
      y[(final_pos[uu]-w$lengths[uu]+1):final_pos[uu]] <- 1
    }
  }
  w <- rle(y)
  # maximum element of each run of '1's
  final_pos <- cumsum(w$lengths)
  out_max_idx <- out_max <- out_clust_len <- out_clust_evts <- out_clust_tot <- 0*1:sum(w$values==1)
  index <- 1
  for (i in 1:length(final_pos)){
    if (w$values[i]==1){
      vec <- v[(final_pos[i]-w$lengths[i]+1):final_pos[i]]
      vec2 <- x[(final_pos[i]-w$lengths[i]+1):final_pos[i]]
      out_max[index] <- max(vec)
      out_max_idx[index] <- which.max(vec)+final_pos[i]-w$lengths[i]
      out_clust_len[index] <- length(vec)
      out_clust_evts[index] <- sum(vec2)
      out_clust_tot[index] <- sum(vec)
      index <- index+1
    }
  }
  return(list(out_max_idx,out_max,out_clust_len,out_clust_tot))
}

#### MAIN FUNCTIONS ####

classify_pr_extremes <- function(pr,t,r=2,quantile_pr=0.99){
  
  ## Inputs ##
  
  # pr: daily precipitation
  # t: time vector ('Date' or 'POSIXct' format)
  # r: declustering time window (days)
  # quantile_pr: daily precipitation quantile to define extremes
  
  doyVec <- mapply(compDoY,year(t),month(t),day(t))
  
  ## Classify extreme events according to their clustering timescale ##
  
  # Runs declustering
  ww <- 0*pr
  for (mm in 1:12){
    ww[month(t)==mm] <- 1*(pr>=quantile(pr[month(t)==mm],quantile_pr,na.rm=T))[month(t)==mm]
  }
  res <- myDecluster(data$pr,ww,r=r)
  x <- 0*pr
  x[res[[1]]] <- 1
  
  # For each extreme event: flag cluster/no cluster
  flag_1 <- flag_2 <- flag_3 <- flag_4 <- flag_6 <- flag_8 <- 0*1:sum(x)
  inds <- which(x==1)
  N <- length(x)
  index <- 1
  for (i in inds){
    flag_1[index] <- 1*(sum(x[max((i-7),1):(i-1)])>=1)
    flag_2[index] <- 1*(sum(x[max((i-14),1):(i-1)])>=1)
    flag_3[index] <- 1*(sum(x[max((i-21),1):(i-1)])>=1)
    flag_4[index] <- 1*(sum(x[max((i-28),1):(i-1)])>=1)
    flag_6[index] <- 1*(sum(x[max((i-42),1):(i-1)])>=1)
    flag_8[index] <- 1*(sum(x[max((i-56),1):(i-1)])>=1)
    index <- index+1
  }
  
  # Non-clustered events
  flag_0 <- 1-flag_8
  
  # Exclude
  flag_8 <- flag_8-flag_6
  flag_6 <- flag_6-flag_4
  flag_4 <- flag_4-flag_3
  flag_3 <- flag_3-flag_2
  flag_2 <- flag_2-flag_1
  
  # Event distribution
  extremes_sum <- c(sum(flag_0),sum(flag_1),sum(flag_2),sum(flag_3),sum(flag_4),sum(flag_6),sum(flag_8))
  
  # Result
  return(list(flag_indices=list(flag_0=inds[flag_0==1],flag_1=inds[flag_1==1],flag_2=inds[flag_2==1],flag_3=inds[flag_3==1],
                                flag_4=inds[flag_4==1],flag_5=inds[flag_6==1],flag_8=inds[flag_8==1]),event_counts=extremes_sum))
}

clustering_discharge <- function(pr,q,t,r=2,quantile_pr=0.99,quantile_q=0.95,baseflow=TRUE,L=31){
  
  ## Inputs ##
  
  # pr: daily precipitation [typically mm day**-1]
  # q: daily discharge [typically m**3 day**-1 or mm day**-1]
  # t: time vector ('Date' or 'POSIXct' format)
  # r: declustering time window [day]
  # quantile_pr: daily precipitation quantile to define extremes
  # quantile_q: daily discharge quantile to define "high discharge"
  # baseflow: if FALSE, remove the baseflow component according to Nathan and McMahon (1990)
  # L: time horizon to calculate distribution of discharge quantile/high discharge probability/... [day]
  
  doyVec <- mapply(compDoY,year(t),month(t),day(t))
  
  ## Calculate daily discharge quantiles ##
  
  if (baseflow){
    qq <- 0*q
    qq[order(q)] <- 1:length(q)/length(q)
  }else{
    q <- EcoHydRology::BaseflowSeparation(q,filter_parameter=0.925,passes=3)$qft
    qq <- 0*q
    qq[order(q)] <- 1:length(q)/length(q)
  }
  
  ## Call extreme precipitation classification ##
  
  flags = classify_pr_extremes(pr,t,r,quantile_pr)
  
  # Daily quantile distribution/flood odds ratio after events
  q0_vec <- p0_vec <- matrix(0,flags$event_counts[1],L)
  index <- 1
  for (z in flags$flag_indices$flag_0){
    q0_vec[index,] <- qq[z:(z+30)]
    p0_vec[index,] <- 1*(qq[z:(z+30)]>=quantile_q)
    index <- index+1
  }
  q1_vec <- p1_vec <- matrix(0,sum(flag_1),31)
  index <- 1
  for (z in inds[flag_1==1]){
    q1_vec[index,] <- qq[z:(z+30)]
    p1_vec[index,] <- 1*(qq[z:(z+30)]>=quantile_q)
    index <- index+1
  }
  q2_vec <- p2_vec <- matrix(0,sum(flag_2),31)
  index <- 1
  for (z in inds[flag_2==1]){
    q2_vec[index,] <- qq[z:(z+30)]
    p2_vec[index,] <- 1*(qq[z:(z+30)]>=quantile_q)
    index <- index+1
  }
  q3_vec <- p3_vec <- matrix(0,sum(flag_3),31)
  index <- 1
  for (z in inds[flag_3==1]){
    q3_vec[index,] <- qq[z:(z+30)]
    p3_vec[index,] <- 1*(qq[z:(z+30)]>=quantile_q)
    index <- index+1
  }
  q4_vec <- p4_vec <- matrix(0,sum(flag_4),31)
  index <- 1
  for (z in inds[flag_4==1]){
    q4_vec[index,] <- qq[z:(z+30)]
    p4_vec[index,] <- 1*(qq[z:(z+30)]>=quantile_q)
    index <- index+1
  }
  q6_vec <- p6_vec <- matrix(0,sum(flag_6),31)
  index <- 1
  for (z in inds[flag_6==1]){
    q6_vec[index,] <- qq[z:(z+30)]
    p6_vec[index,] <- 1*(qq[z:(z+30)]>=quantile_q)
    index <- index+1
  }
  q8_vec <- p8_vec <- matrix(0,sum(flag_8),31)
  index <- 1
  for (z in inds[flag_8==1]){
    q8_vec[index,] <- qq[z:(z+30)]
    p8_vec[index,] <- 1*(qq[z:(z+30)]>=quantile_q)
    index <- index+1
  }
  
  # Climatology
  mat0 <- array(0,c(sum(flag_0==1),31))
  index <- 1
  for (iy in doyVec[inds[flag_0==1]]){
    if (iy<=336){dsub<-iy:(iy+30)}else{dsub<-c(iy:366,1:(30-(366-iy)))}
    df <- data.frame(doy=doyVec[doyVec%in%dsub],index=1*(qq[doyVec%in%dsub]>=quantile_q))
    res <- aggregate(df$index,by=list(df$doy),mean,na.rm=T)
    n <- which(res[,1]==iy)
    if (n>1){mat0[index,] <- res$x[c(n:31,1:(n-1))]}else{mat0[index,] <- res$x}
    index <- index+1
  }
  mat1 <- array(0,c(sum(flag_1==1),31))
  index <- 1
  for (iy in doyVec[inds[flag_1==1]]){
    if (iy<=336){dsub<-iy:(iy+30)}else{dsub<-c(iy:366,1:(30-(366-iy)))}
    df <- data.frame(doy=doyVec[doyVec%in%dsub],index=1*(qq[doyVec%in%dsub]>=quantile_q))
    res <- aggregate(df$index,by=list(df$doy),mean,na.rm=T)
    n <- which(res[,1]==iy)
    if (n>1){mat1[index,] <- res$x[c(n:31,1:(n-1))]}else{mat1[index,] <- res$x}
    index <- index+1
  }
  mat2 <- array(0,c(sum(flag_2==1),31))
  index <- 1
  for (iy in doyVec[inds[flag_2==1]]){
    if (iy<=336){dsub<-iy:(iy+30)}else{dsub<-c(iy:366,1:(30-(366-iy)))}
    df <- data.frame(doy=doyVec[doyVec%in%dsub],index=1*(qq[doyVec%in%dsub]>=quantile_q))
    res <- aggregate(df$index,by=list(df$doy),mean,na.rm=T)
    n <- which(res[,1]==iy)
    if (n>1){mat2[index,] <- res$x[c(n:31,1:(n-1))]}else{mat2[index,] <- res$x}
    index <- index+1
  }
  mat3 <- array(0,c(sum(flag_3==1),31))
  index <- 1
  for (iy in doyVec[inds[flag_3==1]]){
    if (iy<=336){dsub<-iy:(iy+30)}else{dsub<-c(iy:366,1:(30-(366-iy)))}
    df <- data.frame(doy=doyVec[doyVec%in%dsub],index=1*(qq[doyVec%in%dsub]>=quantile_q))
    res <- aggregate(df$index,by=list(df$doy),mean,na.rm=T)
    n <- which(res[,1]==iy)
    if (n>1){mat3[index,] <- res$x[c(n:31,1:(n-1))]}else{mat3[index,] <- res$x}
    index <- index+1
  }
  mat4 <- array(0,c(sum(flag_4==1),31))
  index <- 1
  for (iy in doyVec[inds[flag_4==1]]){
    if (iy<=336){dsub<-iy:(iy+30)}else{dsub<-c(iy:366,1:(30-(366-iy)))}
    df <- data.frame(doy=doyVec[doyVec%in%dsub],index=1*(qq[doyVec%in%dsub]>=quantile_q))
    res <- aggregate(df$index,by=list(df$doy),mean,na.rm=T)
    n <- which(res[,1]==iy)
    if (n>1){mat4[index,] <- res$x[c(n:31,1:(n-1))]}else{mat4[index,] <- res$x}
    index <- index+1
  }
  mat6 <- array(0,c(sum(flag_6==1),31))
  index <- 1
  for (iy in doyVec[inds[flag_6==1]]){
    if (iy<=336){dsub<-iy:(iy+30)}else{dsub<-c(iy:366,1:(30-(366-iy)))}
    df <- data.frame(doy=doyVec[doyVec%in%dsub],index=1*(qq[doyVec%in%dsub]>=quantile_q))
    res <- aggregate(df$index,by=list(df$doy),mean,na.rm=T)
    n <- which(res[,1]==iy)
    if (n>1){mat6[index,] <- res$x[c(n:31,1:(n-1))]}else{mat6[index,] <- res$x}
    index <- index+1
  }
  mat8 <- array(0,c(sum(flag_8==1),31))
  index <- 1
  for (iy in doyVec[inds[flag_8==1]]){
    if (iy<=336){dsub<-iy:(iy+30)}else{dsub<-c(iy:366,1:(30-(366-iy)))}
    df <- data.frame(doy=doyVec[doyVec%in%dsub],index=1*(qq[doyVec%in%dsub]>=quantile_q))
    res <- aggregate(df$index,by=list(df$doy),mean,na.rm=T)
    n <- which(res[,1]==iy)
    if (n>1){mat8[index,] <- res$x[c(n:31,1:(n-1))]}else{mat8[index,] <- res$x}
    index <- index+1
  }
  
  # Two-sided t-test for clustered/non-clustered difference
  ptest_p1<-ptest_p2<-ptest_p3<-ptest_p4<-ptest_p6<-ptest_p8<-0*(1:L)
  for (i in 1:L){
    ptest_p1[i] <- t.test(p0_vec[,i],p1_vec[,i],alternative='two.sided',var.equal=FALSE,conf.level=0.95)$p.value
    ptest_p2[i] <- t.test(p0_vec[,i],p2_vec[,i],alternative='two.sided',var.equal=FALSE,conf.level=0.95)$p.value
    ptest_p3[i] <- t.test(p0_vec[,i],p3_vec[,i],alternative='two.sided',var.equal=FALSE,conf.level=0.95)$p.value
    ptest_p4[i] <- t.test(p0_vec[,i],p4_vec[,i],alternative='two.sided',var.equal=FALSE,conf.level=0.95)$p.value
    ptest_p6[i] <- t.test(p0_vec[,i],p6_vec[,i],alternative='two.sided',var.equal=FALSE,conf.level=0.95)$p.value
    ptest_p8[i] <- t.test(p0_vec[,i],p8_vec[,i],alternative='two.sided',var.equal=FALSE,conf.level=0.95)$p.value
  }
  
  return(list(discharge=data.frame(q_0=apply(q0_vec,2,mean,na.rm=T),q_1=apply(q1_vec,2,mean,na.rm=T),
                                   q_2=apply(q2_vec,2,mean,na.rm=T),q_3=apply(q3_vec,2,mean,na.rm=T),
                                   q_4=apply(q4_vec,2,mean,na.rm=T),q_6=apply(q6_vec,2,mean,na.rm=T),
                                   q_8=apply(q8_vec,2,mean,na.rm=T)),
              probability=data.frame(p0=apply(p0_vec,2,mean,na.rm=T),p1=apply(p1_vec,2,mean,na.rm=T),
                                     p2=apply(p2_vec,2,mean,na.rm=T),p3=apply(p3_vec,2,mean,na.rm=T),
                                     p4=apply(p4_vec,2,mean,na.rm=T),p6=apply(p6_vec,2,mean,na.rm=T),
                                     p8=apply(p8_vec,2,mean,na.rm=T)),
              odds_ratio=data.frame(or0=apply(p0_vec,2,mean,na.rm=T)*(1-colMeans(mat0))/((1-apply(p0_vec,2,mean,na.rm=T))*colMeans(mat0)),
                                    or1=apply(p1_vec,2,mean,na.rm=T)*(1-colMeans(mat1))/((1-apply(p1_vec,2,mean,na.rm=T))*colMeans(mat1)),
                                    or2=apply(p2_vec,2,mean,na.rm=T)*(1-colMeans(mat2))/((1-apply(p2_vec,2,mean,na.rm=T))*colMeans(mat2)),
                                    or3=apply(p3_vec,2,mean,na.rm=T)*(1-colMeans(mat3))/((1-apply(p3_vec,2,mean,na.rm=T))*colMeans(mat3)),
                                    or4=apply(p4_vec,2,mean,na.rm=T)*(1-colMeans(mat4))/((1-apply(p4_vec,2,mean,na.rm=T))*colMeans(mat4)),
                                    or6=apply(p6_vec,2,mean,na.rm=T)*(1-colMeans(mat6))/((1-apply(p6_vec,2,mean,na.rm=T))*colMeans(mat6)),
                                    or8=apply(p8_vec,2,mean,na.rm=T)*(1-colMeans(mat8))/((1-apply(p8_vec,2,mean,na.rm=T))*colMeans(mat8))),
              pvalues=data.frame(p1=ptest_p1,p2=ptest_p2,p3=ptest_p3,p4=ptest_p4,p6=ptest_p6,p8=ptest_p8)))
}


