library(tidyverse)

##make data
sim<-function(b,nr) {
  out<-list()
  for (i in 1:(length(b)-1)) for (j in (i+1):length(b)) {
    test<-runif(1)
    if (test>.5) {
      n<-sample(1:nr,1)
      sc<-rnorm(n)
      s1<-sc+b[i]+rnorm(n,sd=.1) 
      s2<-sc+b[j]+rnorm(n,sd=.1)
      out[[paste(i,j)]]<-cbind(i,j,s1,s2)
    }
  }
  df<-data.frame(do.call("rbind",out))
  names(df)<-c("t1","t2","s1","s2")
  #add name for each node
  test<-df$t1<df$t2
  t1<-df$t1
  t2<-df$t2
  s1<-df$s1
  s2<-df$s2
  df$t1<-ifelse(test,t1,t2)
  df$t2<-ifelse(test,t2,t1)
  df$s1<-ifelse(test,s1,s2)
  df$s2<-ifelse(test,s2,s1)
  df$nm<-paste(df$t1,df$t2,sep="-")
  df
}

propogate<-function(df,
                    niter=50000
) {
  ##1. get network
  net<-table(df$nm)
  nms<-strsplit(names(net),"-")
  nms<-do.call("rbind",nms)
  nbooks<-length(unique(c(nms[,1],nms[,2])))
  ##2. prepare for estimation. 
  node<-1 #always start at the first node
  ##initialize
  bhat<-rep(0,nbooks)
  bhat.N<-rep(0,nbooks)
  for (i in 2:nbooks) {
    tmp<-df[df$nm==paste("1-",i,sep=''),]
    if (nrow(tmp)>0) {
      tt<-c(tmp$t1,tmp$t2)
      ss<-c(tmp$s1,tmp$s2)
      m1<-mean(ss[tt==1])
      mm<-mean(ss[tt==i])
      bhat[i]<-mm-m1
    }
  }
  ##
  node.counter<-rep(0,nbooks)
  b.list<-list()
  b.list[[1]]<-bhat
  ncycle<-counter<-1 #counter counts all iterations, #ncycle counts iterations for a given chain
  del<-numeric()
  s<-sd(c(df$s1,df$s2))
  mdel.max<-s/250 #this is the convergence criteria. s/N of the SD of the scores
  mdel<-s*100 
  ##3. iterate
  while (counter<niter & mdel>mdel.max) {
    node.counter[node]<-node.counter[node]+1
    counter<-counter+1
    bhat<-b.list[[counter-1]]
    bhat[1]<-0
    test1<-which(nms[,1]==node)
    test2<-which(nms[,2]==node)
    subnet<-net[c(test1,test2)]
    ii<-which(rmultinom(1,1,subnet/sum(subnet))[,1]==1)
    tmp<-df[df$nm==names(subnet)[ii],]
    txt<-strsplit(names(subnet)[ii],"-")[[1]]
    nn<-which(txt!=node)
    next.node<-as.numeric(txt[nn])
    test<-tmp$t1[1]==node
    if (test) {
      s1<-tmp$s1
      s2<-tmp$s2
    } else {
      s1<-tmp$s2
      s2<-tmp$s1
    }
    ##
    if (next.node!=1) {
      #opt<-optim(bhat[next.node],obj,s1=s1,s2=s2,b.node=bhat[node],method="Brent",upper=100,lower=-100,...)
      #bhat[next.node]<-opt$par
      zz<-s2-(s1-bhat[node])
      candidate<-mean(zz)#(sum(s2)-sum(s1-bhat[node]))/length(s2)
      update.flag<-rbinom(1,1,length(s2)/(length(s2)+bhat.N[next.node]))
      if (update.flag==1) {
        bhat[next.node]<-candidate
        bhat.N[next.node]<-length(s2)
      }
    } 
    ##
    b.list[[counter]]<-bhat
    del[[counter-1]]<-sqrt(mean((bhat-b.list[[counter-1]])^2))
    if (counter>5000) {
      mdel<-mean(rev(del)[1:1000]) 
    } else mdel<-100
    ##next node
    if (ncycle>4 & counter>50 & runif(1)<.5) {
      #get unsampled nodes
      p<- 1-node.counter/sum(node.counter)
      p<-p/sum(p)
      node<-which(rmultinom(1,1,p)[,1]==1)
      ncycle<-1
    }
    else {
      ncycle<-ncycle+1
      node<-next.node
    }
    print(c(counter,node,mdel/mdel.max))
  }
  est<-do.call("rbind",b.list)
  #list(est=est,node.counter=node.counter)
  b.list
}

data.prep <- function(df){
  n_forms <- max(c(df$t1, df$t2))
  theta <- rep(0, n_forms)
  delta <- matrix(0, ncol=n_forms, nrow=n_forms)
  N <- matrix(0, ncol=n_forms, nrow=n_forms)
  
  for (i in 1:nrow(df)){
    delta[df$t1[i], df$t2[i]] <- delta[df$t1[i], df$t2[i]] + (df$s1[i]-df$s2[i]) 
    delta[df$t2[i], df$t1[i]] <- delta[df$t2[i], df$t1[i]] + (df$s2[i]-df$s1[i]) 
    N[df$t1[i], df$t2[i]] <- N[df$t1[i], df$t2[i]] + 1 
    N[df$t2[i], df$t1[i]] <- N[df$t2[i], df$t1[i]] + 1 
  }
  
  delta <- delta/N
  delta[is.nan(delta)] <- 0
  out <- list(theta = theta, delta = delta, N = N)
  return(out)
}

network.equate <- function(mats, steps=1e5, eps=1e-5){
  theta <- mats$theta
  delta <- mats$delta
  N <- mats$N
  i <- change <- 1
  while (i < steps && change > eps){
    p_theta <- theta
    theta <- rowSums(t(theta - delta)*N) / rowSums(N)
    theta[1] <- 0
    change <- sqrt(sum(theta-p_theta)^2)
    cat(i, change, '\n')
    i <- i + 1
  }
  return(theta)
}

# generate data

set.seed(8675309)
N<-200
sigma<-1
nr<-200
b<-c(0,rnorm(N,mean=.1,sd=sigma))
df<-sim(b,nr=nr)


# Klint's solution

ptm <- proc.time()
mats <- data.prep(df)
theta <- network.equate(mats)
klint_time <- proc.time() - ptm

# Ben's solution

ptm <- proc.time()
b.list<-propogate(df,niter=40000)
est<-do.call("rbind",b.list)
ben_time <- proc.time() - ptm

# klint_time
# ben_time
#
# data.frame(ben = est[nrow(est),],
#            klint = theta) |> 
#   ggplot(aes(x = ben, y = klint)) +
#   geom_point(color='firebrick1', alpha=0.6) +
#   geom_abline(aes(slope=1, intercept=0), lty=2) +
#   theme_bw()
