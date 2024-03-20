######### Genome function include LD information #########
GenerateLDAllele <- function(freq, ld) # return gMat( M*2 ) for one individual two chromosome
{
  gMat = matrix(0, nrow=length(freq), ncol=2)
  for(i in 1:length(freq)) 
  {
    for(j in 1:2)
    {
      idx = ifelse(runif(1, 0, 1) < freq[i], 0, 1)
      if(i == 1) 
      {
        gMat[i,j] = idx 
      }
      else 
      {
        d = runif(1, 0, 1)
        a = gMat[i-1,j]
        f1 = ifelse(a == 0, freq[i-1], 1-freq[i-1])
        f2 = ifelse(a == 0, freq[i], 1-freq[i]) 
        gMat[i,j] = ifelse(d < (f1 * f2 +ld[i-1])/f1, gMat[i-1,j], 1-gMat[i-1,j]) 
      }
    }
  }
  return(gMat)
}
GenerateLDAllele_r <- function(freq, ld, r, c) # return gMat( M*2 ) for one individual two chromosome
{
  ld = ld*(1-c)^r
  gMat = matrix(0, nrow=length(freq), ncol=2)
  for(i in 1:length(freq))
  {
    for(j in 1:2)
    {
      idx = ifelse(runif(1, 0, 1) < freq[i], 0, 1) 
      if(i == 1) 
      {
        gMat[i,j] = idx 
      }
      else 
      {
        d = runif(1, 0, 1)
        a = gMat[i-1,j]
        f1 = ifelse(a == 0, freq[i-1], 1-freq[i-1])
        f2 = ifelse(a == 0, freq[i], 1-freq[i])
        gMat[i,j] = ifelse(d < (f1 * f2 +ld[i-1])/f1, gMat[i-1,j], 1-gMat[i-1,j]) 
      }
    }
  }
  return(gMat)
}
GenerateLDGeno <- function(freq, ld, N) # return unrelated g( N*M )
{
  g = matrix(0, nrow=N, ncol=length(freq))
  if(N==0) return(g)
  for(h in 1:N)
  {
    gMat = matrix(0, nrow=length(freq), ncol=2) 
    for(i in 1:length(freq))
    {
      for(j in 1:2)
      {
        idx = ifelse(runif(1, 0, 1) < freq[i], 0, 1) 
        if(i == 1) 
        {
          gMat[i,j] = idx 
        }
        else # for not first site
        {
          d = runif(1, 0, 1) 
          a = gMat[i-1,j]
          f1 = ifelse(a == 0, freq[i-1], 1-freq[i-1]) 
          f2 = ifelse(a == 0, freq[i], 1-freq[i]) 
          gMat[i,j] = ifelse(d < (f1 * f2 + ld[i-1])/f1, gMat[i-1,j], 1-gMat[i-1,j]) 
        }
      }
    }
    g[h,] = gMat[,1] + gMat[,2]
  }
  return(g)
}
GenerateLDGeno_r <- function(freq, ld, N, r, sibflag=TRUE) # return r-degree relation g = list( N*M , N*M ) 
{
  M = length(freq)
  g = list(matrix(NA, N, M), matrix(NA, N, M))
  if(N==0) return(g) 
  if(r==0)
  {
    g[[1]]=g[[2]]=GenerateLDGeno(freq, ld, N)
    return(g)
  }
  ibdscore = 0.5^r
  for(h in 1:N)
  {
    gMat1 = t(GenerateLDAllele(freq, ld)) #(2*M)
    gMat2 = t(GenerateLDAllele(freq, ld)) #(2*M)
    # sib alike
    if(sibflag)
    {
      for(j in 1:2)
      {
        ibd = sample(1:M, ceiling(M*ibdscore)) 
        gMat2[j,ibd] =  gMat1[j,ibd]
      }
    }
    # father-son alike
    else
    {
      ibd = sample(1:M, ceiling(M*ibdscore*2)) # Is ibd related to ld?
      gMat2[1,ibd] =  gMat1[1,ibd]
    }
    g[[1]][h, ] = apply(gMat1, 2, sum)
    g[[2]][h, ] = apply(gMat2, 2, sum)
  }
  return(g)
}

DprimetoD <- function(freq, Dprime)
{
  M=length(freq)
  f1=freq
  f2=1-freq
  f1r=f1[-1];f1r=c(f1r,0)
  f2r=f2[-1];f2r=c(f2r,0)
  D=apply(cbind(f1*f2r,f2*f1r), 1, min)*Dprime
  return(D)
}
######### Genome function exclude LD information #########
GenerateAllele <- function(freq) # return gMat( 2*M )
{
  M = length(freq)
  gMat = matrix(NA, 2, M)
  for(i in 1:M)
  {
    gMat[, i] = rbinom(2, 1, freq[i])
  }
  return(gMat)
}
GenerateGeno <- function(freq, N) # return unrelated g( N*M )
{
  M = length(freq)
  g = matrix(NA, N, M)
  if(N==0) return(g)
  for(i in 1:M)
  {
    g[,i] = rbinom(N, 2, freq[i])
  }
  return(g)
}
GenerateGeno_r <- function(freq, N, r, sibflag=TRUE) # return r-degree relation g = list( N*M , N*M ) 
{
  M = length(freq)
  g = list(matrix(NA, N, M), matrix(NA, N, M))
  if(N==0) return(g) 
  if(r==0)
  {
    g[[1]]=g[[2]]=GenerateGeno(freq, N)
    return(g)
  }
  ibdscore = 0.5^r
  for(h in 1:N)
  {
    gMat1 = GenerateAllele(freq) #(2*M)
    gMat2 = GenerateAllele(freq) #(2*M)
    # sib alike
    if(sibflag)
    {
      for(j in 1:2)
      {
        ibd = sample(1:M, ceiling(M*ibdscore))
        gMat2[j,ibd] =  gMat1[j,ibd]
      }
    }
    # father-son alike
    else
    {
      ibd = sample(1:M, ceiling(M*ibdscore*2))
      gMat2[1,ibd] =  gMat1[1,ibd]
    }
    g[[1]][h, ] = apply(gMat1, 2, sum)
    g[[2]][h, ] = apply(gMat2, 2, sum)
  }
  return(g)
}

############# Generate two GenoMatrix ##############
GenerateTwoGenoMatrix <- function(m, n1, n2, freq.range = c(0.05,0.5), 
                                  linkage = T, dprime.range = c(0.5,0.9))
{
  freq = runif(m, freq.range[1], freq.range[2])
  
  if(linkage) # linkage disequilibrium
  {
    Dprime = runif(M, dprime.range[1], dprime.range[2])
    ld = DprimetoD(freq, Dprime)
    X1 = GenerateLDGeno(freq, ld, n1)
    X2 = GenerateLDGeno(freq, ld, n2)
  } else {  # linkage equilibrium
    X1 = GenerateGeno(freq, n1)
    X2 = GenerateGeno(freq, n2)
  }
  
  # scale by allele frequency
  mean = apply(X1, 2, mean)/2
  X1 = t(apply(X1, 1, function(x) {(x-2*mean)/sqrt(2*mean*(1-mean))}))
  mean = apply(X2, 2, mean)/2
  X2 = t(apply(X2, 1, function(x) {(x-2*mean)/sqrt(2*mean*(1-mean))}))
  
  return(list(genoMat1 = X1, genoMat2 = X2, freq = freq))
}

GenerateGenoMatrix <- function(n, freq, flg_scale = T)
{
  X = GenerateGeno(freq, n)
  if(flg_scale){
    X = scale(X)
  }
  return(X)
}


############ Calculate me ###############
calculateMe <- function(X) # input standardized genotype matrix (N*M)
{
  M = ncol(X)
  K = tcrossprod(X)/M
  Me = 1/var(K[lower.tri(K)])
  return(Me)
}

calculateMeplink <- function(bfileprefix)
{
  plink = system("which plink",intern=T) # path to plink
  if(!file.exists(paste0(bfileprefix,".grm.gz"))){
    system(paste0(plink, " --silent --bfile ", bfileprefix, " --make-grm-gz --out ", bfileprefix))
  }
  grm = read.table(gzfile(paste0(bfileprefix,".grm.gz")), as.is = T)
  Me = 1/var(grm[grm[,1]!=grm[,2], 4], na.rm = TRUE)
  return(Me)
}

calculateMegear <- function(bfileprefix) {
  gear = system("which gear",intern=T) # path to gear
  if(!file.exists(paste0(bfileprefix,".it.me"))){
    system(paste0(gear, " --me --bfile ",bfileprefix," --iter 100 --out ",bfileprefix))
  }
  itme = read.table(paste0(bfileprefix,".it.me"), header = T)
  Me = itme[nrow(itme),"Me"]
  return(Me)
}

############# deepKin functions ###############

deepKin_minMe <- function(theta, alpha, beta){
  minMe = ( 2/theta/theta ) * (qnorm(1-alpha) + qnorm(1-beta)*(1-theta))^2
  return(minMe)
}

deepKin_deepTheta <- function(me, alpha, beta){
  # deepTheta = sqrt( 2/me ) * (qnorm(1-alpha) + qnorm(1-beta))
  deepTheta = (qnorm(1-alpha) + qnorm(1-beta))/(sqrt( me/2 ) + qnorm(1-beta))
  return(deepTheta)
}

deepKin_alpha <- function(me, theta, beta){
  alpha = pnorm(sqrt(me/2)*theta-qnorm(1-beta), lower.tail = F)
  
  return(alpha)
}

deepKin_logalpha <- function(me, theta, beta){
  logalpha = pnorm(sqrt(me/2)*theta-qnorm(1-beta), lower.tail = F, log.p = T) / log(10)
  return(logalpha)
}


extractPlinkGRM <- function(path, prefix, xcohort = F, n1, n2){
  
  if(!xcohort){  ## single cohort grm calculation
    
    if(!file.exists(paste0(path, prefix,".grm.gz"))){
      stop(paste0("Error:", path, prefix,".grm.gz and .grm.id does not exist!"))
    }
    gzfile = gzfile(paste0(path, prefix,".grm.gz"))
    grm = read.table(gzfile, as.is = T)
    
    grm.diag = as.matrix(grm[grm[,1]==grm[,2], c(4)])
    grm.tri  = as.matrix(grm[grm[,1]!=grm[,2], c(4)])  # n*(n-1)/2 vector
    # grm.id   = read.table(paste0(path, prefix,".grm.id"))

    return(list(diag = grm.diag, tri = grm.tri))
    
  } else { ## cross-cohort grm calculation 
    
    if(!file.exists(paste0(path, prefix,".rel.bin"))){
      stop(paste0("Error:", path, prefix,".rel.bin and .rel.id does not exist!"))
    }
    
    if(is.null(n1)|is.null(n2)){
      stop("Error: n1 and n2 is not specified in a cross-cohort grm extraction")
    }
    n = n1 + n2
    grm = readBin(paste0(path, prefix,".rel.bin"), what="numeric", n=n*(n+1)/2, size=4)
    G = matrix(NA, n, n)
    G[upper.tri(G, diag = T)] = grm
    
    grm.diag = as.matrix(diag(G))
    grm.tri  = G[1:n1,(n1+1):n]   # n1*n2 matrix
    # grm.id   = read.table(paste0(path, prefix,".rel.id"))
    
    return(list(diag = grm.diag, tri = grm.tri))
  }
  
}

deepKin_estimation <- function(grm.diag, grm.tri, xcohort = F, m, me, n1, n2){
  
  if(!xcohort){
    n = nrow(grm.diag)
    a = 0
    KINGX = matrix(NA, n*(n-1)/2, 1)
    for(i in 2:n){
      for(j in 1:(i-1)){
        a = a + 1
        KINGX[a,1] = 1-(grm.diag[i,1]+grm.diag[j,1])/2+grm.tri[a,1]
      }
    }
  } else {
    if(is.null(n1)|is.null(n2)){
      stop("Error: n1 and n2 is not specified in a deepkin estimation")
    }
    n = n1 + n2
    grm.diag1 = as.matrix(grm.diag[1:n1,])
    grm.diag2 = as.matrix(grm.diag[(n1+1):n,])
    
    a = 0
    KINGX = matrix(NA, n1*n2, 1)
    for(i in 1:n1){
      for(j in 1:n2){
        a = a + 1
        KINGX[a,1] = 1-(grm.diag1[i,1]+grm.diag2[j,1])/2+grm.tri[i,j]
      }
    }
  }
  
  KINGX = data.frame(KINGX)
  colnames(KINGX) = c("king")
  KINGX$king = as.numeric(KINGX$king)
  
  ## Calculate p-value
  KINGX$var =  ( 2 * (1-KINGX$king)^2 ) / me
  KINGX$t = KINGX$king / sqrt(KINGX$var)
  KINGX$`-logp` = -pnorm(KINGX$t, lower.tail = F, log.p = T) / log(10)
  
  return(KINGX[,c("king","-logp")])
}

igraph_kinship <- function(df)
{
  library(igraph)
  
  g <- graph_from_data_frame(df, directed = FALSE)
  # plot(g)
  gc <- decompose(g)
  length(gc)
  
  # gc.sub = list()
  for(i in 1:length(gc))
  {
    # gc.sub[[i]] = induced.subgraph(g, gc[[i]])
    gc[[i]] <- gc[[i]] %>% set_vertex_attr("name", value = 1:length(V(gc[[i]])))
  }
  
  output = list(gc[[1]])
  for(i in 2:length(gc))
  {
    tag = T
    for(j in 1:length(output)){
      if(identical_graphs(gc[[i]],output[[j]])){
        tag = F
        break
      }
    }
    if(tag){
      output[[length(output)+1]] = gc[[i]]
    }
  }
  
  length(output)
  
  # pdf("/public2/zqx/kingless/britishSet2/test.pdf", height = 15, width = 24)
  # par(mfrow = c(5,8))
  # for(i in 1:length(output)){
  #   plot(output[[i]])
  # }
  # dev.off()
  
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist(cormat)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
  return(cormat)
}
get_upper_tri <- function(cormat, diag){
  cormat[lower.tri(cormat, !diag)]<- NA
  return(cormat)
}
get_lower_tri <- function(cormat, diag){
  cormat[upper.tri(cormat, !diag)]<- NA
  return(cormat)
}

############# lab function ###############
lab_relations <- function(vec, n1, n2)
{
  lab <- array(9, dim = n1*n2)
  # detected relations set
  set_1 <- which(vec>=0.45 & vec<0.95)
  lab[set_1] <- 1
  return(lab)
}
lab_relations_true <- function(n1, n2, n_cp)
{
  lab <- array(0, dim = n1*n2)
  n1_0 = n1-sum(n_cp)
  n2_0 = n2-sum(n_cp)
  
  set_0 = c()
  if(n_cp[2]>0) set_1 = c()
  if(n_cp[3]>0) set_2 = c()
  if(n_cp[4]>0) set_3 = c()
  
  l0 = 1
  if(n_cp[2]>0) l1 = (n_cp[1])*n2+n_cp[1]+1
  if(n_cp[3]>0) l2 = (n_cp[1]+n_cp[2])*n2+n_cp[1]+n_cp[2]+1
  if(n_cp[4]>0) l3 = (n_cp[1]+n_cp[2]+n_cp[3])*n2+n_cp[1]+n_cp[2]+n_cp[3]+1
  
  r0 = (n_cp[1])*n2
  if(n_cp[2]>0) r1 = (n_cp[1]+n_cp[2])*n2
  if(n_cp[3]>0) r2 = (n_cp[1]+n_cp[2]+n_cp[3])*n2  
  if(n_cp[4]>0) r3 = (n_cp[1]+n_cp[2]+n_cp[3]+n_cp[4])*n2 
  
  # detected relations set
  if(l0<r0) set_0 <- seq(l0, r0, n2+1)
  if(n_cp[2]>0) { if(l1<r1) set_1 <- seq(l1, r1, n2+1) }
  if(n_cp[3]>0) { if(l2<r2) set_2 <- seq(l2, r2, n2+1) }
  if(n_cp[4]>0) { if(l3<r3) set_3 <- seq(l3, r3, n2+1) }
  
  lab[set_0] <- 1
  if(n_cp[2]>0) lab[set_1] <- 2
  if(n_cp[3]>0) lab[set_2] <- 3
  if(n_cp[4]>0) lab[set_3] <- 4
  return(lab)
}
