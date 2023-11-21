library(SNPRelate)
library(ggplot2)

plink_path <- "~/Desktop/project_management/GWASrelatedness/data/CEU_TSI/CEU.2.maf.geno"
out_gds_file <- "~/Desktop/project_management/GWASrelatedness/data/CEU_TSI/CEU.2.maf.geno.gds"

# 读取 PLINK 文件并生成 GDS 文件
snpgdsBED2GDS(bed.fn = paste0(plink_path, ".bed"), bim.fn = paste0(plink_path, ".bim"), fam.fn = paste0(plink_path, ".fam"), out.gdsfn = out_gds_file)
gdsfile <- snpgdsOpen(out_gds_file)

#读取基因型矩阵
summary_data <- snpgdsSummary(gdsfile)
genotype_matrix <- snpgdsGetGeno(gdsfile)
dim(genotype_matrix)
summary_data$sample.id

#填充后矩阵
genotype_matrix1 <- genotype_matrix

#genotype_matrix[,1:10]
#判断NA
points_NA <- which(is.na(genotype_matrix1), arr.ind = TRUE)

m <- ncol(genotype_matrix1)
n <- nrow(genotype_matrix1)

x <- genotype_matrix1[1:n,]
n <- nrow(x)

p_hat = apply(x, 2, sum)/(2*n)
#print(p_hat[1:10])

#scale命令标准化
genotype_matrix2 <- scale(x,center = T,scale = T)

#x[,51212]

#genotype_matrix2 = apply(rbind(x,p_hat), 2, function(x) {if (var(x) != 0) {
#return((x - 2 * x[length(x)]) / sqrt(var(x)))}
#else {
#return(x)
#}}
#)[1:n,]

#除去重复SNP
points_NA <- which(is.na(genotype_matrix2), arr.ind = TRUE)
unique(points_NA[,2])
genotype_matrix2_hat <- genotype_matrix2[,-unique(points_NA[,2])]
dim(genotype_matrix2_hat)
m <- ncol(genotype_matrix2_hat)
n <- nrow(genotype_matrix2_hat)

#genotype_matrix2[,1:10]
#x[,1]

A = genotype_matrix2_hat^2
#A[,11006:11010]
#dim(A)
V <- apply(A, 1, sum)/m
#V

#GRM
G <- (genotype_matrix2_hat %*% t(genotype_matrix2_hat)) /m

theta <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){ 
    theta[i,j] <- 1-0.5*(V[i]+V[j]-2*G[i,j])
  }
}

hist(theta, 100)

G0 <- G[col(G)<row(G)]
G1 <- G[col(G) == row(G)]
me <- 1/var(G0)

#hist(G0)
#hist(G1)
#range(G0)
#range(G1)

Z <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){
    Z[i,j] <-  abs(theta[i,j]/sqrt(abs(2*(1-theta[i,j]*theta[i,j])/me)))
  }
}

P1 <- 2*pnorm(Z,lower.tail=F)

hist(P1)
fdr <- p.adjust(P1, method= "BH")
hist(fdr)

#提取点
log_P <- -log10(P1)
log_P[is.infinite(log_P)] <- 1e-302
plot(theta,log_P,abline(h = -log10(0.05/(n*(n-1)/2)), col = "red", lty = 2), ylim = c(0,6))
points_above <- which(log_P > -log10(0.05/(n*(n-1)/2)), arr.ind = TRUE)
points_above2 <- which(theta > 0.4, arr.ind = TRUE)

######################################################################################################################

#PCA

#write.csv(G, "C:/Users/86180/Desktop/genetic/G.csv")
#write.csv(theta, "C:/Users/86180/Desktop/genetic/theta.csv")
#write.csv(P, "C:/Users/86180/Desktop/genetic/P.csv")

eg <- eigen(G)
barplot(eg$values)
plot(eg$vectors[,1],eg$vectors[,2])

#线性回归
P_value=vector()
SNP=vector()

for(i in 1:m){
  SNP[i] <- i 
}

x_del <- x[,-unique(points_NA[,2])]
for(i in 1:m){
  lm.fit <- lm(eg$vectors[,1]~x_del[,i])
  P_value[i]<- summary(lm.fit)$coefficients[2,4]
}

#plot(SNP,-log10(P_value))

#提取染色体信息
chromosome <- read.gdsn(index.gdsn(gdsfile, "snp.chromosome"))

#chr_AB <- lapply(chromosome, function(x) if (x %% 2 == 0) "A" else "B")

chromosome1 <- chromosome[-unique(points_NA[,2])]

data=data.frame(s1=1:m,s2=-log10(P_value),s3=chromosome1)

chr_gap=vector()
chr=vector()
j <- 1
for (i in 2:m-1){
  if(data$s3[i] != data$s3[i+1]){
    chr_gap[j] <- i
    j <- j+1
  }
}

j <- 2
chr[1] <- floor(chr_gap[1]/2)
chr[22] <- floor((chr_gap[21]+m)/2)
for (i in 1:20){
  chr[j] <- floor((chr_gap[i]+chr_gap[i+1])/2)
  j <- j+1
}

#for (i in 1:m){
#  data$s3[i] <- data$s3[i] %% 2 + 1 
#}

#Manhattan plot
#colors <- rainbow(length(unique(data$s3)))
colors <- NA
colors0 <- c('#fa450f','#242b66')
for (i in 1:22){
  colors[unique(data$s3)[i]] <- colors0[i %% 2 + 1]
}

ggplot(data = data)+geom_point(aes(x =s1, y =s2, color= as.factor(s3)), size = 0.6)+
  scale_color_manual(values = colors)+ 
  scale_y_continuous(expand = c(0,0),limits =c(0,8),breaks=seq(0,8,1))+
  scale_x_discrete(expand = c(0,0),breaks=floor(chr),labels=paste0('chr',unique(data$s3)),limits=as.character(c(1:m)))+
  theme_classic()+ 
  labs(x='chromosome',y='-log10(P_value)')+ 
  theme(legend.position = 'none')+ 
  annotate(geom = 'segment',x=0,xend=nrow(data),y=quantile(data$s2,0.95), 
           yend=quantile(data$s2,0.95),lty=4,color='black')

###################################################################################################################
#SNPselection

effect <- data.frame(P_control=NA, m_hat=NA,me=NA,m_hat_me=NA)
data_ranked <- data[order(data$s2, decreasing = TRUE),]
for(i in 1:50){
  P_control <- i*1e-3
  del <- as.integer(m*P_control)
  del_number <- data_ranked$s1[1:del]
  #data_selected <- data_ranked[-(1:del),]
  #plot(data_selected$s2)
  #plot(data_ranked$s2)
  
  genotype_matrix2_hat2 <- genotype_matrix2_hat[,-del_number]
  dim(genotype_matrix2_hat2)
  m_hat <- ncol(genotype_matrix2_hat2)
  n <- nrow(genotype_matrix2_hat2)
  
  A = genotype_matrix2_hat2^2
  V <- apply(A, 1, sum)/m_hat
  G <- (genotype_matrix2_hat2 %*% t(genotype_matrix2_hat2)) /m_hat
  G0 <- G[col(G)<row(G)]
  me1 <- 1/var(G0)
  effect[i,]$P_control <- P_control

    effect[i,]$m_hat <- m_hat
  effect[i,]$me <- me1
  effect[i,]$m_hat_me <- m_hat/me1
}
plot(effect$P_control,effect$m_hat_me,xlab="P_control",ylab="m/me")

ggplot(effect, aes(x = P_control)) +
  geom_line(aes(y = m_hat/100, color = "m/100")) +
  geom_line(aes(y = me, color = "me")) +
  geom_line(aes(y = m_hat_me*(2.5e2), color = "m/me")) +
  scale_color_manual(values = c("m/100" = "blue", "me" = "red", "m/me" = "green")) +
  labs(x = "P_control",)+scale_y_continuous(name = "m&me",sec.axis = sec_axis(trans = ~./(2.5e2), name = "m/me"))

#####################

#####################

#me_max selection
P_control <- 0.019
del <- as.integer(m*P_control)
data_ranked <- data[order(data$s2, decreasing = TRUE),]
del_number <- data_ranked$s1[1:del]
#data_selected <- data_ranked[-(1:del),]
#plot(data_selected$s2)
#plot(data_ranked$s2)

genotype_matrix2_hat2 <- genotype_matrix2_hat[,-del_number]
dim(genotype_matrix2_hat2)
m_hat <- ncol(genotype_matrix2_hat2)
n <- nrow(genotype_matrix2_hat2)

A = genotype_matrix2_hat2^2
V <- apply(A, 1, sum)/m_hat
G <- (genotype_matrix2_hat2 %*% t(genotype_matrix2_hat2)) /m_hat
G0 <- G[col(G)<row(G)]
me <- 1/var(G0)

theta1 <- matrix(data=NA, nrow=n, ncol=n)
for(i in 1:n){
  for(j in 1:n){
    theta1[i,j] <- 1-0.5*(V[i]+V[j]-2*G[i,j])
  }
}

Z <- matrix(data=NA, nrow=n, ncol=n)
for(i in 1:n){
  for(j in 1:n){
    Z[i,j] <-  abs(theta1[i,j]/sqrt(abs(2*(1-theta1[i,j]*theta1[i,j])/me)))
  }
}
P1 <- 2*pnorm(Z,lower.tail=F)

log_P <- -log10(P1)
log_P[is.infinite(log_P)] <- 1e-302
plot(theta1,log_P,abline(h = -log10(0.05/(n*(n-1)/2)), col = "red", lty = 2))
points_above <- which(log_P > -log10(0.05/(n*(n-1)/2)), arr.ind = TRUE)
points_above2 <- which(theta1 > 0.4, arr.ind = TRUE)
hist(P1)
fdr <- p.adjust(P1, method= "BH")
hist(fdr)

relatedness <- data.frame(sample1=summary_data$sample.id[points_above[,1]],sample2=summary_data$sample.id[points_above[,2]],theta=theta1[points_above])
#unique(theta[points_above])
#summary_data$sample.id[points_above[,1]]
plot(theta,theta1, xlab="theta_fullSNP",ylab="theta_selected")


#####################


#####################


