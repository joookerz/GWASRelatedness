library(SNPRelate)
library(ggplot2)

plink_path <- "~/Desktop/project_management/GWASrelatedness/data/euro_10Ksnp"
out_gds_file <- "~/Desktop/project_management/GWASrelatedness/data/euro_10Ksnp.gds"

# 读取 PLINK 文件并生成 GDS 文件
snpgdsBED2GDS(bed.fn = paste0(plink_path, ".bed"), bim.fn = paste0(plink_path, ".bim"), fam.fn = paste0(plink_path, ".fam"), out.gdsfn = out_gds_file)
gdsfile <- snpgdsOpen(out_gds_file)

#读取基因型矩阵
summary_data <- snpgdsSummary(gdsfile)
genotype_matrix <- snpgdsGetGeno(gdsfile)
dim(genotype_matrix)
summary_data$sample.id
#genotype_matrix <- genotype_matrix[1:20,]
                
#subset_matrix <- genotype_matrix[1:200, 1:10]
#print(subset_matrix)

#填充后矩阵
genotype_matrix1 <- genotype_matrix

m <- ncol(genotype_matrix1)

#缺失数据填充
for(i in 1:m){
  freq_table <- table(genotype_matrix1[,i])
  #print(freq_table)
  element_freq0 <- freq_table["0"]
  element_freq1 <- freq_table["1"]
  element_freq2 <- freq_table["2"]
  if(is.na(element_freq0)){
    element_freq0 = 0
  }
  if(is.na(element_freq1)){
    element_freq1 = 0
  }
  if(is.na(element_freq2)){
    element_freq2 = 0
  }
  total_count <- sum(freq_table)
  freq_value0 <- element_freq0 / total_count #P(AA)
  freq_value1 <- element_freq1 / total_count #P(Aa)
  freq_value2 <- element_freq2 / total_count #P(aa)
  freq_A <- freq_value0 + 0.5* freq_value1
  freq_a <- freq_value2 + 0.5* freq_value1
  freq_AA <- freq_A * freq_A
  freq_Aa <- 2* freq_A * freq_a
  freq_aa <- freq_a * freq_a
  #print(freq_A)
  #print(freq_a)
  #print(freq_value0)
  #print(freq_value1)
  #print(freq_value2)
  for(j in 1:200){
    k <- runif(1,min=0,max=1)
    if(k<=freq_AA){
      genotype_matrix1[j,i][is.na(genotype_matrix1[j,i])] = 0
    }
    if(k>freq_AA && k<=(freq_AA+freq_Aa)){
      genotype_matrix1[j,i][is.na(genotype_matrix1[j,i])] = 1
    }
    if(k>(freq_AA+freq_Aa) && k<=(freq_AA+freq_Aa+freq_aa)){
      genotype_matrix1[j,i][is.na(genotype_matrix1[j,i])] = 2
    }
  }
  
}

#subset_matrix1 <- genotype_matrix1[1:200, 1:10]
#print(subset_matrix1)

x <- genotype_matrix1[1:200,]
n <- nrow(x)


#scale命令标准化
genotype_matrix2 <- scale(x,center = T,scale = T)

#genotype_matrix2[,1]

#p_hat = apply(x, 2, sum)/(2*n)
#genotype_matrix2 = apply(rbind(x,p_hat), 2, function(x) {if (var(x) != 0) {
#  return((x - 2 * x[length(x)]) / sqrt(var(x)))}
#  else {
#  return(x)
# }}
#  )[1:n,]

#genotype_matrix2[,1:10]
#x[,1]

A = genotype_matrix2^2
#A[,11006:11010]
#dim(A)
V <- apply(A, 1, sum)/m
#V

#GRM
G <- (genotype_matrix2 %*% t(genotype_matrix2)) /m


theta <- matrix(data=NA, nrow=200, ncol=200)

for(i in 1:200){
  for(j in 1:200){
    theta[i,j] <- 1-0.5*(V[i]+V[j]-2*G[i,j])
  }
}

#par(mfrow = c(2, 1))
hist(theta, 200, xlim = c(-0.1,0.1))
hist(theta, 200, xlim = c(0.95,1.05))

G0 <- G[col(G)<row(G)]
G1 <- G[col(G) == row(G)]
me <- 1/var(G0)

#hist(G0)
#hist(G1)
#range(G0)
#range(G1)

Z <- matrix(data=NA, nrow=200, ncol=200)

for(i in 1:200){
  for(j in 1:200){
    Z[i,j] <-  abs(theta[i,j]/sqrt(abs(2*(1-theta[i,j]*theta[i,j])/me)))
  }
}


P1 <- 2*pnorm(Z,lower.tail=F)

hist(P1)
fdr <- p.adjust(P1, method= "BH")
hist(fdr)

log_P <- -log10(P1)
log_P[is.infinite(log_P)] <- 1e-302
plot(theta,log_P,abline(h = -log10(0.05/(n*(n-1)/2)), col = "red", lty = 2))
points_above <- which(log_P > -log10(0.05/(n*(n-1)/2)), arr.ind = TRUE)
points_above2 <- which(theta > 0.4, arr.ind = TRUE)

log_fdr <- -log10(fdr)
log_fdr[is.infinite(log_fdr)] <- 1e-302
plot(theta,log_fdr,abline(h = -log10(0.05/(n*(n-1)/2)), col = "red", lty = 2))


######################################################################################################################


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

for(i in 1:m){
  lm.fit <- lm(eg$vectors[,1]~x[,i])
  P_value[i]<- summary(lm.fit)$coefficients[2,4]
}

hist(P_value,100)
#plot(SNP,-log10(P_value))

#提取染色体信息
chromosome <- read.gdsn(index.gdsn(gdsfile, "snp.chromosome"))

#chr_AB <- lapply(chromosome, function(x) if (x %% 2 == 0) "A" else "B")


data=data.frame(s1=1:m,s2=-log10(P_value),s3=chromosome)

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
chr[1] <- chr_gap[1]/2
chr[22] <- (chr_gap[21]+100000)/2
for (i in 1:20){
  chr[j] <- (chr_gap[i]+chr_gap[i+1])/2
  j <- j+1
}

for (i in 1:m){
  data$s3[i] <- data$s3[i] %% 2 + 1 
}
summary_data$sample.id[28]
#Manhattan plot
#colors <- rainbow(length(unique(data$s3)))
colors <- c('#fa450f','#242b66')
ggplot(data = data)+geom_point(aes(x =s1, y =s2, color= as.factor(s3)), size = 0.6)+
  scale_color_manual(values = colors)+ 
  scale_y_continuous(expand = c(0,0),limits =c(0,26),breaks=seq(0,26,2))+
  scale_x_discrete(expand = c(0,0),breaks=floor(chr),labels=paste0('chr',1:22),limits=as.character(c(1:100000)))+
  theme_classic()+ 
  labs(x='chromosome',y='-log10(P_value)')+ 
  theme(legend.position = 'none')+ 
  annotate(geom = 'segment',x=0,xend=nrow(data),y=quantile(data$s2,0.95), 
           yend=quantile(data$s2,0.95),lty=4,color='black')

###################################################################################################################
#SNPselection

effect <- data.frame(s1=NA,s2=NA,s3=NA,s4=NA)
for(i in 1:60){
  P_control <- i*1e-2
  del <- as.integer(m*P_control)
  data_ranked <- data[order(data$s2, decreasing = TRUE),]
  del_number <- data_ranked$s1[1:del]
  #data_selected <- data_ranked[-(1:del),]
  #plot(data_selected$s2)
  #plot(data_ranked$s2)
  
  genotype_matrix2_hat <- genotype_matrix2[,-del_number]
  dim(genotype_matrix2_hat)
  m_hat <- ncol(genotype_matrix2_hat)
  n <- nrow(genotype_matrix2_hat)
  
  A = genotype_matrix2_hat^2
  V <- apply(A, 1, sum)/m_hat
  G <- (genotype_matrix2_hat %*% t(genotype_matrix2_hat)) /m_hat
  G0 <- G[col(G)<row(G)]
  me <- 1/var(G0)
  effect[i,]$s1 <- P_control
  effect[i,]$s2 <- m_hat
  effect[i,]$s3 <- me
  effect[i,]$s4 <- m_hat/me
}
plot(effect$s1,effect$s4,xlab="P_control",ylab="m/me")

ggplot(effect, aes(x = s1)) +
  geom_line(aes(y = s2, color = "m")) +
  geom_line(aes(y = s3, color = "me")) +
  geom_line(aes(y = s4*(3e4/2), color = "m/me")) +
  scale_color_manual(values = c("m" = "blue", "me" = "red", "m/me" = "green")) +
  labs(x = "P_control",)+scale_y_continuous(name = "m&me",sec.axis = sec_axis(trans = ~./(3e4/2), name = "m/me"))


p <- ggplot(effect, aes(x = s1)) +
  geom_point(aes(y = s4, color="m/me"))+
  scale_color_manual(values = c("m/me" = "orange"))+
  labs(x = "P_control",y = "m/me")
q <- ggplot(effect, aes(x = s1)) +
  geom_point(aes(y = s2, color = "m")) +
  geom_point(aes(y = s3*3.5, color = "me")) +
  scale_color_manual(values = c("m" = "blue", "me" = "red"))+
  labs(x = "P_control")+
  scale_y_continuous(name = "m",sec.axis = sec_axis(trans = ~./3.5, name = "me"))
library(cowplot)
plot_grid(p, q, align = 'v', ncol = 1)


ggplot(effect, aes(x = s1)) +
  geom_bar(aes(y = s2, color = "m"), stat = "identity", position = position_dodge(width = 0.4)) +
  geom_bar(aes(y = s3*3, color = "me"), stat = "identity", position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("m" = "blue", "me" = "red"))+
  labs(x = "P_value")+
  scale_y_continuous(name = "m",sec.axis = sec_axis(trans = ~./3, name = "me"))


#####################

#me_max selection
P_control <- 0.085
del <- as.integer(m*P_control)
data_ranked <- data[order(data$s2, decreasing = TRUE),]
del_number <- data_ranked$s1[1:del]
#data_selected <- data_ranked[-(1:del),]
#plot(data_selected$s2)
#plot(data_ranked$s2)

genotype_matrix2_hat <- genotype_matrix2[,-del_number]
dim(genotype_matrix2_hat)
m_hat <- ncol(genotype_matrix2_hat)
n <- nrow(genotype_matrix2_hat)

A = genotype_matrix2_hat^2
V <- apply(A, 1, sum)/m_hat
G <- (genotype_matrix2_hat %*% t(genotype_matrix2_hat)) /m_hat
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

relatedness <- data.frame(s1=summary_data$sample.id[points_above[,1]],s2=summary_data$sample.id[points_above[,2]],s3=theta1[points_above])
#unique(theta[points_above])
#summary_data$sample.id[points_above[,1]]
plot(theta,theta1, xlab="theta_fullSNP",ylab="theta_selected")


#####################


###################################################################################################################
#卡方检验

t_value=vector()
k2=vector()
for(i in 1:m){
  lm.fit <- lm(eg$vectors[,1]~x[,i])
  t_value[i]<- summary(lm.fit)$coefficients[2,3]
}

k2 <- t_value^2 / eg$values[1]


P_K <- pchisq(k2, df = 1, lower.tail = F )
#plot(SNP,-log10(P_K))
hist(P_K)

#chi-square 矫正
P_med <- median(P_K)
q1 <- qchisq(P_med, df = 1, lower.tail = F )
lambda_GC <- q1/0.455
k2_hat <- k2/lambda_GC
P_K_hat <- pchisq(k2_hat, df = 1, lower.tail = F )
hist(P_K_hat)

par(mfrow= c(1,2))
#提取染色体信息
chromosome <- read.gdsn(index.gdsn(gdsfile, "snp.chromosome"))

#chr_AB <- lapply(chromosome, function(x) if (x %% 2 == 0) "A" else "B")


data=data.frame(s1=1:m,s2=-log10(P_K_hat),s3=chromosome)

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
chr[1] <- chr_gap[1]/2
chr[22] <- (chr_gap[21]+100000)/2
for (i in 1:20){
  chr[j] <- (chr_gap[i]+chr_gap[i+1])/2
  j <- j+1
}

for (i in 1:m){
  data$s3[i] <- data$s3[i] %% 2 + 1 
}

#Manhattan plot
#colors <- rainbow(length(unique(data$s3)))
colors <- c('#fa450f','#242b66')
ggplot(data = data)+geom_point(aes(x = s1, y = s2, color = as.factor(s3)), size = 0.6)+
  scale_color_manual(values = colors)+ 
  scale_y_continuous(expand = c(0,0),limits =c(0,20),breaks=seq(0,20,2))+
  scale_x_discrete(expand = c(0,0),breaks=floor(chr),labels=paste0('chr',1:22),limits=as.character(c(1:100000)))+
  theme_classic()+ 
  labs(x='chromosome',y='-log10(P_value)')+ 
  
  theme(legend.position = 'none')+ 
  annotate(geom = 'segment',x=0,xend=nrow(data),y=quantile(data$s2,0.95), 
           yend=quantile(data$s2,0.95),lty=4,color='black')+
  annotate(geom = 'segment',x=0,xend=nrow(data),y=-log10(0.05/100000), 
           yend=-log10(0.05/100000),lty=4,color='green')
