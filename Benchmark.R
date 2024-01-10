library(SNPRelate)
library(ggplot2)

plink_path <- "~/Desktop/project_management/GWASrelatedness/data/euro_10Ksnp"
out_gds_file <- "~/Desktop/project_management/GWASrelatedness/data/euro_10Ksnp.gds"

snpgdsBED2GDS(bed.fn = paste0(plink_path, ".bed"), bim.fn = paste0(plink_path, ".bim"), fam.fn = paste0(plink_path, ".fam"), out.gdsfn = out_gds_file)
gdsfile <- snpgdsOpen(out_gds_file)

summary_data <- snpgdsSummary(gdsfile)
genotype_matrix <- snpgdsGetGeno(gdsfile)
dim(genotype_matrix)
#summary_data$sample.id
#genotype_matrix <- genotype_matrix[1:20,]
                
#subset_matrix <- genotype_matrix[100:103, 2997625:2997635]
#print(subset_matrix)

#填充后矩阵
genotype_matrix1 <- genotype_matrix

m <- ncol(genotype_matrix1)
n <- nrow(genotype_matrix1)

#判断NA
points_NA <- which(is.na(genotype_matrix1), arr.ind = TRUE)

#################################################################################################
#缺失数据填
#

#
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
  for(j in 1:n){
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

x <- genotype_matrix1[1:n,]
n <- nrow(x)

#
freq <- colMeans(x)/2
hist(freq)
#
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

#genotype_matrix2[,1:10]
#x[,1]

A = genotype_matrix2^2
#A[,11006:11010]
#dim(A)
V <- apply(A, 1, sum)/m
#V
################################################################################

#GRM

G <- (genotype_matrix2 %*% t(genotype_matrix2)) /m

G0 <- G[col(G)<row(G)]
G1 <- G[col(G) == row(G)]
me_GRM <- 1/var(G0)

#hist(G0)
#hist(G1)
#range(G0)
#range(G1)

Z_GRM <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){
    Z_GRM[i,j] <- abs(G[i,j]/sqrt(abs((1+G[i,j]^2)/me)))
  }
}

P_GRM <- 2*pnorm(Z_GRM,lower.tail=F)

log_P_GRM <- -log10(P_GRM)
log_P_GRM[is.infinite(log_P_GRM)] <- 1e-302
plot(G,log_P_GRM,abline(h = -log10(0.05/(n*(n-1)/2)), col = "red", lty = 2), ylim = c(0,6))
points_above_GRM <- which(log_P_GRM > -log10(0.05/(n*(n-1)/2)), arr.ind = TRUE)
points_above2_GRM <- which(theta > 0.4, arr.ind = TRUE)
################################################################################

#encGRM
k <- 5000
S <- matrix(rnorm(m * k, mean = 0, sd = 1/sqrt(k)), nrow = m, ncol = k)

genotype_matrix3 <- genotype_matrix2 %*% S
G_enc <- (genotype_matrix3 %*% t(genotype_matrix3)) /m

G0_enc <- G_enc[col(G_enc)<row(G_enc)]
G1_enc <- G_enc[col(G_enc) == row(G_enc)]
me_enc <- 1/var(G0_enc)

#hist(G0)
#hist(G1)
#range(G0)
#range(G1)

Z_G_enc <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){
    Z_G_enc[i,j] <- abs(G_enc[i,j]/sqrt(abs((1+G_enc[i,j]^2)/me+(1+G_enc[i,j]^2)/k)))
  }
}

P_G_enc <- 2*pnorm(Z_G_enc,lower.tail=F)

log_P_G_enc <- -log10(P_G_enc)
log_P_G_enc[is.infinite(log_P_G_enc)] <- 1e-302
plot(G_enc,log_P_G_enc,abline(h = -log10(0.05/(n*(n-1)/2)), col = "red", lty = 2), ylim = c(0,6))
points_above_G_enc <- which(log_P_G_enc > -log10(0.05/(n*(n-1)/2)), arr.ind = TRUE)
points_above2_G_enc <- which(theta > 0.4, arr.ind = TRUE)
################################################################################

#encGRM_reg

G_enc_reg <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){
    lm.fit <- lm(genotype_matrix3[i,]~genotype_matrix3[j,])
    G_enc_reg[i,j]<- summary(lm.fit)$coefficients[2,1]
  }
}

G0_enc_reg <- G_enc_reg[col(G_enc_reg)<row(G_enc_reg)]
G1_enc_reg <- G_enc_reg[col(G_enc_reg) == row(G_enc_reg)]
me_enc_reg <- 1/var(G0_enc_reg)

Z_G_enc_reg <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){
    Z_G_enc_reg[i,j] <- abs(G_enc_reg[i,j]/sqrt(abs((1-G_enc_reg[i,j]^2)/me_enc_reg+(1-G_enc_reg[i,j]^2)/k)))
  }
}

P_G_enc_reg <- 2*pnorm(Z_G_enc_reg,lower.tail=F)

log_P_G_enc_reg <- -log10(P_G_enc_reg)
log_P_G_enc_reg[is.infinite(log_P_G_enc_reg)] <- 1e-302
plot(G_enc_reg,log_P_G_enc_reg,abline(h = -log10(0.05/(n*(n-1)/2)), col = "red", lty = 2), ylim = c(0,6))
points_above_G_enc_reg <- which(log_P_G_enc_reg > -log10(0.05/(n*(n-1)/2)), arr.ind = TRUE)
points_above2_G_enc_reg <- which(theta > 0.4, arr.ind = TRUE)

################################################################################

#GRM_reg

G_reg <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){
    lm.fit <- lm(genotype_matrix2[i,]~genotype_matrix2[j,])
    G_reg[i,j]<- summary(lm.fit)$coefficients[2,1]
  }
}

G0_reg <- G_reg[col(G_reg)<row(G_reg)]
G1_reg <- G_reg[col(G_reg) == row(G_reg)]
me_reg <- 1/var(G0_reg)

Z_G_reg <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){
    Z_G_reg[i,j] <- abs(G_reg[i,j]/sqrt(abs((1-G_reg[i,j]^2)/me_reg)))
  }
}

P_G_reg <- 2*pnorm(Z_G_reg,lower.tail=F)

log_P_G_reg <- -log10(P_G_reg)
log_P_G_reg[is.infinite(log_P_G_reg)] <- 1e-302
plot(G_reg,log_P_G_reg,abline(h = -log10(0.05/(n*(n-1)/2)), col = "red", lty = 2), ylim = c(0,6))
points_above_G_reg <- which(log_P_G_reg > -log10(0.05/(n*(n-1)/2)), arr.ind = TRUE)
points_above2_G_reg <- which(theta > 0.4, arr.ind = TRUE)

################################################################################

#DeepKin

theta <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){ 
    theta[i,j] <- 1-0.5*(V[i]+V[j]-2*G[i,j])
  }
}

hist(theta, 100)

G0 <- G[col(G)<row(G)]
G1 <- G[col(G) == row(G)]
me_DeepKin <- 1/var(G0)

#hist(G0)
#hist(G1)
#range(G0)
#range(G1)

Z <- matrix(data=NA, nrow=n, ncol=n)

for(i in 1:n){
  for(j in 1:n){
    Z[i,j] <-  abs(theta[i,j]/sqrt(abs(2*(1-theta[i,j])^2/me)))
  }
}

P1 <- 2*pnorm(Z,lower.tail=F)

log_P <- -log10(P1)
log_P[is.infinite(log_P)] <- 1e-302
plot(theta,log_P,abline(h = -log10(0.05/(n*(n-1)/2)), col = "red", lty = 2), ylim = c(0,6))
points_above <- which(log_P > -log10(0.05/(n*(n-1)/2)), arr.ind = TRUE)
points_above2 <- which(theta > 0.4, arr.ind = TRUE)
################################################################################

#绘比较图

P_plot <- list(P_GRM, P_G_enc, P_G_enc_reg, P_G_reg, P1)
name_row <- list("","P_GRM", "P_G_enc", "P_G_enc_reg", "P_G_reg", "P_DeepKin")
name_col <- list("P_GRM", "P_G_enc", "P_G_enc_reg", "P_G_reg", "P_DeepKin")
pdf("~/Desktop/benchmark.pdf", width = 12, height = 12) 

par(mfrow = c(6, 6))

for (i in 1:6) {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = name_row[[i]])
  text(0.5, (6 - i - 0.5), name_row[[i]], cex = 1.5, font = 2, adj = 0)
}

for (i in 1:5) {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = name_col[[i]])
  text((i - 0.5), 6, name_col[[i]], cex = 1.5, font = 2, srt = 90, adj = 1)
  for (j in 1:5) {
    plot(P_plot[[i]], P_plot[[j]], xlab = "", ylab = "")
  }
}

dev.off()

#
P_plot <- list(G, G_enc, G_enc_reg, G_reg, theta)
name_row <- list("","P_GRM", "P_G_enc", "P_G_enc_reg", "P_G_reg", "P_DeepKin")
name_col <- list("P_GRM", "P_G_enc", "P_G_enc_reg", "P_G_reg", "P_DeepKin")
pdf("~/Desktop/benchmark2.pdf", width = 12, height = 12) 

par(mfrow = c(6, 6))

for (i in 1:6) {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = name_row[[i]])
  text(0.5, (6 - i - 0.5), name_row[[i]], cex = 1.5, font = 2, adj = 0)
}

for (i in 1:5) {
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = name_col[[i]])
  text((i - 0.5), 6, name_col[[i]], cex = 1.5, font = 2, srt = 90, adj = 1)
  for (j in 1:5) {
    plot(P_plot[[i]], P_plot[[j]], xlab = "", ylab = "")
  }
}

dev.off()

hist(cor)
cor(P_G_enc[col(P_G_enc)<row(P_G_enc)],P_GRM[col(P_GRM)<row(P_GRM)])
cor(P_GRM[col(P_GRM)<row(P_GRM)],P1[col(P1)<row(P1)])
