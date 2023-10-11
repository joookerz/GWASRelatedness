library(SNPRelate)

# 指定 PLINK 文件的路径（不包括文件扩展名）
plink_path <- "C:/Users/86180/Desktop/genetic/euro_10Ksnp"

# 指定输出 GDS 文件的路径和文件名
out_gds_file <- "C:/Users/86180/Desktop/genetic/euro_10Ksnp.gds"

# 读取 PLINK 文件并生成 GDS 文件
snpgdsBED2GDS(bed.fn = paste0(plink_path, ".bed"), bim.fn = paste0(plink_path, ".bim"), fam.fn = paste0(plink_path, ".fam"), out.gdsfn = out_gds_file)

# 打开 GDS 文件
gdsfile <- snpgdsOpen(out_gds_file)

# 查看 GDS 文件的摘要信息
summary_data <- snpgdsSummary(gdsfile)

# 从 GDS 文件中提取基因型矩阵
genotype_matrix <- snpgdsGetGeno(gdsfile)

dim(genotype_matrix)

#genotype_matrix <- genotype_matrix[1:20,]
                
subset_matrix <- genotype_matrix[1:200, 1:10]
print(subset_matrix)

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

subset_matrix1 <- genotype_matrix1[1:200, 1:10]
print(subset_matrix1)

x <- genotype_matrix1[1:200,]
n <- nrow(x)

p_hat = apply(x, 2, sum)/(2*n)
#print(p_hat[1:10])

#scale命令标准化
genotype_matrix2 <- scale(x,center = T,scale = T)

#genotype_matrix2[,1]

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

G <- (genotype_matrix2 %*% t(genotype_matrix2)) /m

G

theta <- matrix(data=NA, nrow=200, ncol=200)

for(i in 1:200){
  for(j in 1:200){
    theta[i,j] <- V[i] + V[j] -2*G[i,j]
  }
}

G0 <- G[col(G)<row(G)]
G1 <- G[col(G) == row(G)]
me <- 1/var(G0)

hist(G0)
hist(G1)
range(G0)
range(G1)

P <- matrix(data=NA, nrow=200, ncol=200)

for(i in 1:200){
  for(j in 1:200){
    P[i,j] <-  abs(theta[i,j]/sqrt(abs(2*(1-theta[i,j]*theta[i,j])/me)))
  }
}

P1 <- pnorm(P)

write.csv(theta, "C:/Users/86180/Desktop/genetic/theta.csv")
write.csv(P, "C:/Users/86180/Desktop/genetic/P.csv")

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

plot(SNP,-log10(P_value))


