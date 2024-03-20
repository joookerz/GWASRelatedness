nsnp <- 50000
nind <- 200
nums <- sample(1:2, nsnp * 2 * nind, replace = T)
snp_matrix <- matrix(nums, nrow = 200)

col_idx <- matrix(1:100000, ncol = 2, byrow = T)

dim(snp_matrix)

base <- sample(c("A", "T", "C", "G"), 2)

for (i in 1:nrow(col_idx)) {
  base <- sample(c("A", "T", "C", "G"), 2)
  snp_matrix[,col_idx[i,]] <- replace(snp_matrix[,col_idx[i,]], which(snp_matrix[,col_idx[i,]] == 1),base[1])
  snp_matrix[,col_idx[i,]] <- replace(snp_matrix[,col_idx[i,]], which(snp_matrix[,col_idx[i,]] == 2),base[2])
}

fid <- rep("pop1", nind)
iid <- paste0("iid", 1:nind)
pid <- rep(0, nind)
mid <- pid
sex <- sample(1:2, nind, replace = T)
phy <- sample(1:2, nind, replace = T)

resultped <- cbind(fid, iid, pid, mid, sex, phy, snp_matrix)

chr <- rep(1, nsnp)
dis <- rep(0, nsnp)
pos <- sort(sample(1:200000000, nsnp))
snpid <- paste0(chr,":", pos)

resultmap <- cbind(chr, snpid, dis, pos)

write.table(resultped, "result.ped", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(resultmap, "result.map", row.names = F, col.names = F, quote = F, sep = "\t")