library(deepKin)

plink_path = "~/Desktop/project_management/tools/plink2"
gear_path = "~/Desktop/project_management/tools/gear.jar"
bfileprefix = "~/Desktop/project_management/GWASrelatedness/data/euro_10Ksnp"

## deepKin Principle I:  
n = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".fam"), intern = T))
m = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".bim"), intern = T))
id = read.table(file = paste0(bfileprefix,".fam"), header = F)[,1]
npairs = n*(n-1)/2
me.min = me.min(theta = (1/2)^(0:4), alpha = 0.05/npairs, beta = 0.1)

## GRM method or RDM method
me = calculate.me(bfileprefix = bfileprefix, method = "GRM", plink_path = plink_path, pop_size = n)  # not suggested for biobank data
me = calculate.me(bfileprefix = bfileprefix, method = "RDM", gear_path = gear_path, pop_size = n)     # suggested for biobank data



## deepKin Principle II:
theta.min = theta.min(me = me, alpha = 0.05/npairs, beta = 0.1)
degree.deep = log(theta.min, base = 1/2)

## deepKin Principle III:
power.max = power.max(me, theta = (1/2)^(0:6), alpha = 0.05/npairs)




## Perform plink GRM
if(!file.exists(paste0(bfileprefix,".rel.bin"))){
  system(paste0(plink_path, " --silent --bfile ", bfileprefix, " --make-rel triangle bin4 --out ", bfileprefix))
}
## Extract plink GRM
grm.rst = extract.plink.grm(bfileprefix, xcohort = F, pop_size1 = n)
grm.diag = grm.rst$diag
grm.tri  = grm.rst$tri

## Perform deepKin
deepkin = deepKin.estimation(grm.diag = grm.diag, grm.tri = grm.tri, xcohort = F, me = me)

## Output estimation and p-values (This file can be large according to the number of pairs)
write.table(deepkin, file = paste0(bfileprefix, ".deepkin"), quote = F, col.names = T, row.names = F)

thrd = (1/2)^seq(0.5, 10.5, 1)
thrd = thrd[which(thrd>theta.min)]
deepkin.qc = deepkin[which(deepkin$king > theta.min), ]
deepkin.qc$degree = cut(deepkin.qc$king,
                        breaks = c(theta.min,thrd,1),
                        labels = c(paste0("Deepest-",length(thrd)-1), (length(thrd)-1):0))
index = as.numeric(rownames(deepkin.qc))
deepkin.qc.id = extract.indi.id(id = id, index = index, xcohort = F)
deepkin.qc.rst = cbind(deepkin.qc.id, deepkin.qc)

write.table(deepkin.qc.rst, file = paste0(bfileprefix, ".related"), quote = F, col.names = T, row.names = F)

