setwd("~/Desktop/project_management/DeepKin_verify/simulation")
source("~/Desktop/project_management/DeepKin_verify/source.R")


set.seed(123456)
n1 = 1000                       # sample size for pop 1
n2 = 1000                       # sample size for pop 2
Mvec = seq(1000, 2000, 250)     # number of markers
Mnum = length(Mvec)
repeats = 10
r = 1

a = 0
dtp =list()
for(i in 1:Mnum)
{
  for(j in 1:repeats)
  {
    a = a + 1
    M = Mvec[i]
    freq = runif(M, 0.25, 0.5)
    ## Parent offsprings ----
    Gr = GenerateGeno_r(freq, n1, r, sibflag = F)
    X1 = Gr[[1]]
    X2 = Gr[[2]]
    fhat = colMeans(rbind(X1,X2))/2
    X1 = t(apply(X1, 1, function(x) {(x-2*fhat)/sqrt(2*fhat*(1-fhat))}))
    X2 = t(apply(X2, 1, function(x) {(x-2*fhat)/sqrt(2*fhat*(1-fhat))}))
    KING.po = 1-rowSums((X1-X2)^2)/2/M
    
    
    ## Full sibs ----
    Gr = GenerateGeno_r(freq, n1, r, sibflag = T)
    X1 = Gr[[1]]
    X2 = Gr[[2]]
    fhat = colMeans(rbind(X1,X2))/2
    X1 = t(apply(X1, 1, function(x) {(x-2*fhat)/sqrt(2*fhat*(1-fhat))}))
    X2 = t(apply(X2, 1, function(x) {(x-2*fhat)/sqrt(2*fhat*(1-fhat))}))
    KING.fs = 1-rowSums((X1-X2)^2)/2/M
    

    dtp[[a]] = data.frame(grp = c("Theoretical", 
                                  "Observed (Parent-offsping)", 
                                  "Observed (Full siblings)"),
                          M = M, rep = j,
                          var = c(2*(1-(0.5)^r)^2/M, var(KING.po), var(KING.fs)))   
    print(a)
  }
}

library(dplyr)
library(ggplot2)
library(ggsci)
library(egg)

# bind rows
dt = dplyr::bind_rows(dtp)
# save r data
save(dt, file = "oneChr.Rdata")
load("oneChr.Rdata")
# calculate mean and sd
dt.sm = dt %>% group_by(M, grp) %>%
  summarise(mean = mean(var), sd = sd(var))
dt.sm$M = as.factor(dt.sm$M)
dt.sm$grp = factor(dt.sm$grp, levels = c("Theoretical", 
                                         "Observed (Parent-offsping)", 
                                         "Observed (Full siblings)"))

df = dt.sm %>% dplyr::mutate(sd = ifelse(sd==0, NA, sd))

pdf(file = "oneChr.pdf", width = 7, height = 5)
ggplot(data = df, aes(x=M, y=mean, color=grp, fill=grp))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin = mean-1.96*sd, ymax = mean+1.96*sd), na.rm = T, color = "black",
                stat="identity", width=.2, position=position_dodge(.9)) +
  scale_color_npg(name = "") +
  scale_fill_npg(name = "", alpha = 0.6) +
  labs(x = "m", y="Variance")+
  theme_article()+
  theme(plot.title = element_text(hjust = 0.5, size=9),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())  
dev.off()
