#[1] "statistic.t"                    "parameter.df"                  
# [3] "p.value"                        "conf.int1"                     
# [5] "conf.int2"                      "estimate.mean of x"            
# [7] "estimate.mean of y"             "null.value.difference in means"
# [9] "alternative"                    "method"                        
#[11] "data.name"                     

#statistic.t
#parameter.df
#p.value
#conf.int1
#conf.int2
#estimate.mean of x
#estimate.mean of y
#null.value.difference in means
#alternative
#method
#data.name

#psych - meso

ttest.pval.all <- NULL

ttest.bonf <- 0.05/1e4


##### lumpy heart

#arg_lys
ttest.save.tmp <- read.csv("ttest.save.arg_lys.csv",row.names=1)
ttest.pval.all <- cbind(ttest.pval.all, ttest.save.tmp[,3])
colnames(ttest.save.tmp) <- names(unlist(t.test(c(0,1),c(1,2))))

ttest.diffs <- ttest.save.tmp[,7] - ttest.save.tmp[,6]
ttest.pval <- -1*log(ttest.save.tmp[,3])

hist(ttest.diffs,breaks=101)
# long left tail
mean(ttest.diffs,na.rm=T)
# [1] -0.2746232
sd(ttest.diffs,na.rm=T)
# [1] 0.4136419
pdf(file="heart.arg_lys.pdf",width=6,height=6)
plot(ttest.pval ~ ttest.diffs,log='y',pch=16,col=rgb(1,.0,.0,.05),ylab="p-value",xlab="difference in means",main="Arg/Lys")
points(ttest.pval[ttest.group3] ~ ttest.diffs[ttest.group3],pch=16,col=rgb(0,0,1,.3))
#lines(c(-3,3), c(-log(ttest.bonf),-log(ttest.bonf)))
#plot(ttest.diffs,na.omit(ttest.save.tmp[,1]),pch=16,col=rgb(0,0,0,0.1))
dev.off()


#aliphaticity
ttest.save.tmp <- read.csv("ttest.save.aliphaticity.csv",row.names=1)
ttest.pval.all <- cbind(ttest.pval.all, ttest.save.tmp[,3])
ttest.diffs <- ttest.save.tmp[,7] - ttest.save.tmp[,6]
ttest.pval <- -1*log(ttest.save.tmp[,3])
hist(ttest.diffs,breaks=101)
mean(ttest.diffs,na.rm=T)
# [1] -0.009030072
sd(ttest.diffs,na.rm=T)
# [1] 0.0198003
pdf(file="heart.aliphaticity.pdf",width=6,height=6)
plot(ttest.pval ~ ttest.diffs,log='y',pch=16,col=rgb(1,.0,.0,.05),ylab="p-value",xlab="difference in means",main="Aliphaticity")
points(ttest.pval[ttest.group3] ~ ttest.diffs[ttest.group3],pch=16,col=rgb(0,0,1,.3))
#lines(c(-3,3), c(-log(ttest.bonf),-log(ttest.bonf)))
dev.off()

#proline
ttest.save.tmp <- read.csv("ttest.save.proline.csv",row.names=1)
ttest.pval.all <- cbind(ttest.pval.all, ttest.save.tmp[,3])
ttest.diffs <- ttest.save.tmp[,7] - ttest.save.tmp[,6]
ttest.pval <- -1*log(ttest.save.tmp[,3])
hist(ttest.diffs,breaks=101)# long left tail
mean(ttest.diffs,na.rm=T)
# [1] -0.00333705
sd(ttest.diffs,na.rm=T)
# [1] 0.006308175
pdf(file="heart.proline.pdf",width=6,height=6)
plot(ttest.pval ~ ttest.diffs,log='y',pch=16,col=rgb(1,.0,.0,.05),ylab="p-value",xlab="difference in means",main="Proline")
points(ttest.pval[ttest.group3] ~ ttest.diffs[ttest.group3],pch=16,col=rgb(0,0,1,.3))
#lines(c(-3,3), c(-log(ttest.bonf),-log(ttest.bonf)))
dev.off()

####### normal heart


#acidic
ttest.save.tmp <- read.csv("ttest.save.acidic.csv",row.names=1)
ttest.pval.all <- cbind(ttest.pval.all, ttest.save.tmp[,3])
ttest.diffs <- ttest.save.tmp[,7] - ttest.save.tmp[,6]
ttest.pval <- -1*log(ttest.save.tmp[,3])
hist(ttest.diffs,breaks=101)
mean(ttest.diffs,na.rm=T)
# [1] 0.00149369
sd(ttest.diffs,na.rm=T)
# [1] 0.01074217
pdf(file="heart.acidic.pdf",width=6,height=6)
plot(ttest.pval ~ ttest.diffs,log='y',pch=16,col=rgb(1,.0,.0,.05),ylab="p-value",xlab="difference in means",main="Acidic")
points(ttest.pval[ttest.group3] ~ ttest.diffs[ttest.group3],pch=16,col=rgb(0,0,1,.3))
#lines(c(-3,3), c(-log(ttest.bonf),-log(ttest.bonf)))
dev.off()

#aliphatic_index
ttest.save.tmp <- read.csv("ttest.save.aliphatic_index.csv",row.names=1)
ttest.pval.all <- cbind(ttest.pval.all, ttest.save.tmp[,3])
ttest.diffs <- ttest.save.tmp[,7] - ttest.save.tmp[,6]
ttest.pval <- -1*log(ttest.save.tmp[,3])
hist(ttest.diffs,breaks=101)# long left tail
mean(ttest.diffs,na.rm=T)
# [1] 0.008941378
sd(ttest.diffs,na.rm=T)
# [1] 0.0396092
pdf(file="heart.aliphatic_index.pdf",width=6,height=6)
plot(ttest.pval ~ ttest.diffs,log='y',pch=16,col=rgb(1,.0,.0,.05),ylab="p-value",xlab="difference in means",main="Aliphatic Index")
points(ttest.pval[ttest.group3] ~ ttest.diffs[ttest.group3],pch=16,col=rgb(0,0,1,.3))
#lines(c(-3,3), c(-log(ttest.bonf),-log(ttest.bonf)))
dev.off()


#gravy
ttest.save.tmp <- read.csv("ttest.save.gravy.csv",row.names=1)
ttest.pval.all <- cbind(ttest.pval.all, ttest.save.tmp[,3])
ttest.diffs <- ttest.save.tmp[,7] - ttest.save.tmp[,6]
ttest.pval <- -1*log(ttest.save.tmp[,3])
hist(ttest.diffs,breaks=101)# bimodal, long left tail
mean(ttest.diffs,na.rm=T)
# [1] 0.009376566
sd(ttest.diffs,na.rm=T)
# [1] 0.0756709
pdf(file="heart.gravy.pdf",width=6,height=6)
plot(ttest.pval ~ ttest.diffs,log='y',pch=16,col=rgb(1,.0,.0,.05),ylab="p-value",xlab="difference in means",main="GRAVY")
points(ttest.pval[ttest.group3] ~ ttest.diffs[ttest.group3],pch=16,col=rgb(0,0,1,.3))
#lines(c(-3,3), c(-log(ttest.bonf),-log(ttest.bonf)))
dev.off()

ttest.pval.max <- apply(ttest.pval.all,1,max)
plot(sort(-log(ttest.pval.max)))

ttest.group1 <- which(ttest.pval.all[,1] < 5e-6 & ttest.pval.all[,2] < 5e-6 & ttest.pval.all[,3] < 5e-6)
ttest.group2 <- which(ttest.pval.all[,1] < 5e-4 & ttest.pval.all[,2] < 5e-4 & ttest.pval.all[,3] < 5e-4)
ttest.group3 <- which(ttest.pval.all[,1] < 5e-2 & ttest.pval.all[,2] < 5e-2 & ttest.pval.all[,3] < 5e-2)


#####################

ttest.proline <- t.test(cluster.add[cluster.add$psy == 0,"proline"], cluster.add[cluster.add$psy == 1,"proline"])
ttest.proline$estimate[2] - ttest.proline$estimate[1] # -0.005511948
#proline: expect fewer; found fewer

ttest.arg_lys <- t.test(cluster.add[cluster.add$psy == 0,"arg_lys"], cluster.add[cluster.add$psy == 1,"arg_lys"])
ttest.arg_lys$estimate[2] - ttest.arg_lys$estimate[1] # -0.4510342
#arg_lys: expect less arg; found less

ttest.aliphaticity <- t.test(cluster.add[cluster.add$psy == 0,"aliphaticity"], cluster.add[cluster.add$psy == 1,"aliphaticity"])
ttest.aliphaticity$estimate[2] - ttest.aliphaticity$estimate[1] # -0.01835353
#aliphatic: expect less; found slightly less (aliphaticity) ; 

###

ttest.acidic <- t.test(cluster.add[cluster.add$psy == 0,"acidic"], cluster.add[cluster.add$psy == 1,"acidic"])
ttest.acidic$estimate[2] - ttest.acidic$estimate[1] # 0.004132139
#acidic: expect fewer; found very slightly more ; 2nd struc?

ttest.aliphatic_index <- t.test(cluster.add[cluster.add$psy == 0,"aliphatic_index"], cluster.add[cluster.add$psy == 1,"aliphatic_index"])
ttest.aliphatic_index$estimate[2] - ttest.aliphatic_index$estimate[1] # 0.0133294
#aliphatic: expect less; found very slightly higher (aliphatic_index) ; murray found same, 2nd struc?

ttest.gravy <- t.test(cluster.add[cluster.add$psy == 0,"gravy"], cluster.add[cluster.add$psy == 1,"gravy"])
ttest.gravy$estimate[2] - ttest.gravy$estimate[1] # 0.00221651
#gravy: expect lower; found very slightly higher ; murray found same, 2nd struc?

hydrophobic: expect more (?);


############################

mp <- merge(med.proline,list.psy,all.x=T)
mp[is.na(mp)] <- 0
mp <- mp[ order(mp$proline) ,]
plot(mp$proline,pch=16,col=rgb(0,0,mp$psy,1),cex=2)

mp <- merge(med.arg_lys,list.psy,all.x=T)
mp[is.na(mp)] <- 0
mp <- mp[ order(mp$arg_lys) ,]
plot(mp$arg_lys,pch=16,col=rgb(0,0,mp$psy,1),cex=2)

mp <- merge(med.aliphaticity,list.psy,all.x=T)
mp[is.na(mp)] <- 0
mp <- mp[ order(mp$aliphaticity) ,]
plot(mp$aliphaticity,pch=16,col=rgb(0,0,mp$psy,1),cex=2)

mp <- merge(med.aliphatic_index,list.psy,all.x=T)
mp[is.na(mp)] <- 0
mp <- mp[ order(mp$aliphatic_index) ,]
plot(mp$aliphatic_index,pch=16,col=rgb(0,0,mp$psy,1),cex=2)

mp <- merge(med.acidic,list.psy,all.x=T)
mp[is.na(mp)] <- 0
mp <- mp[ order(mp$acidic) ,]
plot(mp$acidic,pch=16,col=rgb(0,0,mp$psy,1),cex=2)

mp <- merge(med.gravy,list.psy,all.x=T)
mp[is.na(mp)] <- 0
mp <- mp[ order(mp$gravy) ,]
plot(mp$gravy,pch=16,col=rgb(0,0,mp$psy,1),cex=2)

########


# compare paralogs that are differentially expressed at different temperatures
# to find cold-core and warm-core

# they don't correct for multiple comparisons


disordered region predictions with API:
http://bioinf.cs.ucl.ac.uk/web_servers/web_services/
http://bip.weizmann.ac.il/fldbin/findex

with source:
http://globplot.embl.de/


#Metpally and
#Reddy (41) compared amino acid substitution patterns for
#psychrophilic and mesophilic organisms, and the psychrophilic
#organisms used included Psychrobacter cryohalolentis K5, C.
#psychrerythraea 34H, and P. ingrahamii 37. They concluded that
#psychrophile proteins contained fewer hydrophilic, acidic, and
#proline residues, consistent with our findings (41, 69)

#Between 476 (38%) and 1,074 (84%) of the P.
#arcticus amino acid sequences were statistically significantly
#different from sequences in the Swiss-Prot database for all cold
#adaptation indicators; i.e., they were less hydrophobic, had
#fewer proline residues, were less aliphatic, had fewer acidic
#residues, or had low Arg and increased Lys contents (Fig. 2).

One of the biggest challenges for proteins at low temperatures
is having sufficient flexibility so that they can increase
their interactions with substrates, reducing the required activation
energy. An abundance of proline residues has been
related to increased protein stability due to the rigidity of the
N-C
 bond (19, 51). Hence, the decrease in the amount of
proline observed suggests that there is cold adaptation, which
supports a trend seen in smaller-scale studies (69). Arginine is
also considered an amino acid that stabilizes the structurally
since it forms salt bridges and hydrogen bonds with side chains
(1). Substituting lysine for arginine has been proposed to be a
substitution that results in more flexibility. The negatively
charged (acidic) amino acids glutamic acid and aspartic acid
favor salt bridge formation on protein surfaces, thus favoring a
stable protein structure (23). Removal of acidic residues increases
protein flexibility (20). Finally, an increase in the hydrophobicity
of core amino acids increases protein stability at
higher temperatures (64), while an overall reduction in stability
has been observed in cold-active enzymes (63). Our results
suggest that the adaptation of P. arcticus to low temperatures
involved multiple amino acid substitutions that decreased protein
stability, presumably yielding enzymes that are more active
at low temperatures. The only parameter not consistent with
this conclusion was aliphacity since it showed the opposite
response to cold. This inconsistent result, which has been observed
previously, could be due to the strategy used for analysis,
which did not separate exposed and buried residues (23,
41, 69). In support of this, Metpally and Reddy (41) recently
showed that there is a high frequency of aliphatic residues in
coil-loop regions of proteins encoded by psychrophile genomes.

Acidic amino acids group include D and E; aliphatic: I, L
and V; aromatic: H, F, W and Y; basic: R, H, and K;
charged: R, D, E, H and K; hydrophilic: D, E, K, N, Q and
R; hydrophobic: A, C, F, I, L, M, V, W and Y neutral: G, Q,
H, S and T; non-polar: A, C, G, I, L, M, F, P, V, W and Y;
polar: R, N, D, E, Q, H, K, S and T; small: A, C, D, G, N, P,
S, T, V and tiny: A, C, G, S and T. 

