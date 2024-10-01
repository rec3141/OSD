cluster.dat <- read.csv("cluster_all_sort.csv",sep="\t",header=F)
cluster.dat <- cluster.dat[,-11]
colnames(cluster.dat) <- c("cluster","refseq","id","length","arg_lys","acidic","aliphaticity","aliphatic_index","proline","gravy")

attach(cluster.dat)

plot(cluster.dat[sample(1:nrow(cluster.dat),10000),5:9],pch=16,cex=0.3,col=rgb(1, 0, 0, 0.1))

boxplot(arg_lys ~ cluster, subset = cluster < 40, data=cluster.dat)
boxplot(arg_lys ~ cluster, subset = refseq == "NC_009997", data=cluster.dat)
boxplot(arg_lys ~ cluster, subset = refseq == "NC_009997" & cluster < 20, data=cluster.dat)

###for export to itol

med.arg_lys <- aggregate(arg_lys ~ refseq, data = cluster.dat, median)
write.csv(file="med.arg_lys.csv",c(med.arg_lys,,row.names=F,quote=F)

med.acidic <- aggregate(acidic ~ refseq, data = cluster.dat, median)
write.csv(file="med.acidic.csv",med.acidic,row.names=F,quote=F)

med.aliphaticity <- aggregate(aliphaticity ~ refseq, data = cluster.dat, median)
write.csv(file="med.aliphaticity.csv",med.aliphaticity,row.names=F,quote=F)

med.aliphatic_index <- aggregate(aliphatic_index ~ refseq, data = cluster.dat, median)
write.csv(file="med.aliphatic_index.csv",med.aliphatic_index,row.names=F,quote=F)

med.proline <- aggregate(proline ~ refseq, data = cluster.dat, median)
write.csv(file="med.proline.csv",med.proline,row.names=F,quote=F)

med.gravy <- aggregate(gravy ~ refseq, data = cluster.dat, median)
med.gravy[,2] <- med.gravy[,2] + 1
write.csv(file="med.gravy.csv",med.gravy,row.names=F,quote=F)


cluster.split <- split(cluster.dat, refseq)


# for untransformed
for(j in c("acidic","aliphaticity","aliphatic_index","proline","gravy")) {

 hist.save <- matrix(NA,nrow=length(cluster.split),ncol=100)
 rownames(hist.save) <- names(cluster.split)
 breaks <- seq( min(cluster.dat[,j],na.rm=T), max(cluster.dat[,j],na.rm=T), length.out=101 )

 for(i in 1:length(cluster.split)) {
	 hist.save[i,] <- hist(cluster.split[[i]][,j],breaks=breaks, plot=F)$density

 }
image(hist.save)
print(j)
fivenum(hist.save)
write.csv(file=paste("itol",j,"csv",sep="."),hist.save[,1:30],quote=FALSE)

}

#> fivenum(proline)
#[1]   0.00000000   0.06566224   2.68209832  12.45210728 224.13793103



#for log-transformed
for(j in c("arg_lys")) {

 hist.save <- matrix(NA,nrow=length(cluster.split),ncol=100)
 rownames(hist.save) <- names(cluster.split)

 breaks <- seq( -7, 7, length.out=101)

 for(i in 1:length(cluster.split)) {
	 hist.save[i,] <- hist(log(cluster.split[[i]][,j]),breaks=breaks, plot=F)$density

 }

 image(hist.save)
 write.csv(file=paste("itol",j,"csv",sep="."),hist.save[,30:70],quote=FALSE)
}



# for each cluster
list.psy <- read.csv("blastlist-hima-genomes",header=F)
list.psy <- cbind(list.psy,rep(1,length(list.psy)))
colnames(list.psy) <- c("refseq","psychro")
cluster.add <- merge(cluster.dat,list.psy,all.x=T)
cluster.add$psychro[is.na(cluster.add$psychro)] <- 0
cluster.add$arg_lys <- log(cluster.add$arg_lys)
cluster.add$arg_lys[is.infinite(cluster.add$arg_lys)] <- NA


prot.vars <- colnames(cluster.add)[5:10]
colnames(ttest.save) <- prot.vars

for(j in 1:length(prot.vars)) {
 jc <- prot.vars[j]
 ttest.save <- matrix(NA,nrow=9199,ncol=11)
# hist.save.ref <- matrix(NA,nrow=9199,ncol=100)
# hist.save.psy <- matrix(NA,nrow=9199,ncol=100)
# breaks <- seq( min(cluster.add[,jc],na.rm=T), max(cluster.add[,jc],na.rm=T), length.out=101 )
 print(jc)

 for(i in 1:9199) {
 	 cluster.ref <- cluster.add[ cluster.add$cluster == i & cluster.add$psychro == 0, jc]
 	 cluster.psy <- cluster.add[ cluster.add$cluster == i & cluster.add$psychro == 1, jc]
	 if (i %% 1000 == 0) {print(i)}
	 if (length(na.omit(cluster.ref)) > 1 & length(na.omit(cluster.psy)) > 1) {
		 try(
		 	ttest.run <- t.test(cluster.ref,cluster.psy)
		 	)
		 	ttest.save[i,] <- unlist(ttest.run)
	 }

#	 hist.ref <- hist(cluster.ref,breaks=breaks, plot=T, main=paste(i,jc,ttest.save[i],sep="."))
#	 hist.save.ref[i,] <- hist.ref$density
#	 hist.psy <- hist(cluster.psy,breaks=breaks, plot=T, add=T, col="blue")
#	 hist.save.psy[i,] <- hist.psy$density

#	 if (ttest.save[i] < 0.05/max(cluster.add$cluster)) { readline() }
 }

write.csv(file=paste("ttest.save",jc,"csv",sep="."),ttest.save,quote=FALSE)
 
#image(hist.save,main=jc)
#fivenum(hist.save)
#write.csv(file=paste("ttest.save",jc,"csv",sep="."),ttest.save,quote=FALSE)

}




ttest.my <- ttest.save[1:6126,]
ttest.sort <- sort(apply(ttest.my,1,min),index=T)
ttest.ord <- ttest.my[ ttest.sort$ix, ] 
rownames(ttest.ord) <- ttest.sort$ix


#print out
list.interest <- c(224,216,185,151,11,217,142,224,201,214,33,175,164,140,212,176,107,217,3,114,34,26,136,212,217,94,152,45,158,17,215,221,165,123,147,190,182,192,184,179,194,38,211,51,194,122,188,164,156,132,124,174,128,161,163,94,119,15,223,206,224,128,155,131,197,161,33,17,5,77,213,5,143,185,118,122,121,178,220,181,221,185,146,209,123,96,200,220,94,52,155,128,133,197,8,224,77,159,224,111,113,147,199,177,175,136,120,195,194,209,16,186,89,138,111,110,203,102,170,169,209,90,32,202,158,55,158,61,90,210,204,143,75,139,168,183,166,7,128,190,95,116,22,130,118,18,76,187,182,62,177,67,170,190,163)
for(i in list.interest) {
	for(j in 1:length(prot.vars)) {
	 	jc <- prot.vars[j]
	    breaks <- seq( min(cluster.add[,jc],na.rm=T), max(cluster.add[,jc],na.rm=T), length.out=101 )
	 	cluster.ref <- cluster.add[ cluster.add$cluster == i & cluster.add$psychro == 0, jc]
 	 	cluster.psy <- cluster.add[ cluster.add$cluster == i & cluster.add$psychro == 1, jc]
		png(file=paste('aminos',jc,i,'png',sep='.'),width=640,height=640)
		hist.ref <- hist(cluster.ref,breaks=breaks, plot=T, main=paste(i,jc,ttest.save[i],sep="."))
	 	hist.psy <- hist(cluster.psy,breaks=breaks, plot=T, add=T, col="blue")
		dev.off()
	}
}




