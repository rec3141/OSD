#!/usr/bin/R
rm(list=ls())

library(phytools)
library(geiger)
library(genoPlotR)
library(reshape2)

#pseudocode:

#for each cluster
#  for each genome
#     start at a leaf
#     sort remaining leaves by distance from this leaf
#     for each leaf
#       check for cluster
#       calculate percentage YES
#       next leaf
#      plot distance versus percentage

#inputs: 

#1. taxalist
# added these manually to taxalist-cellular renamed taxalist-cellular-rec
# NZ changed to NCs but not sure why NCs are missing
# NC_009784 Vcho ?
# NC_009457 Vcho ?
# NZ_AAPG000 NC_020802
# NZ_ABSH000 NC_020911
# NZ_ABSK000 NC_020908
# NC_015844 Zgal ?
# NZ_AAPR000 NC_018721
# NZ_AANA000 NC_020830
# NC_014738 Rana ?

#NZ_AAPG00000000
#NZ_ABSH00000000
#NZ_ABSK00000000
#NZ_AAPR00000000
#NZ_AANA00000000

taxalistfile <- "taxalist-cellular-rec"
taxalist.read <- scan(taxalistfile,sep="\n","char")

#2. *.tre -- newick tree file with same genomes as cluster_list
trefile <- "HI_110003.gb.concat.phylip_phyml_tree.tre"
if(!file.exists(trefile)) stop(paste(trefile,"does not exist",sep=" "))
tre.input <- read.tree(trefile)

taxalist.tree <- lapply(tre.input$tip.label, function(x) grep(x,taxalist.read,value=T))
taxalist.tree <- as.data.frame(cbind(tre.input$tip.label,substr(taxalist.tree,1,9)))
colnames(taxalist.tree) <- c("refseq","univ")

tre.input$tip.label <- as.character(taxalist.tree$univ)

#3. *.dist.pair.csv -- distance matrix of distances between genomes, output from tre2table.pl
#should be in same order as tree
distfile <- "list-0.dist.pair.csv"
if(!file.exists(distfile)) stop(paste(distfile,"does not exist",sep=" "))
dist.matrix.read <- read.csv(distfile,header=T,row.names=1,sep="\t")

colnames(dist.matrix.read) <- taxalist.tree$univ
rownames(dist.matrix.read) <- taxalist.tree$univ




## NOT NEEDED HERE SINCE I HAVE med.*.csv NOW

#3. *.cluster_list.r.csv -- list of protein clusters

#clusterfile <- "HI_110001.1e-20_.8_2_1.cluster_list.ren.csv"
#if(!file.exists(clusterfile)) stop(paste(clusterfile,"does not exist",sep=" "))
#cluster.read<-read.csv(clusterfile,sep="\t", header=T)

#4. cold_info.pl output
#infofile <- "cluster_all_sort.csv"
#info.read <- read.csv(infofile,sep="\t")
#info.read <- info.read[,-11]
#colnames(info.read) <- c("cluster","refseq","id","length","arg_lys","acidic","aliphaticity","aliphatic_index","proline","gravy")
#attach(cluster.dat)

#taxalist.info <- lapply(unique(info.read$refseq), function(x) grep(x,taxalist.read,value=T))
#taxalist.info <- cbind(as.character(unique(info.read$refseq)),substr(taxalist.info,1,9))
#colnames(taxalist.info) <- c("refseq","univ")

#info.dat <- merge(info.read,taxalist.info)

#then aggregate by univ

#plot(info.dat[sample(1:nrow(info.dat),10000),5:9],pch=16,cex=0.3,col=rgb(1, 0, 0, 0.1))

#boxplot(arg_lys ~ cluster, subset = cluster < 40, data=info.dat)
#boxplot(arg_lys ~ cluster, subset = refseq == "NC_009997", data=info.dat)
#boxplot(arg_lys ~ cluster, subset = refseq == "NC_009997" & cluster < 20, data=info.dat)

###for export to itol
#info.melt <- melt(info.dat,id=c("univ","cluster","refseq","id"))

##med.arg_lys <- aggregate(arg_lys ~ univ + cluster, data = info.dat, median)
#cast.arg_lys <- cast(info.melt, cluster ~ univ, fun.aggregate=median, fill=NA, add.missing=T, subset=variable=="arg_lys")
#write.csv(file="med.arg_lys.csv",cast.arg_lys,row.names=F,quote=F)

##med.acidic <- aggregate(acidic ~ univ, data = info.dat, median)
#cast.acidic <- cast(info.melt, cluster ~ univ, fun.aggregate=median, fill=NA, add.missing=T, subset=variable=="acidic")
#write.csv(file="med.acidic.csv",cast.acidic,row.names=F,quote=F)

##med.aliphaticity <- aggregate(aliphaticity ~ univ, data = info.dat, median)
#cast.aliphaticity <- cast(info.melt, cluster ~ univ, fun.aggregate=median, fill=NA, add.missing=T, subset=variable=="aliphaticity")
#write.csv(file="med.aliphaticity.csv",cast.aliphaticity,row.names=F,quote=F)

##med.aliphatic_index <- aggregate(aliphatic_index ~ univ, data = info.dat, median)
#cast.aliphatic_index <- cast(info.melt, cluster ~ univ, fun.aggregate=median, fill=NA, add.missing=T, subset=variable=="aliphatic_index")
#write.csv(file="med.aliphatic_index.csv",cast.aliphatic_index,row.names=F,quote=F)

##med.proline <- aggregate(proline ~ univ, data = info.dat, median)
#cast.proline <- cast(info.melt, cluster ~ univ, fun.aggregate=median, fill=NA, add.missing=T, subset=variable=="proline")
#write.csv(file="med.proline.csv",cast.proline,row.names=F,quote=F)

##med.gravy <- aggregate(gravy ~ univ, data = info.dat, median)
##med.gravy[,2] <- med.gravy[,2] + 1  ## FORGOT WHY...
#cast.gravy <- cast(info.melt, cluster ~ univ, fun.aggregate=median, fill=NA, add.missing=T, subset=variable=="gravy")
#write.csv(file="med.gravy.csv",cast.gravy,row.names=F,quote=F)



##choose which data to measure SURPRISE

# rearrange data to match the tree

mat.arg_lys <- read.csv(file="med.arg_lys.csv",row.names=1)
#mat.acidic <- read.csv(file="med.acidic.csv",row.names=1)
#mat.aliphaticity <- read.csv(file="med.aliphaticity.csv",row.names=1)
#mat.aliphatic_index <- read.csv(file="med.aliphatic_index.csv",row.names=1)
#mat.proline <- read.csv(file="med.proline.csv",row.names=1)
#mat.gravy <- read.csv(file="med.gravy.csv",row.names=1)

mat.read <- mat.gravy # to read in arg_lys values

mat <- mat.read[,taxalist.tree$univ]
# trim the tree to account for available data
diff.mat <- name.check(tre.input,t(mat[1,]))
if (diff.mat == "OK") {
	tre.trim <- tre.input
} else {
	tre.trim <- drop.tip(tre.input,diff.mat$tree_not_data)
}


# rearrange the dist matrix to match tree
dist.matrix <- dist.matrix.read[tre.trim$tip.label,tre.trim$tip.label]
max.dist <- max(rowSums(dist.matrix))
sim.matrix <- max(dist.matrix) - dist.matrix


######

testnum <- dim(mat)[1]
testnum <- 2000
startnum <- 1
mat.mini <- t(mat[startnum:testnum,tre.trim$tip.label])
#mat.mini[is.na(mat.mini)] <- 0

m.surp <- matrix(data=NA,nrow=nrow(mat.mini),ncol=ncol(mat.mini))
colnames(m.surp) <- colnames(mat.mini)
rownames(m.surp) <- rownames(mat.mini)

#surpfun <- function(x) sum((n.leaves[x] - n.leaves) * as.double(gen.sims))

for (t.cluster in colnames(mat.mini)) {
	print(t.cluster)
	cluster = mat.mini[,t.cluster]
#	par(mfrow=c(1,2), ask=T)
#	par(mar=c(4.5, 4.5, 4.5, 4.5))
#	plot(c(-1,-1),c(-1,-1),xlim=c(0,1.5),ylim=c(-4,4),main=mat.read[t.cluster,"description"]);
#	plot(c(-1,-1),c(-1,-1),xlim=c(0,5),ylim=c(-10,10),main=t.cluster);
	for (t.genome in names(cluster)) {
		gen.dists <- sort(dist.matrix[t.genome,])
		gen.sims <- (max(gen.dists,na.rm=T) - gen.dists)/max(gen.dists,na.rm=T)
		n.leaves <- cluster[names(gen.dists)]
		m.surp[t.genome,t.cluster] <- sum((n.leaves[t.genome] - n.leaves) * as.double(gen.sims),na.rm=T)
#		points((n.leaves[t.genome] - n.leaves) * as.double(gen.sims) ~ as.numeric(gen.dists))
	}
#	plotTree.wBars(tre.trim,x=mat.mini[,t.cluster]/max(mat.mini[,t.cluster]+1,na.rm=T),scale=0.5)
#	plotTree.wBars(tre.trim,x=m.surp[,t.cluster]/max(m.surp[,t.cluster]+1,na.rm=T),scale=0.5)

#	SIMMAP for stochastic mapping
#	tre.trim.simmap <- make.simmap(tre.trim, cluster, model="SYM", nsim=1)
#	max.states <- as.integer(max(rownames(tre.trim.simmap$Q)))
#	sim.colors <- setNames(colorRampPalette(c('blue','red'))(max.states+1),0:max.states)
#	plotTree.wBars(tre.trim.simmap,x=mat.mini[,t.cluster]/max(mat.mini[,t.cluster]+1),scale=0.5,method="plotSimmap",colors=sim.colors,lwd=3,node.numbers=T)
#	plotTree.wBars(tre.trim.simmap,x=m.surp[,t.cluster]/max(m.surp[,t.cluster]+1),scale=0.5,method="plotSimmap",colors=sim.colors,lwd=3,node.numbers=T)	
	}
write.csv(t(m.surp),"gravy.surprise.csv")
	
stop()





