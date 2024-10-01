#!/usr/bin/R
rm(list=ls())

library(phytools)
library(geiger)
library(genoPlotR)

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
#1. *.cluster_list.r.csv -- list of protein clusters
clusterfile <- "HI_110001.1e-20_.8_2_1.cluster_list.ren.csv"
if(!file.exists(clusterfile)) stop(paste(clusterfile,"does not exist",sep=" "))
mat.read<-read.csv(clusterfile,sep="\t", header=T)

#2. *.tre -- newick tree file with same genomes as cluster_list
trefile <- "HI_110003.gb.concat.phylip_phyml_tree.tre"
if(!file.exists(trefile)) stop(paste(trefile,"does not exist",sep=" "))
tre.input <- read.tree(trefile)

#3. *.dist.pair.csv -- distance matrix of distances between genomes, output from tre2table.pl
distfile <- "list-0.dist.pair.csv"
if(!file.exists(distfile)) stop(paste(distfile,"does not exist",sep=" "))
dist.matrix.read <- read.csv(paste(distfile,sep=""),header=T,row.names=1,sep="\t")


# input protein cluster list
mat.df <- mat.read[,4:(ncol(mat.read)-4)]
mat <- as.matrix(mat.df)
colnames(mat) <- colnames(mat.df)
rownames(mat) <- rownames(mat.df)
rm(list=c("mat.df","mat.read"))

data.mat <- t(mat)[,1:2]
diff.mat <- name.check(tre.input,data.mat)
tre.trim <- drop.tip(tre.input,diff.mat$tree_not_data)
rm(data.mat)
rm(diff.mat)

taxaintree <- tre.trim$tip.label
dist.matrix <- dist.matrix.read[taxaintree,taxaintree]
rm(dist.matrix.read)
max.dist <- max(rowSums(dist.matrix))
sim.matrix <- max(dist.matrix) - dist.matrix


######

testnum <- dim(mat)[1]
testnum <- 3
mat.mini <- t(mat[1:testnum,taxaintree])

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
	for (t.genome in names(cluster)) {
		gen.dists <- sort(dist.matrix[t.genome,])
		gen.sims <- (max(gen.dists) - gen.dists)/max(gen.dists)
		n.leaves <- cluster[names(gen.dists)]
		m.surp[t.genome,t.cluster] <- sum((n.leaves[t.genome] - n.leaves) * as.double(gen.sims))
	}
#	plotTree.wBars(tre.trim,x=mat.mini[,t.cluster]/max(mat.mini[,t.cluster]+1),scale=0.5)
#	plotTree.wBars(tre.trim,x=m.surp[,t.cluster]/max(m.surp[,t.cluster]+1),scale=0.5)

#	SIMMAP for stochastic mapping
#	tre.trim.simmap <- make.simmap(tre.trim, cluster, model="SYM", nsim=1)
#	max.states <- as.integer(max(rownames(tre.trim.simmap$Q)))
#	sim.colors <- setNames(colorRampPalette(c('blue','red'))(max.states+1),0:max.states)
#	plotTree.wBars(tre.trim.simmap,x=mat.mini[,t.cluster]/max(mat.mini[,t.cluster]+1),scale=0.5,method="plotSimmap",colors=sim.colors,lwd=3,node.numbers=T)
#	plotTree.wBars(tre.trim.simmap,x=m.surp[,t.cluster]/max(m.surp[,t.cluster]+1),scale=0.5,method="plotSimmap",colors=sim.colors,lwd=3,node.numbers=T)	

	}
write.csv(m.surp,"surprise.csv")
stop()



################
# read in the clusters
################
### to map taxids to refseq ids
### head -n 1 ./../HO_110001/HO_110001.1e-3_.7_2_1_link.cluster_list.r.csv | tr '\t' '\n' | ~/repo/id2id.pl refseqmulti taxalist-Cyanobacteria >taxids

library(phytools)
library(geiger)
library(genoPlotR)

mm.surp <- read.csv("surprise.csv",row.names=1)
#mm.surp <- mm.surp[,which(colSums(mm.surp)>1)]

rescale <- function(x) as.integer((100-1)/(max(x) - min(x)) * (x - min(x)) + 1)
o.surp <- apply(mm.surp[,which(rowSums(mat)>0)],2,rescale)
o.surp[is.na(o.surp)] <- 0
colnames(o.surp) <- 1:dim(o.surp)[2]
rownames(o.surp) <- rownames(mm.surp)

rescale.all <- function(x) as.integer((100-1)/(max(p.surp) - min(p.surp)) * (x - min(p.surp)) + 1)
p.surp <- scale(mm.surp)
p.surp[is.nan(p.surp)] <- 0
pp.surp <- apply(p.surp[,which(rowSums(mat)>0)],2,rescale.all)
colnames(pp.surp) <- 1:dim(pp.surp)[2]
rownames(pp.surp) <- rownames(mm.surp)


cl_input <- scan("test_cluster_list",what="",sep="\n")
cl_list <- strsplit(cl_input, "[[:space:]]+")
num_list <- lapply(seq_along(cl_list), function(i){ rep.int(i, length( cl_list[[i]] )) })
clusters <- data.frame("cluster"=unlist(num_list),"proteinid"=unlist(cl_list))

nc_list <- list.files(pattern="NC_0\\d+\\.gbk")
nc_list2 <- nc_list2 <- sapply(nc_list,FUN=function(x) {y <- strsplit(x,".",fixed=T);y[[1]][1] } )

taxids <- read.csv("taxids",sep="\t",header=F,row.names=1)
g_list <- vector("list",length(nc_list))
names(g_list) <- taxids[nc_list2,1]
for (gbk_num in 1:length(nc_list)) {
	print(nc_list2[gbk_num])
	taxid <- taxids[gbk_num,1]
	g1 <- read_dna_seg_from_genbank(nc_list[gbk_num],gene_type="arrows")
#	g_list[[length(g_list)+1]] <- merge(g1,clusters)
	g_list[[taxid]] <- merge(clusters,g1,by="proteinid")
	g_list[[taxid]] <- cbind(g_list[[taxid]][,c(-1,-2)],list("proteinid"=g_list[[taxid]][,1],"cluster"=g_list[[taxid]][,2]) )
}

library(RColorBrewer)
library(colorRamps)
s.color <- colorRampPalette(c('blue','white','red'))(101)

h_list <- vector("list", length(g_list))

names(h_list) <- names(g_list)
for (i in 1:length(h_list)) {
	taxid <- names(h_list)[i]
	print(taxid)
	#exclude plasmids for now
	if (taxid %in% rownames(pp.surp)) {
	q.surp <- data.frame("cluster"=1:dim(pp.surp)[2],"surprise"=as.double(t(pp.surp[taxid,])))
#	h_list[[length(h_list)+1]] <- merge(data.frame(g_list[[taxid]]),q.surp)
	h_list[[taxid]] <- merge(q.surp,g_list[[taxid]],by="cluster")
	h_list[[taxid]] <- cbind(h_list[[taxid]][,c(-1,-2)],list("cluster"=h_list[[taxid]][,1],"surprise"=h_list[[taxid]][,2]) )	
	h_list[[taxid]] <- h_list[[taxid]][order(h_list[[taxid]]['start']),]
	h_list[[taxid]][,"col"] <- s.color[1+h_list[[taxid]][,"surprise"]]	
	} else {}
}

h_list <- h_list[lapply(h_list,length)>0]

#### annotations and comparisons
#### annotations and comparisons
#### annotations and comparisons


tre.cyano <- read.tree("CY_110002.phyml.tre")
data.cyano <- matrix(seq_along(names(h_list)))
rownames(data.cyano) <- names(h_list)
diff.cyano <- name.check(tre.cyano,h_list)
tre.trim <- drop.tip(tre.cyano,diff.cyano$tree_not_data)

write.tree(tre.trim,file="CY_110002.phylog.tre")

tree_file <- readLines('CY_110002.phylog.tre')
tree_phylog <- newick2phylog(tree_file)


a_list <- vector("list", length(h_list)) #dna_seg
a_annot <- vector("list", length(h_list))
c_list <- vector("list", length(h_list))
names(a_list) <- names(tree_phylog$leaves)
names(a_annot) <- names(tree_phylog$leaves)
names(c_list) <- names(tree_phylog$leaves)

for (taxid in names(tree_phylog$leaves)) {
	print(taxid)
	a_list[[taxid]] <- dna_seg(h_list[[taxid]])
	a_annot[[taxid]] <- annotation(x1=middle(a_list[[taxid]]), text=as.character(a_list[[taxid]][,"cluster"]), rot=0)
}

#mat.read[grep("nitrate",mat.read[,"description"]),]
#157                      [9|28] nitrate ABC transporter inner membrane subunit
#546                                                  [15|11] nitrate reductase
#1068                            [8|13] nitrate ABC transporter ATPases C and D
#1269                                                   [1|2] nitrate reductase
#1740  [2|13] nitrate/sulfonate/bicarbonate ABC transporter periplasmic protein
#4847                                           [1|4] nitrate transport protein
#16504                        [1|1] nitrate/nitrite porter (NNP) family protein
#170     [52|6] ammonium transporter
#4184    [1|3] ammonia monooxygenase
#6154  [2|2] histidine ammonia-lyase
#12100    [1|1] ammonium transporter

xlims <- vector("list", length(h_list))
names(xlims) <- names(tree_phylog$leaves)
#req_cluster <- c(205,1184,2402,2403,3087,16623)
req_cluster <- c(546,1269)

dist_offset <- 5000
for (taxid in names(tree_phylog$leaves)) {
	print(taxid)
	req_row <- which(a_list[[taxid]][,"cluster"] %in% req_cluster)
	if(length(req_row)>0) {
		req_row_after <- req_row + a_list[[taxid]][req_row,"strand"]
		req_row <- req_row[order(a_list[[taxid]][req_row_after,"cluster"])]
		for (i in req_row) {
			if (a_list[[taxid]][i,"strand"] > 0) {
				startend <- list(start=a_list[[taxid]][i,"start"]-dist_offset,end=a_list[[taxid]][i,"end"]+dist_offset)
			} else {
				startend <- list(start=a_list[[taxid]][i,"end"]+dist_offset,end=a_list[[taxid]][i,"start"]-dist_offset)
			}
			xlims[[taxid]] <- cbind(xlims[[taxid]],startend)			
			xlims[[taxid]][xlims[[taxid]]<0] <- 0
		}
	} else {
		xlims[[taxid]] <- c(0,0)
	}
}

xlims <- lapply(xlims,unlist)
plot.width <- 5+max(unlist(lapply(xlims,length)))*4

taxalist <- read.csv("taxids-Cyanobacteria",sep="\t",row.names=1,header=F)

tre.rename <- tre.trim
tre.rename$tip.label <- taxalist[tre.rename$tip.label,2]
write.tree(tre.rename,file="CY_110002.phylog_rename.tre")

tree_file_rename <- readLines('CY_110002.phylog_rename.tre')
tree_phylog_rename <- newick2phylog(tree_file_rename)
taxalist_rename <- data.frame(names(tree_phylog$leaves),names(tree_phylog_rename$leaves),row.names=1)

p_list <- a_list
p_annot <- a_annot
names(p_list) <- taxalist_rename[names(a_list),1]
names(p_annot) <- taxalist_rename[names(a_list),1]

p_list <- p_list[names(tree_phylog_rename$leaves)]
p_annot <- p_annot[names(tree_phylog_rename$leaves)]
p_title <- paste(as.character(mat[req_cluster,"description"]),collapse="; ")

pdf(file="fig-tree.pdf",width=plot.width,height=18)
plot_gene_map(dna_segs=p_list,annotations=p_annot, tree=tree_phylog_rename, xlims=xlims, tree_branch_labels_cex=0,tree_width=6,main=p_title)
dev.off()








#### plotting heatmaps
stop()

#http://students.washington.edu/bowmanjs/wordpress/?p=562

plot(tre.trim)
root.trim <- root(tre.trim,resolve.root=T,interactive=T)
chr.trim <- chronos(root.trim)
dend.trim <- as.dendrogram(as.hclust.phylo(chr.trim))
clade_order <- order.dendrogram(dend.trim)
clade_name <- labels(dend.trim)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name,row.names(m.surp))
new.surp <- m.surp[new_order,]

library(RColorBrewer)
library(colorRamps)
library(gplots)

color <- colorRampPalette(c('white','blue','red'))(100)

heatmap.2(new.surp[,1:100],Rowv=dend.trim,Colv=F,dendrogram='row',key=T,density.info=c('none'),scale="col",na.color='black',trace='none',margins=c(7,10),sepcolor="white",sepwidth=c(0.01,0.01),col=color,cexCol = 0.7,cexRow = 0.7)






#		lines(1:n.leaf,cumsum(n.leaves)/sum(n.leaves));
#		lines(gen.dists/max(gen.dists),cumsum(n.leaves)/sum(n.leaves));
#		points(gen.dists,cumsum(n.counts)/1:length(n.counts),lty=2)
#		lines(1:length(n.counts),cumsum(n.counts)/1:length(n.counts),lty=1)
#		lines(gen.dists,cumsum(n.counts)/1:length(n.counts),lty=1)
#		lines(gen.dists[1:length(n.counts)-1],cumsum(diff(cumsum(n.counts))),lty=1)
#		lines(gen.dists,cumsum(n.counts))
#		n.surp <- spl$y[1] - spl$y[max(which(spl$x<n.cutoff))]
#		spl <- approx(cumsum(n.counts)~as.double(gen.dists),n=100)
#		n.surp <- spl$y[max(which(spl$x<0.05))] - spl$y[1]
#		m.surp[t.genome,t.cluster] <- n.surp
#		spl <- approx(cumsum(n.leaves)/seq_along(n.leaves)~as.double(gen.dists),n=100)
#		n.surp <- sum(spl$y[1] - spl$y)

#		n.counts <- 1*(n.leaves>0)
#		lines(gen.dists,n.leaves[1] - cumsum(n.leaves)/seq_along(n.leaves))
#		n.surp <- n.leaves[1] - cumsum(n.leaves)/seq_along(n.leaves)
#		w.surp <- sum(n.surp * gen.sims)
#		m.surp[t.genome,t.cluster] <- w.surp
		
#		m.surp <- sum((n.leaves[1] - n.leaves) * gen.sims)

#	plotTree.wBars(tre.trim,x=mat.mini[,t.cluster]/max(mat.mini[,t.cluster]+1),scale=0.5)
#	plotTree.wBars(tre.trim,x=m.surp[,t.cluster]/max(m.surp[,t.cluster]+1),scale=0.5)




