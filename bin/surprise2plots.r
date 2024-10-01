################
################
rm(list=ls())

library(phytools)
library(geiger)
library(genoPlotR)
library(RColorBrewer)
library(colorRamps)

#1. 
surprisefile <- "acidic.surprise.csv"
mm.surp <- read.csv(surprisefile,row.names=1)
#mm.surp <- mm.surp[,which(colSums(mm.surp)>1)]
mm.surp <- mm.surp[1:2000,]

#rescale <- function(x) as.integer((100-1)/(max(x) - min(x)) * (x - min(x)) + 1)
#o.surp <- apply(mm.surp,2,rescale)
#o.surp[is.na(o.surp)] <- 0
#colnames(o.surp) <- 1:dim(o.surp)[2]
#rownames(o.surp) <- rownames(mm.surp)

#rescale to the range 0-100
#rescale.all <- function(x) as.integer((100-1)/(max(p.surp) - min(p.surp)) * (x - min(p.surp)) + 1) 
#rescale to the range 0-100 with mean = 50, stdev=5

p.surp <- scale(mm.surp) #scales to mean=0, sd=1
p.surp[is.nan(p.surp)] <- 0 #assume missing values=mean
#pp.surp <- apply(p.surp,2,rescale.all)
#rescale to mean=50, sd=5, force range to 0--100
pp.surp <- p.surp*5+50
pp.surp[pp.surp<0] <- 0
pp.surp[pp.surp>100] <- 100
rownames(pp.surp) <- 1:dim(pp.surp)[1]
colnames(pp.surp) <- colnames(mm.surp)

#2. taxalist
# manually corrected CL_016196 NC_0039112 NC_006569 to CL_016196 NC_003911 NC_006569
taxalistfile <- "taxalist-cellular-rec"
taxalist.read <- scan(taxalistfile,sep="\n","char")
taxalist.list <- strsplit(taxalist.read, "[[:space:]]+")
taxalist.univ <- lapply(1:length(taxalist.list),
	function(i){
	paste( rep(taxalist.list[[i]][[1]], length(taxalist.list[[i]]))
		 , 0:( length(taxalist.list[[i]]) - 1)
		 , sep="."
		);
	})
taxalist.df <- data.frame("univ"=unlist(taxalist.univ),"refseq"=unlist(taxalist.list))
taxalist.df <- taxalist.df[-grep("CL_", taxalist.df$refseq,fixed=T),]
taxalist.df$refseq <- substr(taxalist.df$refseq,1,10)

#3. *.tre -- newick tree file with same genomes as cluster_list
trefile <- "HI_110003.gb.concat.phylip_phyml_tree.tre"
if(!file.exists(trefile)) stop(paste(trefile,"does not exist",sep=" "))
tre.input <- read.tree(trefile)

taxalist.tree <- lapply(tre.input$tip.label, function(x) grep(x,taxalist.read,value=T))
taxalist.tree <- as.data.frame(cbind(tre.input$tip.label,substr(taxalist.tree,1,9)))
colnames(taxalist.tree) <- c("refseq","univ")

tre.input$tip.label <- as.character(taxalist.tree$univ)

clusterlistfile <- "cluster_list_2000"

	cl_input <- scan(clusterlistfile,what="",sep="\n")
	cl_list <- strsplit(cl_input, "[[:space:]]+")
	num_list <- lapply(seq_along(cl_list), function(i){ rep.int(i, length( cl_list[[i]] )) })
	clusters <- data.frame("cluster"=unlist(num_list),"proteinid"=unlist(cl_list))

	nc_list <- list.files(path="./gbk/",pattern="N._.*\\.gbk")
	nc_list2 <- substr(as.character(sapply(nc_list,FUN=function(x) {y <- strsplit(x,".",fixed=T);y[[1]][1] } )),1,10)
	nc_list2 <- cbind(nc_list,nc_list2)
	colnames(nc_list2) <- c("gbkfile","refseq")

taxids <- as.data.frame(cbind(taxalist.df$refseq,as.character(taxalist.df$univ)))
colnames(taxids) <- c("refseq","univ")
taxids.tre <- taxids[unlist(lapply(taxalist.tree$univ, function(x) grep(x,taxids$univ))),]

tax_nc <- merge(taxids,nc_list2)
tax_nc <- cbind(tax_nc,substr(tax_nc$univ,1,9))
colnames(tax_nc) <- c("refseq","univ","gbkfile","univs")

	g_list <- vector("list",length(tax_nc[,1]))
	names(g_list) <- tax_nc$refseq

	for (gbk_num in 1:length(tax_nc[,1])) {
		print(tax_nc[gbk_num,1])
		taxid <- as.character(tax_nc[gbk_num,1])
		
		g1 <- try(
		{read_dna_seg_from_genbank(paste("./gbk/",tax_nc[gbk_num,3],sep=""),gene_type="arrows")
		g_list[[taxid]] <- merge(clusters,g1,by="proteinid")
	#	g_list[[length(g_list)+1]] <- merge(g1,clusters)
		g_list[[taxid]] <- cbind(g_list[[taxid]][,c(-1,-2)],list("proteinid"=g_list[[taxid]][,1],"cluster"=g_list[[taxid]][,2]) )
			})
	}


#mean(as.numeric(unlist(lapply(h_list,function(x) x[,"surprise"]))))
#sd(as.numeric(unlist(lapply(h_list,function(x) x[,"surprise"]))))

s.color <- colorRampPalette(c('blue','white','red'))(101)

	h_list <- vector("list", length(g_list))
	names(h_list) <- names(g_list)

	for (i in 1:length(h_list)) {
		taxid <- names(h_list)[i]
		print(taxid)
		#exclude plasmids for now
		if (tax_nc$univs[i] %in% colnames(pp.surp)) {
		q.surp <- data.frame("cluster"=1:dim(pp.surp)[1],"surprise"=as.double(t(pp.surp[,as.character(tax_nc$univs[i])])))
	#	h_list[[length(h_list)+1]] <- merge(data.frame(g_list[[taxid]]),q.surp)
		try({
		h_list[[taxid]] <- merge(q.surp,g_list[[taxid]],by="cluster")
		h_list[[taxid]] <- cbind(h_list[[taxid]][,c(-1,-2)],list("cluster"=h_list[[taxid]][,1],"surprise"=h_list[[taxid]][,2]) )	
		h_list[[taxid]] <- h_list[[taxid]][order(h_list[[taxid]]['start']),]
		h_list[[taxid]][,"col"] <- s.color[1+h_list[[taxid]][,"surprise"]]
		})
		} else {}
	}

	h_list <- h_list[lapply(h_list,length)>0]


#### annotations and comparisons
#### annotations and comparisons
#### annotations and comparisons

h_list_names <- tax_nc$univs[as.numeric(lapply(names(h_list), function(x) which(x == tax_nc$refseq)))]
names(h_list) <- h_list_names

#data.cyano <- matrix(seq_along(names(h_list)))
#rownames(data.cyano) <- names(h_list)
diff.input <- name.check(tre.input,h_list)
tre.trim <- drop.tip(tre.input,diff.input$tree_not_data)

write.tree(tre.trim,file="tmp.tre.trim")

tree_file <- readLines('tmp.tre.trim')
tree_phylog <- newick2phylog(tree_file)


a_list <- vector("list", length(h_list)) #dna_seg
a_annot <- vector("list", length(h_list))
c_list <- vector("list", length(h_list))
names(a_list) <- names(tree_phylog$leaves)
names(a_annot) <- names(tree_phylog$leaves)
names(c_list) <- names(tree_phylog$leaves)

for (taxid in names(tree_phylog$leaves)) {
	print(taxid)
	a_list[[taxid]] <- try(dna_seg(h_list[[taxid]]))
	a_annot[[taxid]] <- try(annotation(x1=middle(a_list[[taxid]]), text=as.character(a_list[[taxid]][,"cluster"]), rot=0))
}

	a_list <- a_list[lapply(a_list,length)>0]
	a_list <- a_list[lapply(a_list,class)!="try-error"]
	a_annot <- a_annot[lapply(a_annot,length)>0]
	a_annot <- a_annot[lapply(a_annot,class)!="try-error"]


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

#taxalist <- read.csv("taxids-Cyanobacteria",sep="\t",row.names=1,header=F)

diff.input <- name.check(tre.trim,a_list)
tre.rename <- drop.tip(tre.trim,diff.input$tree_not_data)

tre.names <- read.csv("key.txt",sep="\t")
tre.names[,1] <- substr(tre.names[,1],1,10)
colnames(tre.names) <- c("refseq","taxaname")
tax_nc <- merge(tax_nc,tre.names)
tre.rename$tip.label <- as.character(tax_nc[as.numeric(lapply(tre.rename$tip.label,function(x) which(x==as.character(tax_nc$univs))[1])),"taxaname"])
tre.rename$node.label <- as.numeric(1:length(tre.rename$node.label))+length(tre.rename$tip.label)


# RUN MANUALLY
plot(tre.rename)
nodelabels()
root.rename <- root(tre.rename,resolve.root=T,interactive=T)
#root.rename <- root(tre.rename,node=123)
plot(root.rename)
write.tree(root.rename,file="HI_110003.phylog_rename.tre")

tree_file_rename <- readLines('HI_110003.phylog_rename.tre')
tree_phylog_rename <- newick2phylog(tree_file_rename)
#taxalist_rename <- data.frame(names(tree_phylog$leaves),names(tree_phylog_rename$leaves),row.names=1)

#the top 2000 gene families have 363204/534026 = 68% of the genes

for (req_cluster in 2000:1) {
print(req_cluster)
xlims <- vector("list", length(h_list))
names(xlims) <- names(tree_phylog$leaves)

dist_offset <- 5000
for (taxid in names(tree_phylog$leaves)) {
#	print(taxid)
	if(length(a_list[[taxid]]) > 1) {
	req_row <- which(a_list[[taxid]][,"cluster"] %in% req_cluster)
	if(length(req_row)>0) {
		req_row_after <- req_row + a_list[[taxid]][req_row,"strand"]
		req_row <- req_row[order(a_list[[taxid]][req_row_after,"cluster"])] #weark effort to sort by operon
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
	} else { xlims[[taxid]] <- c(0,0) }
}

xlims <- lapply(xlims,unlist)

plot.width <- 10+max(unlist(lapply(xlims,length)))*4

p_list <- a_list
p_annot <- a_annot
p_list <- p_list[lapply(p_list,length)>0]
p_annot <- p_annot[lapply(p_annot,length)>0]

xlims <- xlims[intersect(names(xlims),names(p_list))]

names(p_list) <- names(tree_phylog_rename$leaves)
names(p_annot) <- names(tree_phylog_rename$leaves)
names(xlims) <- names(tree_phylog_rename$leaves)

#p_title <- paste(as.character(mat[req_cluster,"description"]),collapse="; ")

pdf(file=paste("aminos-figs-trees/acidic.",req_cluster,".pdf",sep=""),width=plot.width,height=64)
try(plot_gene_map(dna_segs=p_list,annotations=p_annot, tree=tree_phylog_rename, xlims=xlims, tree_branch_labels_cex=0,tree_width=12,main=""))
dev.off()


# AMINOS.* MIGHT BE FUCKED UP e.g. no Colwellia?!
}




## STOPPED HERE






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


