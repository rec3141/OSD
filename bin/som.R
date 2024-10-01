require(phytools)
require(kohonen)

coolBlueHotRed <- function(n, alpha = 1) {
rainbow(n, end=4/6, alpha=alpha)[n:1]
}


pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')

taxa <- read.csv("key.txt",sep="\t",row.names=1)

# DATASET 1
# cols are genomes, rows are gene families
mat.arg_lys <- read.csv(file="med.arg_lys.csv",row.names=1)

# DATASET 2
#cols are amino acid properties, rows are proteins
cluster.dat <- read.csv("cluster_all_sort.csv",sep="\t",header=F)
colnames(cluster.dat) <- c("cluster","id","refseq","length","arg_lys","acidic","aliphaticity","aliphatic_index","proline","gravy","blank")

cluster.dat$arg_lys <- log(cluster.dat$arg_lys+1)
cluster.dat$acidic <- 1-cluster.dat$acidic

mat <- cluster.dat[,5:10]
mat[is.na(mat)] <- 0 #better to consider only complete cases

training.indices <- sample(nrow(mat), 20000)
training <- rep(FALSE, nrow(mat))
training[training.indices] <- TRUE

data_train <- as.matrix(scale(mat[training, ]))

som_grid <- somgrid(xdim=32, ydim=1, topo="hexagonal")

som_model <- som(data_train, 
		grid=som_grid, 
		rlen=500, 
		alpha=c(0.05,0.001), 
		keep.data = TRUE,
		n.hood="circular",
		toroidal = T)

#PREDICTION

#pdf(file="genomes-som-60x1.pdf",width=120,height=400)
#par(mfrow=c(length(ids),1))

library(ade4)
library(geiger)
library(phytools)
#tre.input <- read.tree('HI_110003.phylog_rename.tre')
tre.input <- read.tree("HI_110003.gb.concat.phylip_phyml_tree.tre")


pka <- matrix(data=NA,ncol=32,nrow=length(tre.input$tip.label))
dimnames(pka) <- list(tre.input$tip.label, 1:32)

som_save <- list()
ids <- unique(cluster.dat$id)

for (id in tre.input$tip.label) {
	som_map <- map.kohonen(som_model, newdata=as.matrix(scale(mat[cluster.dat$id==id,])))
	som_save <- rbind(som_save,som_map)
	try(
		pka[id,] <- plot(som_model, type="count",classif=som_map, main = taxa[id,1])
		)
}
pka[is.na(pka)] <- 0

#dev.off()

pdf(file="genomes-som-codes.pdf",width=12, height=6)
plot(som_model, type="codes")
dev.off()


tree_plot <- plot.phylo(tre.input,mar=c(1.1,1.1,1.1,1.1))
nodelabels()
tre.reroot <- root(tre.input,resolve.root=T,interactive=T)
write.tree(tre.reroot,file="som.phylo_reroot.tre")
tre.reroot <- read.tree("som.phylo_reroot.tre")
tre.rename <- tre.reroot
tre.rename$tip.label <- as.character(taxa[tre.rename$tip.label,1])


pdf(file="genomes-som-32x1.pdf",width=36,height=112)

layout(matrix(c(1,2),1,2),c(0.5,0.5))
#had to determine 197.5 empirically -- not sure why?
tree_plot <- plot.phylo(tre.rename, cex=3, y.lim=c(1,197.5)) #mar=c(1.1,1.1,1.1,0),

par(mar=c(1.1,0,1.1,1.1))
i<-1
id <- tre.reroot$tip.label[i]
plot(0,0, xlim=c(0,32),ylim=c(0,length(tre.rename$tip.label)), pch=21,bg = heat.colors( max(pka[id,]) )[pka[id,]], cex=7, type="n", bty="n", axes=F)
for (i in 1:length(tre.rename$tip.label)) {
	id <- tre.reroot$tip.label[i]
	points(1:32,rep(i+0.5,32), pch=21, bg = heat.colors( max(pka[id,]) )[pka[id,]], cex=7)
}

dev.off()

pdf(file="amino-acid-som.pdf",width=12,height=12)
par(mfrow=c(3,3))

#plot(som_model, type="changes")
plot(som_model, type="count")
#plot(som_model, type="dist.neighbours")
plot(som_model, type="codes")

#par(mfrow=c(1,2))
numclusters <- 4
#plot(hclust((dist(som_model$codes))))
som_cluster <- cutree(hclust(dist(som_model$codes)), numclusters)
# plot these results:
plot(som_model, type="mapping", bgcol = heat.colors(numclusters)[som_cluster], main = "Clusters") 
#add.cluster.boundaries(som_model, som_cluster)

hist(som_cluster,col=heat.colors(numclusters),breaks=1:(numclusters+1),right=F)

for(i in 1:6) {
	plot(som_model, type = "property", property = som_model$codes[,i], main=names(mat)[i], palette.name=coolBlueHotRed)
}
dev.off()



# how many clusters?
mydata <- som_model$codes 
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var)) 
for (i in 2:50) {
  wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
}
plot(wss)








