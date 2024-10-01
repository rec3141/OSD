library(ggplot2)
library(gplots)
library(RColorBrewer)
options(width=200)

#from Psychrobacter genome paper
#they were less hydrophobic, had fewer proline residues, 
#were less aliphatic, had fewer acidic residues, 
#or had low Arg and increased Lys contents

# acidic residues (14 categories) and proline (10 categories) 
# are the two most frequent cold adaptation indicators across all COG categories

#The categories with higher numbers of adaptations include the following COG categories: 
#no COG designation; replication, recombination, and repair; amino acid transport and 
#metabolism; lipid transport and metabolism; transcription; translation, ribosomal structure, 
#and biogenesis; and signal transduction mechanism

dir1 <- "./../step3-merge"
pval.files <- list.files(dir1,pattern="out.pvals.*")[1:5]
pval.all <- NULL
pval.weight <- data.frame()
#old coldpar.list <-  c("acidic","aliphaticity","aliphatic_index","arg_lys","gravy","proline","mweight")
coldpar.list <- c("acidic","aliphaticity","aliphatic_index","arg_lys","gravy","polar","nonpolar","sulfur","hbond","basic","aromatic","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","nitrogen","mweight")
weights.log <- list()
slopes.log <- list()


## heart plot

for (g in c("temperature","salinity","latitude")) {
	
	weights.coldpar <- matrix(data=NA,nrow=length(pval.files),ncol=length(coldpar.list))
	slopes.coldpar <- matrix(data=NA,nrow=length(pval.files),ncol=length(coldpar.list))
	rownames(weights.coldpar) <- as.character(sapply(pval.files, function(x) paste(strsplit(x,'[.]')[[1]][3])))
	colnames(weights.coldpar) <- coldpar.list
	rownames(slopes.coldpar) <- as.character(sapply(pval.files, function(x) paste(strsplit(x,'[.]')[[1]][3])))
	colnames(slopes.coldpar) <- coldpar.list


for (my.file in pval.files) {
	
	my.taxa <- strsplit(my.file,'[.]')[[1]][3]
	if ( length(grep("Psychrobacter",my.taxa))>0) { next }
	print(c(g, my.taxa))
	if (!file.exists(paste(dir1,my.file,sep="/"))) { next }
	if (file.info(paste(dir1,my.file,sep="/"))["size"] < 100) { next }
	pval.in <- read.csv(paste(dir1,my.file,sep="/"),head=T)
	pval.in <- pval.in[,c("my.coldpar", "my.taxa","my.protfam", "latitude.intercept", "latitude.slope", "latitude.t", "latitude.p", "temperature.intercept", "temperature.slope", "temperature.t", "temperature.p", "salinity.intercept", "salinity.slope", "salinity.t", "salinity.p")]
 
	if(nlevels(pval.in$my.protfam) < 30) { print("skipping... too few"); next }

	if(is.null(pval.all)) {
		pval.all <- rbind(pval.all,pval.in)
	} else {
		if(ncol(pval.in) == ncol(pval.all)) {
			pval.all <- rbind(pval.all,pval.in)
		} else { next }
	}

	for (f in levels(pval.in$my.coldpar)) {
	

	#get subset for each cold parameter
	d <- pval.in[which(pval.in$my.coldpar==f),]

#	if (file.exists(paste("./figs/",g,"/",f,"/","fig.",d$my.taxa[1],".pdf",sep=""))) { next}
	if (any(is.na(d))) { print("skipping... NAs"); next }


	gg <- data.frame(c(d[paste(g,"p",sep=".")],d[paste(g,"slope",sep=".")]))
	colnames(gg) <- c("coldpar.p","coldpar.slope")

	#replace p-values of 0 with smallest nonzero p-value
	 gg$coldpar.p[gg$coldpar.p == 0] <- 10^min(log10(gg$coldpar.p[gg$coldpar.p>0]))
	#calculate weighted mean based on least significant hits
#	 d.mean <- mean(gg$coldpar.slope[gg$coldpar.p > (.05/length(gg$coldpar.p)/100) ])
	 d.mean <- 0
	#calculate standard deviations from mean
	 d.sd <- sd(gg$coldpar.slope)
	#calculate maximum value of slope
#	 d.abs <- max(abs(gg$coldpar.slope),na.rm=T)
	 d.abs <- d.sd
	#calculate normalized slope values
	 d.norm <- (gg$coldpar.slope)/d.abs

	#calculate normalized p-value weighted mean 
#	 d.weight <- sum(-log10(gg$coldpar.p) * d.norm,na.rm=T)/dim(d)[1]
	 d.weight <- sum(d.norm,na.rm=T)/dim(d)[1]
	#color weighted mean line blue if above zero, red if below
	 d.color = "blue"; if (d.weight < 0) d.color="red"
	#set ylim to 1 because it has been normalized
	 d.ylim <- 1
#	 d.ylim <- max(d.abs,abs(d.intercept)/10)*1.1
#	 if (d.ylim > 0.05) d.ylim <- 0.05
	 weights.coldpar[my.taxa, f] <- d.weight #?d.mean
	 slopes.coldpar[my.taxa, f] <- sum(gg$coldpar.slope,na.rm=T)/dim(d)[1]
	 
	 pval.weight <- rbind(pval.weight,cbind("file"=my.file,"my.coldpar"=f, "weight"=d.weight))
	 
 	if (file.exists(paste("./figs/",g,"/",f,"/","fig.",d$my.taxa[1],".pdf",sep=""))) { next} 	

if (f %in% c("acidic","aliphaticity","polar","nonpolar","sulfur","hbond","basic","aromatic","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","nitrogen","mweight")) {

	 ggplot(as.data.frame(gg) , aes(x = -log10(coldpar.p), y = d.norm)) + #, colour = f, size=2
	 	geom_abline(intercept = 0, slope = 0, colour="black", alpha=0.3, size=1.5) + 
	 	geom_abline(intercept = d.weight, slope = 0, colour=d.color, alpha=0.3, size=1.5) + 
	 	geom_point(alpha = 0.4, colour="black", size=3) + 
	 	ylim(-4*d.ylim,4*d.ylim) + 
	 	scale_x_log10(limits=c(5e-3,1e2), breaks=c(0.004364805, 0.301029996, 1.301029996, 6, 15, 60, 180), labels=c("0.99","0.5","0.05", "1e-6", "1e-15", "1e-60", "1e-180")) +
		labs(x="p-value",y="stdev-normalized slope") +  ggtitle(paste(my.taxa,g,f,sep="\n")) + 
		annotate("text", x=80, y=1, label=signif(15*300*d.sd,2)) +
		annotate("text", x=80, y=2, label=2*signif(15*300*d.sd,2)) +
		annotate("text", x=80, y=3, label=3*signif(15*300*d.sd,2))

} else {

	 ggplot(as.data.frame(gg) , aes(x = -log10(coldpar.p), y = d.norm)) + #, colour = f, size=2
	 	geom_abline(intercept = 0, slope = 0, colour="black", alpha=0.3, size=1.5) + 
	 	geom_abline(intercept = d.weight, slope = 0, colour=d.color, alpha=0.3, size=1.5) + 
	 	geom_point(alpha = 0.4, colour="black", size=3) + 
	 	ylim(-4*d.ylim,4*d.ylim) + 
	 	scale_x_log10(limits=c(5e-3,1e2), breaks=c(0.004364805, 0.301029996, 1.301029996, 6, 15, 60, 180), labels=c("0.99","0.5","0.05", "1e-6", "1e-15", "1e-60", "1e-180")) +
		labs(x="p-value",y="stdev-normalized slope") +  ggtitle(paste(my.taxa,g,f,sep="\n")) + 
		annotate("text", x=80, y=1, label=signif(d.sd,2)) +
		annotate("text", x=80, y=2, label=2*signif(d.sd,2)) +
		annotate("text", x=80, y=3, label=3*signif(d.sd,2))

}

#		annotate("text", x=(max(-log10(gg$coldpar.p))+10), y=d.norm[which(-log10(gg$coldpar.p) == max(-log10(gg$coldpar.p)))], label=signif(gg$coldpar.slope[which(-log10(gg$coldpar.p) == max(-log10(gg$coldpar.p)))],3))
 
#		scale_colour_discrete(drop=TRUE,limits = levels(factor(d$my.coldpar))) + 
 
 		if (!file.exists(paste("./figs/",g,"/",f,sep=""))) {
 			dir.create(paste("./figs/",g,"/",f,sep=""),recursiv=T)
 		}	    
 
 	 	ggsave(paste("./figs/",g,"/",f,"/","fig.",d$my.taxa[1],".pdf",sep=""))
 
	} 
		
	}
	
	
	weights.log[[g]] <- weights.coldpar
	write.csv(weights.log[[g]],file=paste("pvals.weights",g,"csv",sep="."))
	h <- as.matrix(weights.log[[g]][!is.na(weights.log[[g]])[,1],])
	heatmap(h,scale="none",mar=c(20,20),revC=T,col=brewer.pal(11,"RdBu"),main=g)
	heatmap.2(h,scale="none",mar=c(20,20),revC=T,col=brewer.pal(11,"RdBu"),main=g)

	slopes.log[[g]] <- slopes.coldpar
	write.csv(slopes.log[[g]],file=paste("pvals.slopes",g,"csv",sep="."))
	j <- as.matrix(slopes.log[[g]][!is.na(slopes.log[[g]])[,1],])
	heatmap(j,scale="none",mar=c(20,20),revC=T,col=brewer.pal(11,"RdBu"),main=g)
	
}

for (g in c("temperature","salinity","latitude")) {
	pdf(file=paste("heatmap.normalized",g,"pdf",sep="."),width=20,height=20)
	h <- as.matrix(weights.log[[g]][!is.na(weights.log[[g]])[,1],])
#	h <- h[-grep("Psychro",rownames(h)),]
	heatmap(h,scale="none",mar=c(20,60),revC=T,col=brewer.pal(11,"RdBu"),main=g)

	h.genus <- h[grep("^g__",rownames(h)),]
	heatmap(h.genus,scale="none",mar=c(20,60),revC=T,col=brewer.pal(11,"RdBu"),main=g)

	h.species <- h[grep("^s__",rownames(h)),]
	heatmap(h.species,scale="none",mar=c(20,60),revC=T,col=brewer.pal(11,"RdBu"),main=g)

	h.order <- h[grep("^o__",rownames(h)),]
	heatmap(h.order,scale="none",mar=c(20,60),revC=T,col=brewer.pal(11,"RdBu"),main=g)

	h.higher <- h[c(grep("^c__",rownames(h)),grep("^p__",rownames(h)),grep("^k__",rownames(h))),]
	heatmap(h.higher,scale="none",mar=c(20,60),revC=T,col=brewer.pal(11,"RdBu"),main=g)

	dev.off()

	pdf(file=paste("heatmap.raw",g,"pdf",sep="."),width=20,height=20)
	h <- as.matrix(slopes.log[[g]][!is.na(slopes.log[[g]])[,1],])
#	h <- h[-grep("Psychro",rownames(h)),]
	heatmap(h,scale="col",mar=c(20,60),revC=T,col=brewer.pal(11,"RdBu"),main=g)
	dev.off()

}

