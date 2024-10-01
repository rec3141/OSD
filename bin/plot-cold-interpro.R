library(ggplot2)
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
pval.files <- list.files(dir1,pattern="out.pvals.*")
pval.all <- NULL
pval.weight <- data.frame()
coldpar <-  c("acidic","aliphaticity","aliphatic_index","arg_lys","gravy","proline","mweight")
weights.log <- list()

## heart plot

for (g in c("temperature","salinity","latitude")) {
	
	weights.coldpar <- matrix(data=NA,nrow=length(pval.files),ncol=length(coldpar))
	rownames(weights.coldpar) <- pval.files
	colnames(weights.coldpar) <- coldpar

for (my.file in pval.files) {
	pval.in <- read.csv(paste(dir1,my.file,sep="/"),head=T)
	pval.all <- rbind(pval.all,pval.in)
	
	for (f in levels(pval.in$my.coldpar)) {
	
	#get subset for each cold parameter
	 d <- pval.in[which(pval.in$my.coldpar==f),]

	gg <- data.frame(c(d[paste(g,"p",sep=".")],d[paste(g,"slope",sep=".")]))
	colnames(gg) <- c("coldpar.p","coldpar.slope")

	#replace p-values of 0 with smallest nonzero p-value
	 gg$coldpar.p[gg$coldpar.p == 0] <- 10^min(log10(gg$coldpar.p[gg$coldpar.p>0]))
	#calculate weighted mean based on least significant hits
	 d.mean <- mean(gg$coldpar.slope[gg$coldpar.p > (.05/length(gg$coldpar.p)/100) ]) #bonferroni
	#calculate standard deviations from mean
	 d.sd <- sd(abs(gg$coldpar.slope))
	#calculate maximum value of slope
#	 d.abs <- max(abs(gg$coldpar.slope),na.rm=T)
	 d.abs <- d.sd
	#calculate normalized slope values
	 d.norm <- (gg$coldpar.slope)/d.abs

	#calculate normalized p-value weighted mean 
	 d.weight <- sum(-log10(gg$coldpar.p) * d.norm,na.rm=T)/dim(d)[1]
	#color weighted mean line blue if above zero, red if below
	 d.color = "blue"; if (d.weight < 0) d.color="red"
	#set ylim to 1 because it has been normalized
	 d.ylim <- 1
#	 d.ylim <- max(d.abs,abs(d.intercept)/10)*1.1
#	 if (d.ylim > 0.05) d.ylim <- 0.05
	 weights.coldpar[my.file, f] <- d.weight #?d.mean
	
	 pval.weight <- rbind(pval.weight,cbind("file"=my.file,"my.coldpar"=f, "weight"=d.weight))
	 ggplot(as.data.frame(gg) , aes(x = -log10(coldpar.p), y = d.norm, colour = f, size=3)) + 
	 	geom_point(alpha = 0.3) + scale_colour_discrete(drop=TRUE,limits = levels(factor(d$my.coldpar))) + 
	 	geom_abline(intercept = d.weight, slope = 0, colour=d.color, size=1.5) + 
	 	geom_abline(intercept = 0, slope = 0, colour="black", size=1.5) + 
	 	ylim(-4*d.ylim,4*d.ylim)	+ scale_x_log10(limits=c(1e-3,1e2)) # + scale_x_log10()
	 ggsave(paste("./fig",g,d$my.taxa[1],f,"pdf",sep="."))

	} 
		
	}
	
	weights.log[[g]] <- weights.coldpar

}

for (g in c("temperature","salinity","latitude")) {
	pdf(file=paste("heatmap",g,"pdf",sep="."),width=10,height=20)
	heatmap((as.matrix(weights.log[[g]])),scale="none",mar=c(20,20),revC=T,col=brewer.pal(11,"RdBu"),main=g)
	dev.off()
}

##write.csv(pval.weight,file="pvals.weights.latitude.taxa.csv")
