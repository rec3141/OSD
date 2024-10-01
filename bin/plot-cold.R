#library(ggplot2)

pval.files <- list.files(".",pattern="pvals.new*")
pval.all <- NULL
pval.weight <- data.frame()

for (my.file in pval.files) {
	pval.in <- read.csv(my.file,head=T)
	pval.all <- rbind(pval.all,pval.in)

#	if(dim(pval.in)[1] < 250) { next }	

	for (f in levels(pval.in$my.coldpar)) {
	 d <- pval.in[which(pval.in$my.coldpar==f),]
	d.weight <- sum(-log(d$temperature.p) * d$temperature.slope,na.rm=T)/dim(pval.in)[1]/5
	pval.weight <- rbind(pval.weight,cbind("file"=my.file,"my.coldpar"=f, "weight"=d.weight))
#	 ggplot(d , aes(x = -log(temperature.p), y = temperature.slope, colour = my.coldpar, size=3)) + geom_point(alpha = 0.3) + scale_colour_discrete(drop=TRUE,limits = levels(factor(d$my.coldpar))) # + scale_x_log10()
#	 ggsave(paste("./figs/fig",d$my.taxa[1],f,"pdf",sep="."))
	} 
}
write.csv(pval.weight,file="pvals.weights.taxa.csv")


pval.weight.int <- NULL

for (g in levels(pval.all$my.coldpar)) {
	for (h in levels(pval.all$my.protfam)) {

	d <- pval.all[which(pval.all$my.coldpar==g & pval.all$my.protfam==h),]
	d <- d[grepl("g__*",d$my.taxa),] #just look at genera level

	if(dim(d)[1] < 30) { cat("."); next }
	 d.abs <- max(abs(d$temperature.slope),na.rm=T)
	 d.intercept <- sum(-log(d$temperature.p) * d$temperature.slope,na.rm=T)/dim(d)[1]
	 d.color = "blue"; if (d.intercept < 0) d.color="red"
	 d.ylim <- max(d.abs,abs(d.intercept))*1.1
	pval.weight.int <- rbind(pval.weight.int,cbind("file"=h,"my.coldpar"=g, "weight"=d.intercept))

#	 ggplot(d , aes(x = -log(temperature.p), y = temperature.slope, colour = my.coldpar, size=3)) + geom_point(alpha = 0.3)  + scale_colour_discrete(drop=TRUE,limits = levels(factor(d$my.coldpar)))  	+ geom_abline(intercept = d.intercept, slope = 0, colour=d.color, size=1.5) + geom_abline(intercept = 0, slope = 0, colour="black", size=1.5) + ylim(-1*d.ylim,d.ylim)	+ scale_x_log10(limits=c(1e-3,1e2)) 
#	 ggsave(paste("./figs/figs",g,h,"pdf",sep="."))
		
	}
}
write.csv(pval.weight.int,file="pvals.weights.interpro.csv")

#	 d[order(d$temperature.p)[1:20],] 
#	 plot(pval.in$temperature.p ~ pval.in$temperature.slope,log="y",col=rgb(1, 0, 0, 0.3), pch=16,cex=1.5)
#	 plot.id <- identify(pvals$temperature.slope,pvals$temperature.p)
#	 pvals[plot.id,]
#	 ggplot(d , aes(x = reorder(my.protfam,as.numeric(temperature.p)), y = temperature.p, colour = my.tax)) + geom_point() + scale_colour_discrete(drop=TRUE,limits = levels(pval.in$my.tax)) + scale_y_log10()


