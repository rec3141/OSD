#OSD9_2014-06-18_0m_NPL022.protparam.csv #this file holds the calculated parameters from cold_info_pfam.pl
#out.OSD9_2014-06-18_0m_NPL022.pfam.tax #this file holds the metadata for each gene
#osd-metadata.csv #this file holds the metadata for each site

#rm(list=ls())

args<-commandArgs(TRUE)

file.count.taxa <- args[1] #count3.taxa
file.count.interpro <- args[2] #count3.interpro
outdir <- args[3]

file.count.taxa <- "~/Dropbox/OSD/2016-02-15/count3.genera"
file.count.interpro <- "~/Dropbox/OSD/2016-02-15/count3.interpro"
outdir <- "."

#curdir <- getwd()
mydir <- paste(outdir,"/",sep="")

osd.in.saved <- "osd.less.in.saved.Rdata"
osd.metadata <- read.csv("./osd-metadata.csv",sep="\t")

#if run interactively with source(), dont reload the whole dataset

if(exists("osd.less.in")) {
		print("using loaded database")
} else {
	if(file.exists(osd.in.saved)) {
		#load previously saved table
		cat("loading saved data from",osd.in.saved,"\n")
		#import 'osd.less.in'
		load(osd.in.saved)
	}
}

#access as osd.less.in[["OSD99_2014-06-21_1m"]]

###look for protein and taxa of interest

all.uniq.taxa <- read.csv(file.count.taxa,header=F,sep="~")
all.uniq.taxa <- sapply(all.uniq.taxa,function(x) gsub('^\\s*[[:digit:]]+\\s+','',x))
all.uniq.taxa <- as.character(sapply(all.uniq.taxa,function(x) gsub('[[:punct:]]','_',x)))
all.uniq.interpro <- read.csv(file.count.interpro,header=F,sep=" ",skip=1)
many.uniq.interpro <- all.uniq.interpro[all.uniq.interpro[,1]>=100,]
many.uniq.interpro <- many.uniq.interpro[order(many.uniq.interpro[,2]),]

#coldpar.list <- c("acidic","aliphaticity", "aliphatic_index", "arg_lys","gravy","proline","mweight") # ,"gc.pct"
coldpar.list <- c("acidic","aliphaticity","aliphatic_index","arg_lys","gravy","polar","nonpolar","sulfur","hbond","basic","aromatic","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","nitrogen","mweight")
		
if (!file.exists("figs")) { dir.create("figs")}
for (coldpar in coldpar.list) {
	if (!file.exists(paste("figs",coldpar,sep="/"))) {
	dir.create(paste("figs",coldpar,sep="/"))
	}
}

#	my.taxa <- "k__Bacteria"
for (my.taxa in all.uniq.taxa) {
	hitcount <- 0
	my.outfile <- paste(mydir,"out.pvals.",my.taxa,".csv",sep="")
	if (length(grep("__$",my.taxa)) > 0) { next }
#	if (file.exists(my.outfile)) { next }

#	file.remove(my.outfile)
	fh <- file(my.outfile,open="at")
	fhc <- 0

	cat("selecting ",my.taxa,"\n")
#	osd.taxa.in <- lapply(osd.less.in, subset, tax.domain==my.taxa | tax.phylum==my.taxa | tax.class==my.taxa | tax.order==my.taxa | tax.family==my.taxa | tax.genus==my.taxa | tax.species==my.taxa, select=c("interpro",coldpar.list))

	if (length(grep('k__',my.taxa))) { osd.taxa.in <- lapply(osd.less.in, subset, tax.domain==my.taxa, select=c("interpro",coldpar.list)) }
	if (length(grep('p__',my.taxa))) { osd.taxa.in <- lapply(osd.less.in, subset, tax.phylum==my.taxa, select=c("interpro",coldpar.list)) }
	if (length(grep('c__',my.taxa))) { osd.taxa.in <- lapply(osd.less.in, subset, tax.class==my.taxa, select=c("interpro",coldpar.list)) }
	if (length(grep('o__',my.taxa))) { osd.taxa.in <- lapply(osd.less.in, subset, tax.order==my.taxa, select=c("interpro",coldpar.list)) }
	if (length(grep('f__',my.taxa))) { osd.taxa.in <- lapply(osd.less.in, subset, tax.family==my.taxa, select=c("interpro",coldpar.list)) }
	if (length(grep('g__',my.taxa))) { osd.taxa.in <- lapply(osd.less.in, subset, tax.genus==my.taxa, select=c("interpro",coldpar.list)) }
	if (length(grep('s__',my.taxa))) { osd.taxa.in <- lapply(osd.less.in, subset, tax.species==my.taxa, select=c("interpro",coldpar.list)) }
	
	taxa.protfams <- sort(as.character(unique(unlist(lapply(osd.taxa.in,function(x) as.character(unique(x$interpro)))))))
	num.protfams <- length(taxa.protfams)
	cat("found\t",num.protfams, "\tinterpro domains in\t",my.taxa,"\n")
	many.protfams <-  merge(taxa.protfams,many.uniq.interpro,by.x=1,by.y="V2")
	many.protfams <- many.protfams[order(-many.protfams[,2]),]
	num.protfams <- length(many.protfams[,1])
	cat("analyzing\t",num.protfams, "\tinterpro domains in\t",my.taxa,"\n")
	

if (1 == 1) {
# calculates linear regressions of protein domain richness/diversity against T, S, lat

options(stringsAsFactors = FALSE)
	require(vegan)
	div.mat <- as.data.frame(matrix(data = 0, nrow = length(names(osd.less.in)), ncol = 1000))
	colnames(div.mat) <- as.character(all.uniq.interpro[1:1000,2])
	rownames(div.mat) <- names(osd.less.in)

#	div.df <- data.frame(stringsAsFactors=FALSE)

		for (osd.label in names(osd.less.in)) {
		
		  osd.label.row <- which(osd.metadata$label == osd.label)
			div.ln <- dim(osd.taxa.in[[osd.label]])[1]
			
#			div.df <- rbind( div.df, cbind(
#			"osd.label"=osd.label,
#			"latitude"=osd.metadata$latitude[osd.label.row],
#			"longitude"=as.numeric(osd.metadata$longitude[osd.label.row]),
#			"temperature"=as.numeric(osd.metadata$temperature[osd.label.row]),
#			"salinity"=as.numeric(osd.metadata$salinity[osd.label.row]),
#			"prot.richness"=as.numeric(div.ln)
#			) )			
		

		div.rle <- rle(sort(as.character(osd.taxa.in[[osd.label]][,1])))
		if (length(div.rle$values)) {
		for (j in 1:length(div.rle$values)) {
			if (div.rle$values[j] == "") { next} 
			if( length(grep(div.rle$values[j],colnames(div.mat))) ) {
				div.mat[osd.label,div.rle$values[j]] <- div.rle$lengths[j]
			}
		}
		 		
		}
		}
		
#	div.df[,2:6] <- apply(div.df[,2:6],2,as.numeric)
#	q3 <- quantile(div.df$prot.richness)[3]
#	div.q3 <- div.df[div.df$prot.richness > q3,]

	rc.mid <- signif(max(rowSums(div.mat>0))/2,digits=1)+1
	rc.rows <- rowSums(div.mat>0) > rc.mid
	rc.midname <- paste("N",rc.mid,sep="")
	rc <- rarecurve(div.mat[rc.rows,], step = 10, sample = rc.mid, col = "blue", cex = 0.6, label=FALSE)
	rc.div <- as.numeric(unlist(lapply(rc,function(x) x[rc.midname])))
	rc.osd <- rownames(div.mat[rc.rows,])
	rc.df <- as.data.frame(cbind(rc.osd,rc.div))
	colnames(rc.df) <- c("label",rc.midname)
	rc.merge <- merge(osd.metadata[,c("label","temperature","salinity","latitude")],rc.df)
	
	rc.mergemid <- rc.merge[[eval(rc.midname)]]
	
	pdf(file=paste("fig",coldpar,my.taxa,"pdf",sep="."))
	rc.lat <- summary(lm(rc.mergemid ~ rc.merge$latitude))
	plot(rc.mergemid ~ rc.merge$latitude,main=my.taxa)
	abline(lm(rc.mergemid ~ rc.merge$latitude), main=my.taxa)

	rc.sal <- summary(lm(rc.mergemid ~ rc.merge$salinity))
	plot(rc.mergemid ~ rc.merge$salinity, main=my.taxa)
	abline(lm(rc.mergemid ~ rc.merge$salinity), main=my.taxa)

	rc.temp <- summary(lm(rc.mergemid ~ rc.merge$temperature))
	plot(rc.mergemid ~ rc.merge$temperature, main=my.taxa)
	abline(lm(rc.mergemid ~ rc.merge$temperature), main=my.taxa)
	dev.off()	
	next;
}

# merges taxonomy and protein parameters, calculates linear regressions against T, S, lat
#		my.protfam <- "IPR000719"; {
	for (my.protfam in many.protfams[,1]) {
		if (my.protfam=="") { next }
		for (coldpar in coldpar.list) {
			if (!file.exists(paste("figs",coldpar,my.protfam,sep="/"))) {
			dir.create(paste("figs",coldpar,my.protfam,sep="/"))
			}
		}

		aa.all <- lapply(osd.taxa.in, subset, interpro==my.protfam, select=coldpar.list)
		
		#skip if empty
		num.all <- sum(unlist(lapply(aa.all,function(x) dim(x)[1])))
		cat(num.protfams,"left / found",num.all,"peptides from",my.protfam," in ",my.taxa)
		num.protfams <- num.protfams - 1
		if ( num.all < 50 ) { 
			cat(" ... skipped\n")
			next
		}

		#summarize and plot
		aa.df <- data.frame()
		for (osd.label in names(osd.less.in)) {
		  osd.label.row <- which(osd.metadata$label == osd.label)
			aa.ln <- dim(aa.all[[osd.label]])[1]
			aa.df <- rbind( aa.df, cbind(
			"osd.label"=rep(osd.label,aa.ln),
			"latitude"=rep(osd.metadata$latitude[osd.label.row],aa.ln),
			"longitude"=rep(osd.metadata$longitude[osd.label.row],aa.ln),
			"temperature"=rep(osd.metadata$temperature[osd.label.row],aa.ln),
			"salinity"=rep(osd.metadata$salinity[osd.label.row],aa.ln),
			aa.all[[osd.label]]))
		  }
		  
#		if(dim(aa.df)[1] < 2) { next }

		all.p.vals <- data.frame()
		for (my.coldpar in coldpar.list) {
		#	my.coldpar <- "aliphatic_index"

		if(all(is.nan(aa.df[,my.coldpar]))) { next }
		if(sum(aa.df[,my.coldpar],na.rm=T)==0) { next }

#		all.lm <- lm(aa.df[,my.coldpar] ~ aa.df$latitude + aa.df$temperature + aa.df$salinity) # + aa.df$gc.pct)
#		summary.all <- summary.lm(all.lm)

		latitude.lm <- summary(lm(aa.df[,my.coldpar] ~ aa.df$latitude))
		if (dim(latitude.lm$coefficients)[1] < 2) { next }
		longitude.lm <- summary(lm(aa.df[,my.coldpar] ~ aa.df$longitude))
		if (dim(longitude.lm$coefficients)[1] < 2) { next }
		temperature.lm <- summary(lm(aa.df[,my.coldpar] ~ aa.df$temperature))
		if (dim(temperature.lm$coefficients)[1] < 2) { next }
		salinity.lm <- summary(lm(aa.df[,my.coldpar] ~ aa.df$salinity))
		if (dim(salinity.lm$coefficients)[1] < 2) { next }
	
		all.p.vals <- rbind(
			all.p.vals,
			cbind(
				"my.coldpar"=my.coldpar, 
				"my.taxa"=my.taxa,
				"my.protfam"=my.protfam,
				
				"latitude.intercept"=latitude.lm$coefficients[1,1],
				"latitude.slope"=latitude.lm$coefficients[2,1],
				"latitude.t"=latitude.lm$coefficients[2,3],
				"latitude.p"=latitude.lm$coefficients[2,4],	
				
				"longitude.intercept"=longitude.lm$coefficients[1,1],
				"longitude.slope"=longitude.lm$coefficients[2,1],
				"longitude.t"=longitude.lm$coefficients[2,3],
				"longitude.p"=longitude.lm$coefficients[2,4],	

				"temperature.intercept"=temperature.lm$coefficients[1,1],
				"temperature.slope"=temperature.lm$coefficients[2,1],
				"temperature.t"=temperature.lm$coefficients[2,3],
				"temperature.p"=temperature.lm$coefficients[2,4],
				
				"salinity.intercept"=salinity.lm$coefficients[1,1],
				"salinity.slope"=salinity.lm$coefficients[2,1],
				"salinity.t"=salinity.lm$coefficients[2,3],
				"salinity.p"=salinity.lm$coefficients[2,4]
				

#				"intercept.estimate"=summary.all$coefficients[1,1],
#				"intercept.t"=summary.all$coefficients[2,3],
#				"intercept.p"=summary.all$coefficients[2,4],
#,
#				"gc.pct.slope"=summary.all$coefficients[5,1],
#				"gc.pct.t"=summary.all$coefficients[5,3],
#				"gc.pct.p"=summary.all$coefficients[5,4]
				)
			)


  f.cols <- function(x) {
  if (x < 1e-4) { ur.col <- "#FF000066" } #blue
  else if (x < 1e-3) { ur.col <- "#FFEEB466" } #goldenrod
  else if (x < 1e-2) { ur.col <- "#FF8B3A66" } #red
  else { ur.col <- "#AAAAAA66" } #grey
  ur.col
	}

#plotting
		if (1 == 1) {
			pdf(file=paste(mydir,"figs/",my.coldpar,"/",my.protfam,"/","fig",".",my.taxa,".pdf",sep=""),width=10,height=10)
			par(mfrow=c(1,3))  

			#prevent boxplot not plotting
			aa.df[sapply(aa.df,is.nan)] <- 0

			boxplot(aa.df[,my.coldpar] ~ aa.df$latitude, col=f.cols(latitude.lm$coefficients[2,4]), at=sort(unique(aa.df$latitude)), pch=20,lty=0,xlab="latitude",ylab=my.coldpar)
			abline(coef(latitude.lm)[1:2,1])

			boxplot(aa.df[,my.coldpar] ~ aa.df$longitude, col=f.cols(longitude.lm$coefficients[2,4]), at=sort(unique(aa.df$longitude)), pch=20,lty=0,xlab="longitude",ylab=my.coldpar)
			abline(coef(longitude.lm)[1:2,1])

			boxplot(aa.df[,my.coldpar] ~ aa.df$temperature, col=f.cols(temperature.lm$coefficients[2,4]), at=sort(unique(aa.df$temperature)), pch=20,lty=0, xlab="temperature",ylab=my.coldpar)
			abline(coef(temperature.lm)[1:2,1])

			boxplot(aa.df[,my.coldpar] ~ aa.df$salinity, col=f.cols(salinity.lm$coefficients[2,4]), at=sort(unique(aa.df$salinity)), pch=20,lty=0, xlab="salinity",ylab=my.coldpar)
			abline(coef(salinity.lm)[1:2,1])
			
			dev.off()
		}
			
	} #end my.coldpar

	if (fhc < 1) {
		write.table(all.p.vals, file=fh, append = TRUE, sep = ",", dec=".", qmethod="double",col.names=T,row.names=F)
		fhc <- 1
	} else {
		write.table(all.p.vals, file=fh, append = TRUE, sep = ",", dec=".", qmethod="double",col.names=F,row.names=F)	
	}

	hitcount <- hitcount +1
	print(c(" ... ",hitcount))
	
	} #end my.protfam
	flush(fh)
    close(fh)
} #end my.taxa

print("complete")
