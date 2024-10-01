#OSD9_2014-06-18_0m_NPL022.protparam.csv #this file holds the calculated parameters from cold_info_pfam.pl
#out.OSD9_2014-06-18_0m_NPL022.pfam.tax #this file holds the metadata for each gene
#osd-metadata.csv #this file holds the metadata for each site

#rm(list=ls())

args<-commandArgs(TRUE)

file.count.taxa <- args[1] #count3.taxa
file.count.interpro <- args[2] #count3.interpro
outdir <- args[3]

file.count.taxa <- "./../count3.order"
file.count.interpro <- "./../count3.interpro"
outdir <- "./"

curdir <- getwd()
mydir <- paste(curdir,"/",outdir,"/",sep="")

#osd.all.in.saved <- "osd.small.in.saved.Rdata"
osd.all.in.saved <- "osd.all.in.saved.Rdata"
osd.metadata <- read.csv("./../osd-metadata.csv",sep="\t")


#if run interactively with source(), don't reload the whole dataset
if(!exists("osd.all.in")) {
	osd.all.in <- list()

	if(file.exists(osd.all.in.saved)) {

	#load previously saved table from osd.all.in.saved.Rdata
	cat("loading saved data from",osd.all.in.saved,"\n")
	load(osd.all.in.saved)
	} else {
	cat("error\n")
	quit()
}
}

#access as osd.all.in[["OSD99_2014-06-21_1m"]]

###look for protein and taxa of interest

all.uniq.taxa <- read.csv(file.count.taxa,header=F,sep="~")
all.uniq.taxa <- sapply(all.uniq.taxa,function(x) gsub('^\\s*[[:digit:]]+\\s+','',x))
#all.uniq.taxa <- sapply(all.uniq.taxa,function(x) gsub('([[:punct:]])|\\s+','_',x))
all.uniq.interpro <- read.csv(file.count.interpro,header=F,sep=" ",skip=1)
many.uniq.interpro <- all.uniq.interpro[all.uniq.interpro[,1]>=100,]
many.uniq.interpro <- many.uniq.interpro[order(many.uniq.interpro[,2]),]

coldpar.list <- c("acidic","aliphaticity", "aliphatic_index", "arg_lys","gravy","proline","mweight") # ,"gc.pct"

#		my.protfam <- "IPR000719"; {

	for (my.protfam in many.uniq.interpro) {
		if (my.protfam=="") { next }


		#	my.taxa <- "k__Bacteria"
	  for (my.taxa in as.character(all.uniq.taxa[,1])) {
		my.outfile <- paste(mydir,"out.pvals.",my.taxa,".csv",sep="")
		if (length(grep("__$",my.taxa)) > 0) { next }
		if (file.exists(my.outfile)) { next }

	#	file.remove(my.outfile)
		fh <- file(my.outfile,open="at")
		fhc <- 0

		cat("selecting ",my.taxa,"\n")
		osd.taxa.in <- lapply(osd.all.in, subset, tax.domain==my.taxa | tax.phylum==my.taxa | tax.class==my.taxa | tax.order==my.taxa | tax.family==my.taxa | tax.genus==my.taxa | tax.species==my.taxa, select=c("interpro",coldpar.list))

		taxa.protfams <- sort(as.character(unique(unlist(lapply(osd.taxa.in,function(x) as.character(unique(x$interpro)))))))
		num.protfams <- length(taxa.protfams)
		cat("found\t",num.protfams, "\tprotein families in\t",my.taxa,"\n")

		many.protfams <-  merge(taxa.protfams,many.uniq.interpro,by.x=1,by.y="V2")
		many.protfams <- many.protfams[order(-many.protfams[,2]),]
		num.protfams <- length(many.protfams[,1])
		cat("analyzing\t",num.protfams, "\tprotein families in\t",my.taxa,"\n")
	
		aa.all <- lapply(osd.taxa.in, subset, interpro==my.protfam, select=coldpar.list)
		
		#skip if empty
		num.all <- sum(unlist(lapply(aa.all,function(x) dim(x)[1])))
		cat(num.protfams,"left / found",num.all,"proteins from",my.protfam," in ",my.taxa)
		num.protfams <- num.protfams - 1
		if ( num.all < 50 ) { 
			cat(" ... skipped\n")
			next
		}

		#summarize and plot
		aa.df <- data.frame()
		for (osd.label in names(osd.all.in)) {
		  osd.label.row <- which(osd.metadata$label == osd.label)
			aa.ln <- dim(aa.all[[osd.label]])[1]
			aa.df <- rbind( aa.df, cbind(
			"osd.label"=rep(osd.label,aa.ln),
			"latitude"=rep(osd.metadata$latitude[osd.label.row],aa.ln),
			"temperature"=rep(osd.metadata$temperature[osd.label.row],aa.ln),
			"salinity"=rep(osd.metadata$salinity[osd.label.row],aa.ln),
			aa.all[[osd.label]]))
		  }
		  
#		if(dim(aa.df)[1] < 2) { next }

		all.p.vals <- data.frame()
		for (my.coldpar in coldpar.list) {
		#	my.coldpar <- "aliphatic_index"

		if(all(is.nan(aa.df[,my.coldpar]))) { next }

#		all.lm <- lm(aa.df[,my.coldpar] ~ aa.df$latitude + aa.df$temperature + aa.df$salinity) # + aa.df$gc.pct)
#		summary.all <- summary.lm(all.lm)

		latitude.lm <- summary(lm(aa.df[,my.coldpar] ~ aa.df$latitude))
		if (dim(latitude.lm$coefficients)[1] < 2) { next }
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
			pdf(file=paste(mydir,"fig",".",my.coldpar,".",my.taxa,".",my.protfam,".pdf",sep=""),width=10,height=10)
			par(mfrow=c(1,3))  

			#prevent boxplot not plotting
			aa.df[sapply(aa.df,is.nan)] <- 0

			boxplot(aa.df[,my.coldpar] ~ aa.df$latitude, col=f.cols(latitude.lm$coefficients[2,4]), at=sort(unique(aa.df$latitude)), pch=20,lty=0,xlab="latitude",ylab=my.coldpar)
			abline(coef(latitude.lm)[1:2,1])

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

	cat(" ... done\n")
	
	} #end my.protfam
	flush(fh)
    close(fh)
} #end my.taxa

print("complete")
