#OSD9_2014-06-18_0m_NPL022.protparam.csv #this file holds the calculated parameters from cold_info_pfam.pl
#out.OSD9_2014-06-18_0m_NPL022.pfam.tax #this file holds the metadata for each gene
#osd-metadata.csv #this file holds the metadata for each site

#rm(list=ls())


args<-commandArgs(TRUE)

file.count.taxa <- args[1]
#file.count.taxa <- "list.k__"
#my.taxa <- "k__Bacteria"

outdir <- args[2]
#outdir <- "kraken-out-ensemblgenomes-80G"

curdir <- getwd()
mydir <- paste(curdir,"/",outdir,"/",sep="")
tmpdir <- ./ #"/Volumes/ramdisk/"

osd.all.in <- list()

osd.all.in.saved <- "" #paste(tmpdir,"20150811.osd.all.in.Rdata",sep="")
osd.metadata <- read.csv("osd-metadata.csv",sep="\t")

if(file.exists(osd.all.in.saved)) {

	#load previously saved table from osd.all.in.Rdata
	cat("loading saved data from ",osd.all.in.saved,"\n")
	load(osd.all.in.saved)

	for (osd.label in names(osd.all.in)) {
		osd.all.in[[osd.label]]$tax.species <- gsub('([[:punct:]])|\\s+','_',osd.all.in[[osd.label]]$tax.species)
	}

} else {

#do this once, then just load file above
#read data into dataframes
	print("reading in data")
	for (osd.label in osd.metadata$label) {
		cat(osd.label)
		osd.cold.file <- paste(mydir,osd.label,"_NPL022.protparam.csv",sep="")
		if( ! file.exists( osd.cold.file )) { cat("\tcouldn't find file",osd.cold.file,"\n"); next }
		osd.tax.file <- paste(mydir,"out.",osd.label,"_NPL022.pfam.tax",sep="")
		if( ! file.exists( osd.tax.file )) { cat("\tcouldn't find file ",osd.tax.file,"\n"); next }
		osd.gcs.file <- paste(osd.label,"_NPL022_readsWithMatches.gz.gcs",sep="")
		if( ! file.exists( osd.gcs.file)) { cat("couldn't find file ",osd.gcs.file,"\n"); next }
			
		osd.cold.in <- read.csv( osd.cold.file, sep="\t", row.names = NULL)
		tmp.name <- data.frame(do.call('rbind', strsplit(as.character(osd.cold.in$read_id),'_',fixed=TRUE))) #what's this for?
		osd.cold.in <- cbind(tmp.name,osd.cold.in)
		osd.cold.in <- osd.cold.in[,-c(2,3,4,5)]
		colnames(osd.cold.in)[1:2] <- c("read.part","read.full")
		osd.cold.in$OSD.label <- rep(osd.label,dim(osd.cold.in)[1])

		tmp.gc <- read.csv(osd.gcs.file,sep="\t",header=F)
		colnames(tmp.gc) <- c("read.part","gc.pct")

		osd.cold.in <- merge(osd.cold.in,tmp.gc,by="read.part")	
		osd.tax.in <- read.csv( osd.tax.file, sep="\t", row.names = NULL, header=F)
		#only works with my edited kraken
		tmp.tax <- data.frame(do.call('rbind', strsplit(as.character(osd.tax.in$V15),'|',fixed=TRUE)))
		tmp.name <- data.frame(paste(as.character(osd.tax.in$V1),as.character(osd.tax.in$V2),sep=""))

		osd.tax.in <- cbind(osd.tax.in,tmp.name,tmp.tax)

		osd.tax.in <- osd.tax.in[,-c(1,2,3,4,5,6,13,14,15)]
		colnames(osd.tax.in) <- c("pfam","pfam.def","interpro","interpro.def","other.db","other.db.def","read.full","tax.domain","tax.phylum","tax.class","tax.order","tax.family","tax.genus","tax.species")
		osd.all.in[[osd.label]] <- merge(osd.cold.in,osd.tax.in,by="read.full")
		osd.all.in[[osd.label]]$tax.species <- gsub('([[:punct:]])|\\s+','_',osd.all.in[[osd.label]]$tax.species)

		cat("\n")
	}


	save(osd.all.in,file=osd.all.in.saved)


}

#access as osd.all.in[["OSD99_2014-06-21_1m"]]


###look for protein and taxa of interest

#sort by most abundant
# cut -f15 *.tax | tr '|' '\n' | sort | uniq -c | sort -rn > count.taxa
# cut -f9 *.tax | tr '|' '\n' | sort | uniq -c | sort -rn > count.interpro

all.uniq.taxa <- read.csv(file.count.taxa,header=F,sep="~")
all.uniq.taxa <- sapply(all.uniq.taxa,function(x) gsub('^\\s*[[:digit:]]+\\s+','',x))
#all.uniq.taxa <- sapply(all.uniq.taxa,function(x) gsub('([[:punct:]])|\\s+','_',x))
#all.uniq.interpro <- read.csv("count3.interpro",header=F,sep=" ")

coldpar.list <- c("acidic","aliphaticity", "aliphatic_index", "arg_lys","gravy","proline","gc.pct")
#	my.taxa <- "k__Bacteria"
for (my.taxa in as.character(all.uniq.taxa[,1])) {
	my.outfile <- paste(mydir,"pvals.new.",my.taxa,".csv",sep="")
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
	for (my.protfam in taxa.protfams) {
		if (my.protfam=="") { next }
#		my.protfam <- "IPR000456"

		aa.all <- lapply(osd.taxa.in, subset, interpro==my.protfam, select=coldpar.list)
		
		#skip if empty
		num.all <- sum(unlist(lapply(aa.all,function(x) dim(x)[1])))
		cat(num.protfams,"\tleft / found\t",num.all,"\tproteins from\t",my.protfam," in ",my.taxa)
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

		all.lm <- lm(aa.df[,my.coldpar] ~ aa.df$latitude + aa.df$temperature + aa.df$salinity + aa.df$gc.pct)
		summary.all <- summary.lm(all.lm)
		all.p.vals <- rbind(
			all.p.vals,cbind(
				"my.coldpar"=my.coldpar, 
				"my.taxa"=my.taxa,
				"my.protfam"=my.protfam,
				"intercept.estimate"=summary.all$coefficients[1,1],
				"intercept.t"=summary.all$coefficients[1,3],
				"intercept.p"=summary.all$coefficients[1,4],
				"latitude.slope"=summary.all$coefficients[2,1],
				"latitude.t"=summary.all$coefficients[2,3],
				"latitude.p"=summary.all$coefficients[2,4],				
				"temperature.slope"=summary.all$coefficients[3,1],
				"temperature.t"=summary.all$coefficients[3,3],
				"temperature.p"=summary.all$coefficients[3,4],
				"salinity.slope"=summary.all$coefficients[4,1],
				"salinity.t"=summary.all$coefficients[4,3],
				"salinity.p"=summary.all$coefficients[4,4],
				"gc.pct.slope"=summary.all$coefficients[5,1],
				"gc.pct.t"=summary.all$coefficients[5,3],
				"gc.pct.p"=summary.all$coefficients[5,4]
				)
			)
			
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





#for plotting

#f.cols <- function(x) {
#  if (x < 0.0001) { ur.col <- "#FF0000FF" }
#  else if (x < 0.01) { ur.col <- "#BB0000BB" }
#  else if (x < 0.05) { ur.col <- "#88000088" }
#  else { ur.col <- "#44000044" }
#  ur.col
#}

#plotting
#			pdf(file=paste(mydir,"fig",".",my.coldpar,".",my.taxa,".",my.protfam,".pdf",sep=""),width=10,height=10)
#			par(mfrow=c(1,3))  
#			summary(latitude.lm)
#			plot(aa.df[,my.coldpar] ~ aa.df$temperature, col=f.cols(lmp(latitude.lm)) )
#			if(latitude.lmp!=1) abline(latitude.lm)
#			summary(temperature)
#			plot(aa.df[,my.coldpar] ~ aa.df$temperature, col=f.cols(lmp(temperature.lm)) )
#			if(temperature.lmp!=1) abline(temperature.lm)
#			summary(salinity.lm)
#			plot(aa.df[,my.coldpar] ~ aa.df$salinity, col=f.cols(lmp(salinity.lm)) )
#			if(salinity.lmp!=1) abline(salinity.lm)
#			dev.off()


#some functions
#lmp returns the p value for an lm object

#lmp <- function (modelobject) {
#    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
#    p <- list()
#    for (i in length(modelobject$coefficients)) {
#		if (is.na(modelobject$coefficients[i])) {warning("Univariate input "); return(1)}
#		f <- summary(modelobject)$fstatistic
#		p[i] <- pf(f[1],f[2],f[3],lower.tail=F)
#		attributes(p[i]) <- NULL
#   }
#    return(p)
#}
