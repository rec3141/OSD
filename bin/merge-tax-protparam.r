#R script
# run as Rscript ./merge-tax-protparam.r ./step3-merge

#OSD9_2014-06-18_0m_NPL022.protparam.csv #this file holds the calculated parameters from cold_info_pfam.pl
#out.OSD9_2014-06-18_0m_NPL022.pfam.tax #this file holds the metadata for each gene
#osd-metadata.csv #this file holds the metadata for each site

#rm(list=ls())

#args<-commandArgs(TRUE)
#indir <- args[1] #where the data is found

dir1 <- "./../step1-pfam.tax/"
dir2 <- "./../step2-pfam.protparam/"

osd.all.in <- list()
osd.less.in <- list()
osd.all.in.saved <- "osd.all.in.saved.Rdata"
osd.less.in.saved <- "osd.less.in.saved.Rdata"

osd.metadata <- read.csv("./osd-metadata.csv",sep="\t")

#read data into dataframes
	print("reading in data")
	for (osd.label in osd.metadata$label) {
		cat(osd.label)

		osd.tax.file <- paste(dir1, "out.",osd.label,"_NPL022.pfam.tax",sep="")
		if( ! file.exists( osd.tax.file )) { cat("\tcouldn't find file ",osd.tax.file,"\n"); next }

		osd.cold.file <- paste(dir2,"out.",osd.label,"_NPL022.protparam.csv",sep="")
		if( ! file.exists( osd.cold.file )) { cat("\tcouldn't find file",osd.cold.file,"\n"); next }

#		osd.gcs.file <- paste(osd.label,"_NPL022_readsWithMatches.gz.gcs",sep="")
#		if( ! file.exists( osd.gcs.file)) { cat("couldn't find file ",osd.gcs.file,"\n"); next }
			
		osd.cold.in <- read.csv( osd.cold.file, sep="\t", row.names = NULL, header=F)
#		colnames(osd.cold.in) <- c("method","read_id","length","acidic","aliphaticity","aliphatic_index","arg_lys","gravy","proline","mweight")
#		method	read_id	length	acidic	aliphaticity	aliphatic_index	arg_lys	gravy	polar	nonpolar	sulfur	hbond	basic	aromatic	A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	nitrogen	mweight
		colnames(osd.cold.in) <- c("method","read_id","length","acidic","aliphaticity","aliphatic_index","arg_lys","gravy","polar","nonpolar","sulfur","hbond","basic","aromatic","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","nitrogen","mweight")
		
		tmp.name <- data.frame(do.call('rbind', strsplit(as.character(osd.cold.in$read_id),'_',fixed=TRUE))) #whats this for?
		osd.cold.in <- cbind(tmp.name,osd.cold.in)
		osd.cold.in <- osd.cold.in[,-c(2,3,4,5)]
		colnames(osd.cold.in)[1:2] <- c("read.part","read.full")
		osd.cold.in$OSD.label <- rep(osd.label,dim(osd.cold.in)[1])

#		tmp.gc <- read.csv(osd.gcs.file,sep="\t",header=F)
#		colnames(tmp.gc) <- c("read.part","gc.pct")

#		osd.cold.in <- merge(osd.cold.in,tmp.gc,by="read.part")	
		osd.tax.in <- read.csv( osd.tax.file, sep="\t", row.names = NULL, header=F)
		
		#only works with my edited kraken
		tmp.tax <- data.frame(do.call('rbind', strsplit(as.character(osd.tax.in$V15),'|',fixed=TRUE)))
		tmp.name <- data.frame(paste(as.character(osd.tax.in$V1),as.character(osd.tax.in$V2),sep=""))

		osd.tax.in <- cbind(osd.tax.in,tmp.name,tmp.tax)

		osd.tax.in <- osd.tax.in[,-c(1,2,3,4,5,6,13,14,15)]
		colnames(osd.tax.in) <- c("pfam","pfam.def","interpro","interpro.def","other.db","other.db.def","read.full","tax.domain","tax.phylum","tax.class","tax.order","tax.family","tax.genus","tax.species")
#they should be comparable/sorted
#		osd.all.in[[osd.label]] <- merge(osd.cold.in,osd.tax.in,by="read.full")
		osd.all.in[[osd.label]] <- cbind(osd.cold.in,osd.tax.in)
		osd.all.in[[osd.label]]$tax.species <- gsub('([[:punct:]])|\\s+','_',osd.all.in[[osd.label]]$tax.species)

		osd.less.in[[osd.label]] <- osd.all.in[[osd.label]][,c("read.part","length","acidic","aliphaticity","aliphatic_index","arg_lys","gravy","polar","nonpolar","sulfur","hbond","basic","aromatic","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","nitrogen","mweight","OSD.label","interpro","tax.domain","tax.phylum", "tax.class","tax.order","tax.family","tax.genus","tax.species")]
	#	osd.genus.in <- osd.all.in[[osd.label]][c("read.part","read.full","acidic","aliphaticity","aliphatic_index","arg_lys","gravy","proline","mweight","interpro","tax.genus")]
		cat("\n")
	}

	print("writing osd.all.in file")
	save(osd.all.in,file=osd.all.in.saved)
	Sys.chmod(osd.all.in.saved, mode = "0444")
	print("writing osd.less.in file")
	save(osd.less.in,file=osd.less.in.saved)
	Sys.chmod(osd.less.in.saved, mode = "0444")


	

