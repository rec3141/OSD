save.ttest <- read.csv("HI_110001.1e-20_.8_2_1.cluster_list.ttest.csv",row.names=1,header=T)
COG.letter <- unlist(lapply(save.ttest$COG.function, function(x) unlist(strsplit(as.character(x)," ",fixed=T)[[1]][2])))
save.ttest <- cbind(save.ttest,COG.letter)

COG.rank.save <- NULL

pdf(file="COG.ranks.pdf",height=8,width=8)
#acidic
subdf <- subset(save.ttest, nchar(as.character(save.ttest$COG.letter))==1)
subdf$COG.letter <- factor(subdf$COG.letter)
COG.med <- aggregate(rank(acidic) ~ COG.letter, FUN=function(x) median(x), data=subdf)
COG.rank.save <- rank(COG.med[,2])
COG.med <- COG.med[ order(COG.med[,2]) , ]
subdf$COG.letter <- factor(subdf$COG.letter,COG.med[,1])
boxplot(rank(acidic) ~ COG.letter, data=subdf,main="acidic")
boxplot(-log(acidic) ~ COG.letter, data=subdf,main="acidic",ylim=c(0,8))

#arg_lys
subdf <- subset(save.ttest, nchar(as.character(save.ttest$COG.letter))==1)
subdf$COG.letter <- factor(subdf$COG.letter)
COG.med <- aggregate(rank(arg_lys) ~ COG.letter, FUN=function(x) median(x), data=subdf)
COG.rank.save <- cbind(COG.rank.save,rank(COG.med[,2]))
COG.med <- COG.med[ order(COG.med[,2]) , ]
subdf$COG.letter <- factor(subdf$COG.letter,COG.med[,1])
boxplot(rank(arg_lys) ~ COG.letter, data=subdf,main="arg_lys")
boxplot(-log(arg_lys) ~ COG.letter, data=subdf,main="arg_lys",ylim=c(0,10))

#aliphaticity
subdf <- subset(save.ttest, nchar(as.character(save.ttest$COG.letter))==1)
subdf$COG.letter <- factor(subdf$COG.letter)
COG.med <- aggregate(rank(aliphaticity) ~ COG.letter, FUN=function(x) median(x), data=subdf)
COG.rank.save <- cbind(COG.rank.save,rank(COG.med[,2]))
COG.med <- COG.med[ order(COG.med[,2]) , ]
subdf$COG.letter <- factor(subdf$COG.letter,COG.med[,1])
boxplot(rank(aliphaticity) ~ COG.letter, data=subdf,main="aliphaticity")
boxplot(-log(aliphaticity) ~ COG.letter, data=subdf,main="aliphaticity",ylim=c(0,10))


#aliphatic_index
subdf <- subset(save.ttest, nchar(as.character(save.ttest$COG.letter))==1)
subdf$COG.letter <- factor(subdf$COG.letter)
COG.med <- aggregate(rank(aliphatic_index) ~ COG.letter, FUN=function(x) median(x), data=subdf)
COG.rank.save <- cbind(COG.rank.save,rank(COG.med[,2]))
COG.med <- COG.med[ order(COG.med[,2]) , ]
subdf$COG.letter <- factor(subdf$COG.letter,COG.med[,1])
boxplot(rank(aliphatic_index) ~ COG.letter, data=subdf,main="aliphatic_index")
boxplot(-log(aliphatic_index) ~ COG.letter, data=subdf,main="aliphatic_index",ylim=c(0,10))

#proline
subdf <- subset(save.ttest, nchar(as.character(save.ttest$COG.letter))==1)
subdf$COG.letter <- factor(subdf$COG.letter)
COG.med <- aggregate(rank(proline) ~ COG.letter, FUN=function(x) median(x), data=subdf)
COG.rank.save <- cbind(COG.rank.save,rank(COG.med[,2]))
COG.med <- COG.med[ order(COG.med[,2]) , ]
subdf$COG.letter <- factor(subdf$COG.letter,COG.med[,1])
boxplot(rank(proline) ~ COG.letter, data=subdf,main="proline")
boxplot(-log(proline) ~ COG.letter, data=subdf,main="proline",ylim=c(0,20))

#gravy
subdf <- subset(save.ttest, nchar(as.character(save.ttest$COG.letter))==1)
subdf$COG.letter <- factor(subdf$COG.letter)
COG.med <- aggregate(rank(gravy) ~ COG.letter, FUN=function(x) median(x), data=subdf)
COG.rank.save <- cbind(COG.rank.save,rank(COG.med[,2]))
COG.med <- COG.med[ order(COG.med[,2]) , ]
subdf$COG.letter <- factor(subdf$COG.letter,COG.med[,1])
boxplot(rank(gravy) ~ COG.letter, data=subdf,main="gravy")
boxplot(-log(gravy) ~ COG.letter, data=subdf,main="gravy",ylim=c(0,4))

dev.off()
