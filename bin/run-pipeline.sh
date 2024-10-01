#!/usr/bin/bash

#run kraken
mkdir step1-pfam.tax
cd step1-pfam.tax

#remove blank lines
for file in out.*.tax; do sed -i.bak '/^$/d' $file; done;

#tax file can have multiple instances of each peptide because different parts can match db
#tax file can have missing annotations
#to remove missing and remove duplicates
for file in out.*.tax; do sort -k1,2 -u $file | grep '\tk__' > tmp.$$; mv tmp.$$ $file; done;

cd ..

#run cold_info
mkdir step2-pfam.protparam
cd step2-pfam.protparam

ln -s ~/Dropbox/OSD/2016-02-15/step1-pfam.tax/*.tax ./
ln -s ~/Dropbox/OSD/2016-02-15/amino-acid-info.txt ./
for file in ./../step1-pfam.tax/*.tax; do ~/Dropbox/OSD/2016-02-15/bin/cold_info_pfam_no_bioperl.pl $file; done;
#for file in out.*.csv; do sort -k2 -u $file > tmp.$$; mv tmp.$$ $file; done;

cd ..

#run R merge
mkdir step3-merge
cd step3-merge

ln -s ~/Dropbox/OSD/2016-02-15/osd-metadata.csv ./
Rscript ~/Dropbox/OSD/2016-02-15/bin/merge-tax-protparam.r

#sort by most abundant
# cut -f15 *.tax | tr '|' '\n' | sort | uniq -c | sort -rn > count.taxa
# cut -f9 *.tax | tr '|' '\n' | sort | uniq -c | sort -rn > count.interpro

Rscript ~/Dropbox/OSD/2016-02-15/bin/select-summarize.r

#R osd.genus.in <- osd.all.in[["OSD1_2014-06-21_0m"]][c("read.part","read.full","acidic","aliphaticity","aliphatic_index","arg_lys","gravy","proline","mweight","interpro","tax.genus")]

cd ..
mkdir step4-heart
cd step4-heart

Rscript ~/Dropbox/OSD/2016-02-15/bin/plot-cold-alt.R
