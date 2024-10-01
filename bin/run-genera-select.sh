#awk -F'\t' '{print $1 "\t" $15}' kraken-out-ensemblgenomes-80G/out.*.pfam.tax | tr ' ' '_' > tmp.tax
#awk -F'\t' '{print $1 "\t" $15}' unipept_tax/out.*.pfam.tax | tr ' ' '_' >> tmp.tax
#sort tmp.tax > tmp.sort
#cut -f1-6 -d'|' tmp.sort | grep -v '_$' | cut -f2 | sort | uniq -d > tmp.genera

#ls */out.*.pfam.tax | parallel --gnu "grep -f tmp.genera {} > {}.genera"
#ls */out.*.genera | parallel --gnu "perl ./cold_info_pfam.pl {}"
#Rscript --vanilla osd-aa.R list.s__ kraken-out-ensemblgenomes-80G
ls kraken-out-ensemblgenomes-80G/list* | parallel --gnu "Rscript --vanilla osd-aa.R {} kraken-out-ensemblgenomes-80G"

mail -s "osd-aa.R done" rec3141@gmail.com <<EOF
`date`
EOF
