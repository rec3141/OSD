OUTDIR=/work/rec3141/OSD-analysis/unipept_tax
GZFILE=$1

#to run in Parallel
#for GZFILE in `ls out.OSD*.pfam.gz`; do
	PFAMFILE=`basename $GZFILE .gz`
	if [ -e $OUTDIR/$PFAMFILE.tax ]; then
 	echo "$PFAMFILE.tax exists"
	else
	touch $OUTDIR/$PFAMFILE.tax
	echo "processing $GZFILE"
	OUTFILE=`basename $( echo $GZFILE ) .pfam.gz | cut -f2-3 -d'.' | cut -f1-4 -d'_'`
	gunzip -c $GZFILE | sed -E 's#(_[0-9]+_[0-9]+_[+-])#	\1#g' | sort -b > $OUTDIR/$PFAMFILE.sort
	echo "writing $OUTDIR/$PFAMFILE.tax"
#kraken	join -1 1 -2 1 -t'	' $OUTDIR/$PFAMFILE.sort $OUTDIR/kraken.output."$OUTFILE"_readsWithMatches.tax.sort | sort -k7 > $OUTDIR/$PFAMFILE.tax
	join -1 1 -2 1 -t'	' $OUTDIR/$PFAMFILE.sort $OUTDIR/$OUTFILE.tax | sort -k7 > $OUTDIR/$PFAMFILE.tax
	rm $OUTDIR/$PFAMFILE.sort
	fi
#done
