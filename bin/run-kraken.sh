#!/bin/bash


KDB="/Volumes/ramdisk/"
INDIR="/work/rec3141/seq_data/OSD2014"
OUTDIR="kraken-out-ensemblgenomes-80G"
mkdir $OUTDIR

#su admin; sudo diskutil erasevolume HFS+ 'ramdisk' `hdiutil attach -nomount ram://262144000`
#su admin; sudo cp -a $KDB /Volumes/ramdisk/
#exit

for INPUTFILE in `find $INDIR -name "OSD*readsWithMatches.gz"`; do
 OUTFILE=`basename $INPUTFILE .gz`
 if [ -e $OUTDIR/kraken.taxonomy.$OUTFILE.sort ]; then echo "$OUTFILE exists"; else
 echo "running $OUTFILE"
 #kraken --preload --db ~/Applications/kraken/minikraken_20141208 --threads 16 --classified-out kraken.classified.$INPUTFILE --unclassified-out kraken.unclassified.$INPUTFILE --output kraken.output.$INPUTFILE --fastq-input $INPUTFILE
 #first run on minikraken# kraken --preload --db ~/Applications/kraken/minikraken_20141208 --threads 15 --classified-out $OUTDIR/kraken.classified.$OUTFILE --only-classified-output --output $OUTDIR/kraken.output.$OUTFILE  $INPUTFILE
# kraken --db $KDB --threads 14 --classified-out $OUTDIR/kraken.classified.$OUTFILE --only-classified-output --output $OUTDIR/kraken.output.$OUTFILE  $INPUTFILE
 kraken-translate-rec --db $KDB --complete-format $OUTDIR/kraken.output.$OUTFILE | sort -b > $OUTDIR/kraken.taxonomy.$OUTFILE.sort
 fi
done;


find ./kraken-out-ensemblgenomes-80G -name "kraken.output.*" | 
  parallel "echo {};
            kraken-translate-rec --db /work/rec3141/apps/kraken/ensemblgenomes-80G --complete-format {} | sort -b > {}.tax.sort"

#            OUTDIR=`dirname {}`;
#            echo $OUTDIR;
#            OUTFILE=`basename {} .gz | cut -f3-9 -d'.'`;
#            echo $OUTFILE;

#FYI
#make a ramdisk at /Volumes/ramdisk
#su admin; sudo diskutil erasevolume HFS+ 'ramdisk' `hdiutil attach -nomount ram://262144000`
#Started erase on disk8
#Unmounting disk
#Erasing
#Initialized /dev/rdisk8 as a 125 GB case-insensitive HFS Plus volume
#Mounting disk
#Finished erase on disk8 ramdisk

#ln $KDB /Volumes/ramdisk/
