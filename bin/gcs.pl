while(<>) {
if ($_ =~ s/^\>//) {chomp; print $_ . "\t"} else { $x=length($_); $y=tr/[GC]//; print $y/$x . "\n"}
}
