
if($#ARGV < 1){
  die "USAGE: $0 seq_file gff3_file offset[= 0] \n"
}

open $SEQ_FH, "<$ARGV[0]";

$null = <$SEQ_FH>;
while(<$SEQ_FH>){
  last if />/;
  chomp;
  $seq .= $_;
}

close $SEQ_FH;

open $TAB_FH, "<$ARGV[1]";

while(<$TAB_FH>){
  next if /^#/;
  chomp;
  @data = split "\t";
  $start = $data[3] - $ARGV[2];
  $end = $data[4]  + $ARGV[2];
  $out = substr($seq, $start - 1 , $end - $start + 1);
  $out =~ tr/ATGCatgc/TACGtacg/;
  $out = reverse($out);
  print ">SC3.182_2_rev_${start}-${end}_${data[8]}_homolog\n$out\n";
}

close $TAB_FH;
