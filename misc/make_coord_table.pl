
open $FA, "<$ARGV[0]";

my $seq_n = 0;
$null = <$FA>;

die "ERROR: Could not read fasta.\n" if not $null =~/^>/;

while(<$FA>){
  if(/^>/){
    $seq_n++;
  }else{
    chomp;
    $SEQ[$seq_n] .= $_;
  }
}

close $FA;

die "There are less or more than 2 sequenses.\n" if scalar @SEQ != 2;

my $seq1 =  $SEQ[0];
my $seq2 =  $SEQ[1];

if(length($seq1) != length($seq2)){
  die "Error:Two sequencies were different in size.\n"
}

my $n1 = 0;
while($seq1 =~ s/^(.)//){
  if($1 ne "-"){
    $n1++;
    push @SEQ1_COORD, $n1;
  }else{
    push @SEQ1_COORD, "-";
  }
}

my $n2 = 0;
while($seq2 =~ s/^(.)//){
  if($1 ne "-"){
    $n2++;
    push @SEQ2_COORD, $n2;
  }else{
    push @SEQ2_COORD, "-";
  }
}

while(@SEQ1_COORD){
 $c1 = shift @SEQ1_COORD;
 $c2 = shift @SEQ2_COORD;

 print "$c1\t.\t$c2\n";

}
