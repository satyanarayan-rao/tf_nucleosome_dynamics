# Convert a bed file to a wig file for paired end sequencing : fill in the reads
#	Extend midpoint of the reads to 30 bp on each side - assuming the central region of a nucleosome
# For calling nucleosome positions - use only nucleosome sized fragments.


BEGIN { push @INC, '/beevol/home/zukowski/ER_CUTnRUN/code' }

use genome_size;
use ngs;

die "Usage: perl bed2wig.pl <BED FILE> <WIG FILE NAME> <MIN> <MAX>\n" if(!$ARGV[3]);


$min=$ARGV[2];
$max=$ARGV[3];

$genome = 'hg38';
if($ARGV[4]){
	$genome = $ARGV[4];
}

$tgs = &genome_size::getSize($genome);
%genome_size = %{$tgs};

open(LIST,$ARGV[0]) || die "INPUT $!\n";
while(chomp($l=<LIST>)){
	open(FILE,$l) || die "Bed file: $l $!\n";
	while(chomp($line=<FILE>)){
		$lno++;
		@temp = split /[\ \s\n\t]+/, $line;
		$frag_length = $temp[2]-$temp[1];
		#Assuming 6 column bed file - change if different
		$temp[0]=~s/chr//;
		if($#temp != 5){
			print STDERR "Not regular BED line?\n$line\n";
		}elsif($frag_length>=$min && $frag_length<=$max && exists($genome_size{$temp[0]})){
			$mp = int( ($temp[1] + $temp[2])/20 + 0.5 );
			$lower=$mp-2;
			$upper=$mp+2;
			for($i=$lower;$i<=$upper;$i++){
				if($i<= int($genome_size{$temp[0]}/10)){
					$read{$temp[0]}{$i}++;
					$count++;
				}
			}
		}
		print STDERR "Count:$count\n" if($count%1000000==0);
		print STDERR "Line No:$lno\n" if($lno%10000==0);
	}
}
print STDERR "Finished reading bed file\nCount=$count\n";
close(FILE);

my $GN = 139712364; # Genome size
if($genome eq 'mm10'){
	$GN = 2800000000;
}elsif($genome eq 'hg38'){
	$GN = 3300000000;
}elsif($genome eq 'sacCer3'){
	$GN = 12157105;
}

open(OUT,">$ARGV[1]") || die "OUT $!\n";
print OUT "track type=wiggle_0\n";
foreach $i (keys (%read) ){
	print OUT "variableStep  chrom=chr$i span=10\n";
	%thash = %{$read{$i}};
	for($j=0;$j<=int($genome_size{$i}/10);$j++){
		if(exists $thash{$j}){
			$normval = $thash{$j}*$GN/$count ;		
			print OUT ($j*10+1)." $normval\n";
		}
	}
}
close(OUT);
