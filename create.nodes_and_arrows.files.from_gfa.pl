#!/usr/bin/perl
#
# Apologies for any clunkiness in this code
#
my $homeDir = "/Users/srmcgee/graph.map.Shiny";

my $gfa = "$homeDir/roi17q21.31.minigraph.baseAln.gfa";

my $out = "$homeDir/gfa.for_ggplot.txt";
open(OUT,">$out") || die "can't open out file: $out\n";
print OUT "id\trank\tsubrank\tstart\tstop\tsize\n";


my $aout = "$homeDir/gfa.arrows.for_ggplot.txt";
open(AOUT,">$aout") || die "can't open out file: $aout\n";
print AOUT "id\trank\tx1\ty1\tx2\ty2\n";

my $lengthTag = "LN:i:";
my $nameTag = "SN:Z:";
my $offsetTag = "SO:i:";
my $rankTag = "SR:i:";

my %rankHash = ();
my %startStopByRankHash = ();
my %scaledStartStopHash = ();
my %positionTranslationHash = (); # true position to scaled position
my $gapSize = 1;

my $totalRelength = 0;
my %idRankHash = ();
my $lastStop = 0;
my $lastRank = 0;

open(GFA,$gfa) || die "can't open GFA: $gfa\n";
while(<GFA>) {
	my $line = "$_";
	chomp($line);
	my @lineArray = split(/\t/,$line);
	if ($lineArray[0] eq "S") {

		$id = $lineArray[1];
		$length = $lineArray[3];
		$name = $lineArray[4];
		$offset = $lineArray[5];
		$rank = $lineArray[6];
		$length =~ s/$lengthTag//;
		$name =~ s/$nameTag//;
		$offset =~ s/$offsetTag//;
		#$offset += 1; # to deal with log 
		$rank =~ s/$rankTag//;

		$relength = resize($length);
		$loglength = sprintf("%6.2f",log($length));
		if ($loglength == 0) {
			$loglength = sprintf("%6.2f",0.5);;
		}
		

		if ($rank == 0) {
			#
			# first get the horizontal scale from summing up lengths and offset
			#  and scale each section by size and units of graph
			#  -- this assumes a sorted file
			#			
			for (my $i = $offset; $i < ($offset + $length); $i++) {
				my $scale = 0;
				if ($length > 0) { # this should never happen
					$scale = $relength/$length;
				} else {
					print STDOUT "*** segment length is 0!\n";
				}
				my $repos = int(($i - $offset)*$scale) + $totalRelength;
				#print STDOUT "$i\t$scale\t$repos\n";
				
				$positionTranslationHash{$i} = $repos;
			}
			
			print STDOUT "$rank\t$offset\t$relength\t$totalRelength\t$positionTranslationHash{$offset}\n";
			$totalRelength += $relength+$gapSize;
			
			$startStopByRankHash{$rank}{$id}{'start'} = $positionTranslationHash{$offset};
			$startStopByRankHash{$rank}{$id}{'size'} = $relength;
			$startStopByRankHash{$rank}{$id}{'stop'} = $startStopByRankHash{$rank}{$id}{'start'} + $relength;
			$startStopByRankHash{$rank}{$id}{'subrank'} = 0;
			$idRankHash{$id} = $rank;
			
		} else {
			
			
			$startStopByRankHash{$rank}{$id}{'start'} = $positionTranslationHash{$offset};
			$startStopByRankHash{$rank}{$id}{'size'} = $relength;
			$startStopByRankHash{$rank}{$id}{'stop'} = $startStopByRankHash{$rank}{$id}{'start'} + $relength;
			print STDOUT "$rank\t$offset\t$relength\t$totalRelength\t$positionTranslationHash{$offset}\n";
			
			#
			# alternate rank>0 segments above and below stable sequence
			#
			my $subrank = $rank;
			if ($rank%2 == 0) {
				#$rankHash{$rank} = 
				$subrank = -$rank/2;
			} else {
				$subrank = 1+($rank-1)/2;
			}
			#
			# watch out for overlapping segments
			#
			$rankLoc = $subrank;
			if ($lastRank == $rank){
				if ($startStopByRankHash{$rank}{$id}{'start'} <= $lastStop) {
					$rankHash{$rank}++;
					$rankLoc = $subrank + 0.2*$rankHash{$rank};
				} 
			} else {
				$lastRank = 0;
			}
			$startStopByRankHash{$rank}{$id}{'subrank'} = $rankLoc;
			$idRankHash{$id} = $rank;
			
			$lastStop = $startStopByRankHash{$rank}{$id}{'stop'};
			$lastRank = $rank;
		}
		
	} elsif($lineArray[0] eq "L") {
		
		$localrank = $lineArray[6];
		$localrank =~ s/$rankTag//;
		$id1 = $lineArray[1];
		$strand1 = $lineArray[2];
		$id2 = $lineArray[3];
		$strand2 = $lineArray[4];
		
		
		$id = "$id1\_to\_$id2";
		$xstart1 = $startStopByRankHash{$idRankHash{$id1}}{$id1}{'stop'};
		$xstop1 = $startStopByRankHash{$idRankHash{$id2}}{$id2}{'start'};
		#
		# this needs fixing, not sure how to incorporate strand
		#
		#if ($strand1 eq "-") {
		#	$xstart1 = $startStopByRankHash{$idRankHash{$id2}}{$id2}{'start'};
		#	$xstop1 = $startStopByRankHash{$idRankHash{$id1}}{$id1}{'stop'};
		#}
		$ystart1 = $startStopByRankHash{$idRankHash{$id1}}{$id1}{'subrank'};
		$ystop1 = $startStopByRankHash{$idRankHash{$id2}}{$id2}{'subrank'};
		#if ($strand1 eq "-") {
		#	$ystart1 = $startStopByRankHash{$idRankHash{$id2}}{$id2}{'subrank'};
		#	$ystop1 = $startStopByRankHash{$idRankHash{$id1}}{$id1}{'subrank'};
		#}
		
		#
		# dump line info to a file
		#
		print AOUT "$id\t$localrank\t$xstart1\t$ystart1\t$xstop1\t$ystop1\n";
	}	
}
close GFA;
close AOUT;

#
# dump node info to a file
#
foreach my $rank (sort {$a<=>$b} keys(%startStopByRankHash)) {
	foreach my $id (sort keys(%{$startStopByRankHash{$rank}})) {
		if ($rank == $maxRank) {
			$segmentScale = 0;
			if ($startStopByRankHash{$rank}{$id}{'size'} > 0){
				$segmentScale = ($startStopByRankHash{$rank}{$id}{'stop'} - $startStopByRankHash{$rank}{$id}{'start'})/$startStopByRankHash{$rank}{$id}{'size'};
			}
			$start = $stopStopByRankHash{$rank}{$id}{'start'};
			$stop = $startStopByRankHash{$rank}{$id}{'stop'};
		} else {
			
		}
		$scaledStart = $scaledStartStopHash{$startStopByRankHash{$rank}{$id}{'start'}};
		$scaledStop = $scaledStart+$startStopByRankHash{$rank}{$id}{'size'};
		$stop = $startStopByRankHash{$rank}{$id}{'start'} + $startStopByRankHash{$rank}{$id}{'size'};
		print OUT "$id\t$rank\t$startStopByRankHash{$rank}{$id}{'subrank'}\t$startStopByRankHash{$rank}{$id}{'start'}\t$stop\t$startStopByRankHash{$rank}{$id}{'size'}\n";
	}
}





close OUT;

sub resize {
    my $num = shift;
 	if ($num <= 1e1) {
 		$out = 1;
 	}
 	if ($num > 1e1 && $num <= 1e2) {
 		$out = 2;
 	}
 	if ($num > 1e2 && $num <= 1e3) {
 		$out = 3;
 	}
 	if ($num > 1e3 && $num <= 1e4) {
 		$out = 4;
 	}
 	if ($num > 1e4 && $num <= 1e5) {
 		$out = 5;
 	}
 	if ($num > 1e5 && $num <= 1e6) {
 		$out = 6;
 	}
 	if ($num > 1e6 && $num <= 1e7) {
 		$out = 7;
 	}
    return $out;
}
