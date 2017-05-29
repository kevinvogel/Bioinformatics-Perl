#!/usr/bin/perl -w

use strict;

# Program description: Perl program designed to find all the open reading frames (ORFs) in the input fasta sequence.

# Add a “Welcome” message to the beginning of the program that describes to the user what the program does.
print "Welcome! This program is designed to find all the open reading frames (ORFs) in a input fasta DNA sequence file.\n\n";

print "Please enter the file name in fasta format (e.g. <test.fasta>):\n";
chomp(my $fileName = <STDIN>);

# Check if the file has .fasta extention
isAFasta($fileName);

# Capture sequence in a hash 
my %forward_sequence = readFile($fileName);

# Check to see if the fasta file contains a DNA sequence
isDNA(\%forward_sequence);

# Prompts user for ORF nucleotide length
my $ORF_LENGTH = set_ORF_parameters();

# reverse complement the DNA.Key is the id and value is the reversed sequence
my %reversed_sequence = reverseDNA(\%forward_sequence);

# Key has the DNA ID, value has the array which contains all 6 frames.
my %framed_Sequence = frame(\%forward_sequence, \%reversed_sequence);

# printing the data structure
orf(\%framed_Sequence, $ORF_LENGTH);

sub readFile {
	# Author: Eunice Onyekweli 
	# function: to capture the sequence information in a hash
	# input: .fasta file name
	# return: hash (keys are name of the sequence,values are DNA sequence)
	# variable for the file name
	my $file = shift;
	chomp($file);
	# hash to contain sequence information from the file. keys are the headers, values are the sequences
	my %sequence = ();
	open(SEQ, $file) or die "Could not open file '$file' $!";
	my $header;
	while(my $line = <SEQ>) {
	chomp $line;
		if($line =~ /^\s*$/) {
		} elsif($line=~ /^>/) {
			$header = $line;
			if(!exists $sequence{$header}) {
				$sequence{$header} = ();
			}
		} else {
			if ( defined $sequence{$header}) {
				$sequence{$header} .= $line;		
			} else {
				$sequence{$header} = $line;
			}					
		}
	}
	return %sequence;
}

sub isAFasta {
	# Author: Eunice Onyekweli 
	# function: to check if the correct file extention is given and the file exists in the directory
	# input: fasta file name
	# return: n/a. if wrong file extention is given the sub will prompt user to enter correct file name and die.
		#variable contains the user given file name
		my $subextention = shift;
		chomp $subextention;
		# check if the file is in fasta format
		if($subextention =~ /.fasta$/) {
			print "The input file is in .fasta \n\n";
			# check if file exists in the directory
			if (-e $subextention) {
				print "File found in the directory!\n\n";
			}
		} else {
			die "User did not enter .fasta file. Please use protein file in fasta format.\n";
		}	
}

sub isDNA {
	# Author: Kevin Vogel
	# function: subroutine allows the user to know if the fasta file contains DNA sequences
	# input: hash reference 
	# return: print statement
	my %forward_DNA = %{$_[0]};
	my @dna_nucleotides = ("A", "G", "T", "C");
	foreach my $key (keys %forward_DNA) {
	    my @temp_array = split("", uc $forward_DNA{$key});
	    foreach my $letter (@temp_array) {
		if(grep({$_ eq $letter} @dna_nucleotides)) {
		    print "$key file contains DNA!\n\n";
		    last;
		} else {
		    print "$key file does not contain DNA.\n\n";
		    last;
		}
	    }
	}
}
	
sub set_ORF_parameters {
	# Author: Kevin Vogel
	# function: prompts users for the minimum number of ORF nucleotides to search for (minimum of 50 nt).
	# input: integer >= 50 nt bases
	# return: integer
	print "Please enter the minimum ORF length (minimum of 50):\n";
	chomp(my $input_orfs_len = <STDIN>);
	if ($input_orfs_len < 50) {
	    die "Error. Must search for a minimum of 50 ORF nucleotides.\n";
	    exit;
	} else {
	    return $input_orfs_len;
	}
}

sub reverseDNA {
	# Author: Kevin Vogel
	# function: To reverse the DNA
	# input: hash reference
	# return: hash
	my (%forwardDNA) = %{$_[0]};
	# variable to capture reversed sequence.keys will contain sequence id,values will contain reverse complement sequence of values from %forwardDNA
	my %reversedDNA;
	# looping through the original DNA to reverse it
	foreach my $forwardkey(keys %forwardDNA) {
		chomp($forwardkey);
		chomp($forwardDNA{$forwardkey});
		#Capturing reverse DNA in $temp
		my $temp = reverse($forwardDNA{$forwardkey});
		#Capturing reversed complemented DNA in $rDNA
		my $rDNA = reverseComplement($temp);
		chomp($rDNA);
		$reversedDNA{$forwardkey} = $rDNA;
	}
	# hash contains ids as key and reverse complemented DNA as the value
	return %reversedDNA;
}

sub reverseComplement {
	# Author: Fatma Onmus-Leone
	# function: Creates the complement of a given sequence
	# input: DNA string
	# return: complement of a DNA string
	my ($dnaToReverse) = shift;
	chomp($dnaToReverse);
	my $reverseComp = $dnaToReverse;
	$reverseComp =~ tr/ACGT/TGCA/;
	$reverseComp =~ s/\W*//;
	return $reverseComp;
}

sub frame {
	# Author: Fatma Onmus-Leone
	# function: to find the sequence for all 6 frames
	# input: forward and reverse DNA hash
	# return: return hash of DNA with all 6 frames
	my %forward_DNA = %{$_[0]};
	my %reverse_DNA = %{$_[1]};
	my (@allFrames, @reversed_frames, @forward_frames);
	# hash to capture all 6 frames as a value in an array.$framedDNA{sequenceID}=[forward_frame1,forward_frame2,forward_frame3,reverse_sequence1,reverse_sequence2,reverse_sequence3]
	my (%framedDNA);
		# capturing forward frames 
		foreach my $f(keys %forward_DNA) {
			chomp($f); 
			chomp($forward_DNA{$f});
			my $forwardframe_1 = $forward_DNA{$f};
			my $forwardframe_2 = substr($forwardframe_1, 1, length($forwardframe_1));
			my $forwardframe_3 = substr($forwardframe_1, 2, length($forwardframe_1));
			# forward frame 1
			$framedDNA{$f}[0] = $forwardframe_1;
			# forward frame 2
			$framedDNA{$f}[1] = $forwardframe_2;
			# forward frame 3
			$framedDNA{$f}[2] = $forwardframe_3;
		}
		# capturing reverse frames
		foreach my $r(keys %reverse_DNA) {
			chomp($r); 
			chomp($reverse_DNA{$r});
			my $reverseframe_1 = $reverse_DNA{$r};
			my $reverseframe_2 = substr($reverseframe_1, 1, length($reverseframe_1));
			chomp ($reverseframe_2);
			my $reverseframe_3 = substr($reverseframe_1, 2, length($reverseframe_1));
			chomp ($reverseframe_3);
			# reverse frame 1
			$framedDNA{$r}[3] = $reverseframe_1;
			# reverse frame 2
			$framedDNA{$r}[4] = $reverseframe_2;
			# reverse frame 3
			$framedDNA{$r}[5] = $reverseframe_3;
		}
	
	#Key has the DNA ID, value has the array which contains all 6 frames.e.g $reversed_sequence{sequenceID}=[forward_frame1,forward_frame2,forward_frame3,reverse_sequence1,reverse_sequence2,reverse_sequence3]	
	return %framedDNA;
}

sub codon {
	# Author: Fatma Onmus-Leone
	# input: A string of a sequence
	# function: splits the given sequence into codons.
	# return: an array of a sequence where each element of an array is a codon
	my $d = shift;
	chomp($d);
	my $position = 0;
	my @CODON_DNA;
	my $codonedDNA;
		for (my $count = 0; $count < length($d);) {
			my $codon = substr($d,$position,3);
			push(@CODON_DNA, $codon);
			if(defined $codonedDNA) {
				$codonedDNA = $codonedDNA." ";
				$codonedDNA .= $codon;
			} else {
				$codonedDNA=$codon;
			}
			$count = $count + 3;
			$position = $position + 3;
		}
	return @CODON_DNA;
}

sub orf {
	# Author: Fatma Onmus-Leone
	# input: Takes two variables. A Hash containing sequence id as the key and 6 frames where each frame is an element of the array.A string of orf length.
	# function: Finds the open reading frame and prints out information such as Sequence ID,Frame number,Length of the orf and sequence of the orfs.
	# returns: Doesn't return anything to main. However sends sequence information to Sub codon as a string and takes an array of codons in return.
	my %frames = %{$_[0]};
	my $Size = $_[1];
	my @Codons = "";
	my $frame = 0;
	# sequence length
	my $SSize;
	my $FrameNumber;
	my $sequence;
	# start locations captured here
	my @FramesPos = ();
	my @Orfs = ();
	my $Orf = "";
	my $OrfBegin = 0;
	# individual codon
	my $Codon = "";
	# Codon number within the @Codon
	my $CodonNumber;
	my $j = 0;
	my $Length = 0;
	my $Offset = 0;
	# looping through the hash 
	foreach my $Header(keys %frames) {
	chomp($Header);
		# looping through the frames
		for (my $index = 0; $index < 6; $index++) {
		@Orfs = "";
		$sequence = $frames{$Header}[$index];
		$FrameNumber = $index + 1;
		# retrieved the frame in an array(one codon per array position)
		my @formatDNA = codon($frames{$Header}[$index]);
		$SSize = length($frames{$Header}[$index]);
		my $CodonNumber = 1;
		my $baseNum = 1;
		@FramesPos = "";
		foreach $Codon(@formatDNA) {
		# if codon has start the OrfBegin is set to 1 otherwise it is 0. When OrfBegin is 1 keep concatinating until stop codon is found.
			if($Codon eq "ATG") {
				if($OrfBegin == 1) {
					$Orf = $Orf.$Codon;
				} else {
					push(@FramesPos, $baseNum);
					$Orf = $Codon;
					$OrfBegin = 1;
				}
			# if OrfBegin is 1, which means start found, look for stop. When stop is found set orfbegin to 0 so new search can start. Add the orf to @orfs.
			} elsif($Codon eq "TAA" or $Codon eq "TAG" or $Codon eq "TGA") {
				if($OrfBegin == 1){
					$Orf = $Orf.$Codon;
					push(@Orfs, $Orf);
					$OrfBegin = 0;
				}
			} else {
				if($OrfBegin == 1) {
					$Orf = $Orf.$Codon;
				}
			}
		$CodonNumber++;
		$baseNum = $baseNum + 3;
		}
		my $count = 1;
		# looping through found orf to print out.
		if(scalar @Orfs > 0){
			for($j = 0; $j < scalar @Orfs; $j++) {
				$Length = length($Orfs[$j]);
				if($Length >= $Size) {
					if($FrameNumber == 1) {	
						print "$Header | FRAME = $FrameNumber POS = $FramesPos[$j] LEN = $Length\n" ;
					} elsif($FrameNumber == 2){
						print "$Header | FRAME = $FrameNumber POS = ",($FramesPos[$j]+1)," LEN = $Length\n" ;
					} elsif($FrameNumber == 3){
						print "$Header | FRAME = $FrameNumber POS = ",($FramesPos[$j]+2)," LEN = $Length\n" ;
					} elsif($FrameNumber == 4){
						print "$Header | FRAME = $FrameNumber POS = ",($FramesPos[$j]*-1)," LEN = $Length\n" ;
					} elsif($FrameNumber == 5){
						print "$Header | FRAME = $FrameNumber POS = ",(($FramesPos[$j]+1)*-1)," LEN = $Length\n" ;
					} elsif($FrameNumber == 6){
						print "$Header | FRAME = $FrameNumber POS = ",(($FramesPos[$j]+2)*-1)," LEN = $Length\n" ;
					} else {
					}
					$count = 1;	
					for ($Offset = 0; $Offset < length($Orfs[$j]) ; $Offset += 3) {
						if($count<15) {
							print substr($Orfs[$j], $Offset, 3)." ";
						} else {
							print substr($Orfs[$j], $Offset, 3)."  \n";
							# count for printing 15 codon in a line
							$count = 0;
						}
					$count++;
					}
					print "\n\n";
					}
				}
			}	
		}	
	}
}
