#!/usr/bin/perl -w
use strict;
use FileHandle;
use List::Util qw(first);
use List::Util qw(min max);

# Morag Lewis
# morag.lewis@kcl.ac.uk

#VEP FiVar
#Perl script to filter VEP-annotated variants


##############################################################################
## Copyright 2024 Morag Lewis                                               ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##   http://www.apache.org/licenses/LICENSE-2.0                             ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
##                                                                          ##
##############################################################################


#this script takes a vep-annotated file and outputs the variants passing the input filter settings
#filter settings: minor allele frequency (0-1), impact ("high", "low", "none")
#if no MAF filter is required, use 1
#The MAF is taken from the gnomAD AF grpmax - this is the maximum allele frequency across all populations in gnomAD (genomes)
#The pathogenicity filter is whether a variant is predicted to be a high impact variant (according to CADD (>25), SpliceAI (>0.5), UTRAnnotator (scoring based on https://academic.oup.com/bioinformatics/article/37/8/1171/5905476 and https://www.nature.com/articles/s41467-019-10717-9) or AlphaMissense (>0.564), a low-impact variant (anything within a gene or which has a high ReMM score (>0.95) outside a gene), or whether no impact filter is applied

#A private allele frequency test is carried out to check the input cohort and filter out any variants present at a frequency over MAF+0.4. This can be deactivated by using the "ukbb" setting, although that also means the calls won't be output
#A biotype filter test is carried out to remove any transcripts clsssified as either a pseudogene or nonsense-mediated decay

#The script requires a header with the IDs in the right order and very specific VEP annotation (see the associated file for an example header output from the VEP, and the annotation sources)
#The output will be vcf by default; add "tab" to command line to get tabulated file
#The "ukbb" setting means the script will only output the variants, it won't output all the individual genotypes (important for very large datasets with hundreds of thousands of participants)
#The "ukbb" setting will also prevent the script from carrying out the private frequency test

#Gene-specific CADD threshold scores come from PMID: 26820543; use "set-cadd" to keep cadd threshold at 25 for all genes. 
#Note that getting gene-specific CADD scores means that only the ~2950 genes with HGMD-generated CADD scores will be included in the output
#A different set of scores can be used (in the same format: tab-separated lines of 11 columns; column 9 has MSC, column 10 has "HGMD", column 11 has the Ensembl gene ID), in which case the output will include only those genes

($#ARGV > -1) or die("Usage: process.pl combined-snvs.vcf msc-gene-score-file/set-cadd maximum-MAF impact (filter/tabulate/ukbb)\n");
my($snvs, $gcadds, $maxmaf, $impact, $setting) = @ARGV;

if ( $#ARGV < 4 ) { 
	$setting = "filt"; 
}

#set for debugging:
my $debug = 0;

my %genecadds; #store of gene-specific CADD thresholds, hashed on gene ID
if ( $gcadds ne "set-cadd" ) {
	my $gcfh = new FileHandle('<' . $gcadds);
	while(<$gcfh>) {
		my ( $stupidname, $nmut, $mean, $med, $cil, $ciu, $min, $max, $msc, $source, $geneids ) = split(/\s+/);
		if ( $source eq "HGMD" ) { #only want the HGMD ones
			chomp($geneids);
			if ( $geneids !~ /N\/A/) {
				my ( @gids ) = split( /;/, $geneids );
				for my $pip ( 0 .. $#gids ) {
					if ( exists( $genecadds{ $gids[ $pip ] } ) ) {
						#if there's a duplicate, take the higher score
						if ( $msc > $genecadds{ $gids[ $pip ] } ) { $genecadds{ $gids[ $pip ] } = $msc; }
					} else {
						$genecadds{ $gids[ $pip ]} = $msc;
					}
				}
			}
		}
	}
	$gcfh->close();
}



print STDERR "Filtering variants from $snvs\n";
#lists of consequences to pass (see https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html for description)
my $consequences = "transcript_ablation\tsplice_acceptor_variant\tsplice_donor_variant\tstop_gained\tframeshift_variant\tstop_lost\tstart_lost\ttranscript_amplification\tinframe_insertion\tinframe_deletion\tmissense_variant\tprotein_altering_variant\tsplice_donor_5th_base_variant\tsplice_region_variant\tsplice_donor_region_variant\tsplice_polypyrimidine_tract_variant\tincomplete_terminal_codon_variant\tstart_retained_variant\tstop_retained_variant\tmature_miRNA_variant";
my $geneconsq = $consequences."\tcoding_sequence_variant\tsynonymous_variant\t5_prime_UTR_variant\t3_prime_UTR_variant"; #note that intron_variant is NOT in this list
my $spliceconsq = $geneconsq."\tintron_variant"; #intron_variant is in this list because novel splice variants can be intronic

#tracker
my $chrom=0;

my $snvfh = new FileHandle('<' . $snvs);
while(<$snvfh>) {
	my(@columns) = split(/\s+/);
	if ( $columns[0] =~ /^##/) { 		#skip headers
		next;
	} elsif ( $columns[0] =~ /^#/) {		#print header with IDs; skipped unless "filter" or "tab" set
		if ( $setting =~ /tab/ ) {
			print "Chrom\tPos\tID\tRef\tAlt\tGene ID\tGene name\tTranscript\tConsequence\tMax MAF\tCADD (Phred-scaled)\tAlphaMissense\tReMM\t5UTR consequence\tAcceptor gain score\tAcceptor loss score\tDonor gain score\tDonor loss score\t";
			my $samples=join("\t",@columns[9..($#columns)]);
			print $samples."\n"; #print the header for the sample IDs               
		} elsif ( $setting =~ /filt/ ) {
			print;
		}
		next;
	} else  { 
		if ( $columns[0] ne $chrom) { #tracking
			$chrom = $columns[0];
			print STDERR "Chr $chrom\n";
		}

######################################################## Call-specific filters - not for large datasets #################################################################################

		my @callarray;
		my $calls;
		my $private=0;

		if ( $setting !~ /ukbb/ ) {
			@callarray = splice(@columns,9); #get calls
			$calls=join("\t", @callarray); #get calls as string
			if ( $calls !~ /1/ ) { next; } #skip if there isn't at least one call for this alternative allele

	                my $pip;
       		        my $total = 0;
      	       	  		my $mac = 0;
              	 		for $pip ( 0 .. $#callarray ) {
              	         		if ( $callarray[ $pip ] =~ /1\/1/ ) { #if it's a homozygous call, add 2 to the allele count
                               		$mac = $mac + 2;
                        	} elsif ( $callarray[ $pip ] =~ /\/1/ ) { #else if there's an alternative allele (*/1), add 1 to the allele count
       	                        	$mac = $mac + 1;
                        	} elsif ( $callarray[ $pip ] =~ /1\// ) { #else if there's an alternative allele (1/*), add 1 to the allele count
       	                        	$mac = $mac + 1;
               	        	}
                       		$total = $total + 2; #either way, add 2 to the total allele count
                	}
       	        	$private = $mac/$total; #calculate private allele count for later filtering once MAF is known
		}

######################################################################################################################################################################

		my $toprint = "$columns[0]\t$columns[1]\t$columns[2]\t$columns[3]\t$columns[4]\t";

		my @annots = split( /\;/, $columns[7]); #get annotations
		my @csqs = split( /\,/, $annots[$#annots]); #get consequences (the final annotation), which are repeated for each transcript	

		my $printthis = 0; #assume not printing this variant
		my $printvcf = $toprint.$columns[5]."\t".$columns[6]."\t"; #get the line ready if outputting a vcf - only want the consequences which pass the filters
		my $lastbutone = $#annots - 1;
		for my $ctr ( 0 .. $lastbutone ) {
			$printvcf = $printvcf.$annots[ $ctr ].";"; #add back the annotations before the consequences
		}
		$printvcf = $printvcf."CSQ=";
	
		for my $fred ( 0 .. $#csqs ) {
                        my $input = $csqs[ $fred ];
                        chomp($input);
                        $input = $input."|end"; #added to prevent empty fields at the end of the line causing problems in split
			my @scores = split( /\|/, $input );
			if ( $#scores < 51 ) { next; } #print STDERR "Fewer fields than expected for $toprint\n"; next; } #check the correct number of fields are generated by the split; if not, skip this variant
			my $consq = $scores[1];
			my $gene = $scores[3];
			my $ensg = $scores[4];
			my $transcript = $scores[6];
			my $biotype = $scores[7];
			my $cadd = $scores[43];
			my $fputrcons = $scores[30];
			my (@fputrannots) = split(/\&/, $scores[29]); #split the annotations by '&' for the cases where there are multiple uORFs affected
			my $fputrscore = 0; #score for 5'UTR variants, set based on the annotations (see below) - basically, 1 means interesting and 0 is not
			my $remm = $scores[48];
			my @splicearray = @scores[38,39,40,41];
			my $alphamissense = $scores[46];
			my $afmax = $scores[51]; #gnomAD AF grpmax - maximum allele frequency across all populations in gnomAD (genomes)

##################### CHECK FOR MULTIPLE MAF ENTRIES ###################################################################################################################

			if ( $afmax =~ /&/ ) { print STDERR "Multivalue grpmax for $toprint: $afmax\n"; }

################################ First filter - not interested in pseudogene or NMD transcripts ########################################################################

			if ( ( $biotype =~ /pseudogene/) || ( $biotype =~ /nonsense/ ) ) {
				if ($debug ) { print STDERR "Biotype filter failed: $biotype\n"; }
				next; 
			}

################################ Setting up MAF #################################################################################
			my $maf = 0; #assume this is an unknown variant with allele frequency=0
			if ( $afmax =~ /[0-9]/ ) { $maf = $afmax; } #if there is a max frequency value, use it
#			for my $mfc ( 0 .. $#afs ) { #go through allele frequencies
#				if ( $afs[ $mfc ] =~ /[0-9]/ ) { #if there is a numerical value
#					if ( $afs[ $mfc ] > $maf ) {
#						$maf = $afs[ $mfc ]; #if the allele frequency > stored MAF, reset stored MAF
#					}
#				} 
#			}

#			#Skip if the private allele count is more than the MAF + 0.4 (always passes for large datasets where the private allele frequency is not calculated)
			if ( $private > ( $maf + 0.4 )) { 
				if ($debug ) { print STDERR "Private allele frequency filter failed (private frequency, MAF): $private $maf\n"; }
				next; 
			}

################################ Setting up CADD, ReMM, SpliceAI, UTRAnnotator and AlphaMissense pathogenicity scores #########################################################################
			if ( $remm !~ m/[0-9]/ ) { $remm=0; } #if there isn't a remm score, assume it is 0 - no regulatory effect

			if ( $remm =~ /\&/ ) { #if there are multiple ReMM values, take the lower - don't permit false positives
				#if there are some numbers, remove the .s, they shouldn't be considered for this ("NA" is not the same as "0" for ReMM if there is another score available)
				$remm =~ s/\&\.//g; #remove .s next to &
				$remm =~ s/\.\&//g; #remove .s next to &
				my @multiremms = split( /\&/, $remm );
				$remm = min @multiremms;
			}

			my $threshold = 25; #CADD threshold if there's no defined MSC score
			if ( $gcadds ne "set-cadd" ) { #if want gene-specific CADD scores
				if ( exists( $genecadds{ $ensg } )) { #if there is an HGMD-generated CADD threshold
					$threshold = $genecadds{ $ensg }; #use it
				} else {
					$threshold = 100000; #use very high number to prevent anything passing this threshold which doesn't have an HGMD-generated CADD threshold
				}
			}

			if ( $cadd !~ /[0-9]/ ){	#if cadd score is not defined, assume it is an indel and potentially pathogenic
				$cadd = 100;
			}

			if ( $alphamissense !~ /[0-9]/ ){	#if alpha missense score is not defined, this is probably not a missense variant; set to 0
				$alphamissense = 0;
			}

			#Set up SpliceAI scores - take the maximum score
			for my $deb ( 0 .. $#splicearray ) {
				if ( $splicearray[ $deb ] =~ /^$/ ) {
					$splicearray[ $deb ] = 0; #set any nonexistent scores to 0
				}
			}
			my $spliceai = max(@splicearray);
			my $spliceaistring = join("\t", @splicearray);

			#UTRAnnotator; fputrscore previously set to 0, which means it isn't interesting. Go through the UTRAnnotator annotations to see if that should be changed
			if ( $fputrcons !~ /^$/ ) { #not interested if there's no 5'UTR annotation (!)
				my $evidence = 0; #whether this is a uORF known to be translated (assume it is not)
				my $type = "NA"; #type of uORF (only matters for start gain)
				my $kozak = "NA"; #Kozak consensus sequence strength
				my $altstop = 1; #whether there is an alternative stop codon or not (assume there is, because that means the stop-loss doesn't matter - ie default to not pathogenic; only matters for stop loss)
				#go through annotations to set these; UTRAnnotator uses '&' where there are multiple uORFs disrupted by a variant so split by that first
				for my $fpa ( 0 .. $#fputrannots) {
					my (@theseannots) = split(/\:/, $fputrannots[ $fpa ]); #get the individual annotations
					for my $fpv ( 0 ..$#theseannots) {
						my ( $ann, $fpval ) = split(/\=/, $theseannots[ $fpv ]);
						if ( $ann eq "KozakStrength" ) {
							$kozak = $fpval;
						} elsif ( $ann eq "AltStop" ) {
							if ( $fpval eq "False" ) { $altstop = 0; }
						} elsif ($ann eq "Evidence" ) {
							if ( $fpval eq "True" ) { $evidence = 1; }
						} elsif ( $ann eq "type" ) {
							$type = $fpval;
						}
					}
					#Scoring based on https://academic.oup.com/bioinformatics/article/37/8/1171/5905476 and https://www.nature.com/articles/s41467-019-10717-9
					if ( $fputrcons =~ /5_prime_UTR_uORF_frameshift_variant/ ) {
						if ( $evidence == 1 ) { $fputrscore = 1; }
					} elsif ( $fputrcons =~ /5_prime_UTR_stop_codon_gain_variant/ ) {
						if ( $evidence == 1 ) { $fputrscore = 1; }
					} elsif ( $fputrcons =~ /5_prime_UTR_premature_start_codon_loss_variant/ ) {
						if ( $evidence == 1 ) { $fputrscore = 1; }
					} elsif ( $fputrcons =~ /5_prime_UTR_stop_codon_loss_variant/ ) {
						if ( ( $altstop == 0 ) && ( ( $evidence == 1 ) || ( $kozak =~ /Moderate|Strong/) ) ) { $fputrscore = 1; }
					} elsif ( $fputrcons =~ /5_prime_UTR_premature_start_codon_gain_variant/ ) {
						if ( ( $kozak  =~ /Moderate|Strong/ ) && ( $type =~ /oORF/ ) ) { $fputrscore = 1; }
					}
					
					#if there's a positive fputrscore at this point, exit the loop - only one likely-disrupted uORF is required to make this an interesting variant
					if ( $fputrscore == 1 ) { last; }
				}
			} else {
				$fputrcons = "."; #set for output
			}
				
################################ Filter step #################################################################################
			#have to cycle through consequences if there are lots of them
			my @indivcsqs = split(/\&/, $consq );
			my $jill;
			for $jill ( 0 .. $#indivcsqs ) {
				my $includethis = 0; #assume not keeping this consequence

				#MAF filter
				if ( $maf < $maxmaf ) {
					if ( $impact eq "high" ) { # either severe consequences and high pathogenicity, or missense variant and AlphaMissense likely pathogenic and cadd > 20, or 5'UTR variants and UTRAnnotator score = 1, or splice-affecting variants and spliceAI score > 0.5
						if ( ( ( $consequences =~ /$indivcsqs[ $jill ]/) && ( $cadd > $threshold ) ) || ( ( $indivcsqs[ $jill ] =~ /missense_variant/ ) && ( $alphamissense > 0.564) && ( $cadd > 20) )  || ( ($indivcsqs[ $jill ] =~ /5_prime_UTR_variant/) && ( $fputrscore == 1 ) ) || ( ( $spliceconsq =~ /$indivcsqs[ $jill ]/ ) && ( $spliceai > 0.5 ) ) ) {
							$includethis = 1;
						} else {
							if ($debug ) { print STDERR "High impact filter failed (consequence, CADD, AlphaMissense, UTRAnnotator, SpliceAI): $indivcsqs[ $jill ] $cadd $alphamissense $fputrcons $spliceai\n"; }
						}
					} elsif ( $impact eq "low" ) { #either consequence is within a gene or variant has a high ReMM score (<0.95) or variant is intron variant with high spliceai score
						if ( ( $geneconsq =~ /$indivcsqs[ $jill ]/) || ( $remm > 0.95 ) || ( ($indivcsqs[ $jill ] =~ /intron_variant/ ) && ( $spliceai > 0.5 ) ) ) {
							$includethis = 1;
						} else {
							if ($debug ) { print STDERR "Low impact filter failed (consequence, ReMM, SpliceAI): $indivcsqs[ $jill ] $remm $spliceai\n"; }
						}
					} elsif ( $impact eq "none" ) {
						$includethis = 1;
					} else { print STDERR "Input impact not recognised: $impact\n"; } 
				} else {
					if ($debug ) { print STDERR "MAF filter failed ( $maxmaf ): $maf\n"; }
				}

################################## Printing (for tab-separated output, which is one line per consequence) #####################
				if ( $includethis == 1 ) {
					if ( $setting =~ /tab/ ) {
						if ( $setting =~ /ukbb/ ) { #for outputting biobank scores in tabulated format - typically only used for testing; just doesn't include calls
							print $toprint.$ensg."\t".$gene."\t".$transcript."\t".$consq."\t".$maf."\t".$cadd."\t".$alphamissense."\t".$remm."\t".$fputrcons."\t".$spliceaistring."\n";
						} else {
							print $toprint.$ensg."\t".$gene."\t".$transcript."\t".$consq."\t".$maf."\t".$cadd."\t".$alphamissense."\t".$remm."\t".$fputrcons."\t".$spliceaistring."\t".$calls."\n";
						}
					} else {
						$printthis=1;
						my $consline = $csqs[ $fred ];
						$consline =~ s/^CSQ\=//; #get rid of the initial CSQ=, if it's there
						$printvcf = $printvcf.$consline.",";
					}
				}
			}
		}

################################## Printing (for vcf output, which is one line per variant) #################################

		if ( ( $setting !~ /tab/ ) && ( $printthis != 0 )) { 
			$printvcf =~ s/\,$//; #remove trailing comma
			if ( $setting =~ /ukbb/ ) {
				print $printvcf."\t".$columns[8]."\n"; #print without calls for ukbb setting
			} else {
				print $printvcf."\t".$columns[8]."\t".$calls."\n"; #print with calls
			}
		}
	}
}
$snvfh->close();

