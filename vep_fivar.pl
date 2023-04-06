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
## Copyright 2020 Morag Lewis                                               ##
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


#this script takes a vep-annotated file and outputs the variants passing whichever filter is chosen
#filter settings:
#0 - no filter, just output (default)
#1 - MAF < 10%, no pathogenicity filter
#2 - MAF < 0.5%, no pathogenicity filter)
#3 - MAF < 10%, pathogenicity filter
#4 - MAF < 0.5%, pathogenicity filter
#5 - no MAF filter, pathogenicity filter

#The script requires a header with the IDs in the right order and very specific VEP annotation (see the associated file for an example header output from the VEP, and the annotation sources)
#The output will be vcf by default; add "tab" to command line to get tabulated file
#The "ukbb" setting means the script will only output the variants, it won't output all the individual genotypes (important for very large datasets with hundreds of thousands of participants)
#The "ukbb" setting will also prevent the script from carrying out the private frequency test
#The "ukbb" setting will use the maximum MAF from all populations and all databases, otherwise the script will assume the Non-Finnish European setting and will check GnomAD first, then the 1000 Genomes project, then TopMed, then ESP6500
#Gene-specific CADD threshold scores come from PMID: 26820543; use "set-cadd" to keep cadd threshold at 25 for all genes. 
#Note that getting gene-specific CADD scores means that only the ~2950 genes with HGMD-generated CADD scores will be included in the output
#Pathogenicity prediction filter settings are hardwired: cadd > 25, |sutr| > 1, spliceai > 0.5


($#ARGV > -1) or die("Usage: process.pl combined-snvs.vcf msc-gene-score-file/set-cadd filternumber (filter/tabulate/ukbb)\n");
my($snvs, $gcadds, $fnum, $setting) = @ARGV;

if ( $#ARGV < 3 ) { 
	$setting = "filt"; 
}

if ( $fnum > 5 ) { #set the $filter number (assume 0 if the number is over 5, and if it is not a number this will error and exit)
	$fnum = 0;
}

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
my $consequences = "transcript_ablation\tsplice_acceptor_variant\tsplice_donor_variant\tstop_gained\tframeshift_variant\tstop_lost\tstart_lost\ttranscript_amplification\tinframe_insertion\tinframe_deletion\tmissense_variant\tprotein_altering_variant\tsplice_region_variant\tincomplete_terminal_codon_variant\tstart_retained_variant\tstop_retained_variant\tmature_miRNA_variant";
my $geneconsq = $consequences."\tcoding_sequence_variant\tsynonymous_variant\t5_prime_UTR_variant\t3_prime_UTR_variant"; #note that intron_variant is NOT in this list
my $spliceconsq = $geneconsq."\tintron_variant"; #intron_variant is in this list because splice variants can be intronic

#MAF filters
my $commonmaf=0.1;
my $raremaf=0.005;

#tracker
my $chrom=0;

my $snvfh = new FileHandle('<' . $snvs);
while(<$snvfh>) {
	my(@columns) = split(/\s+/);
	if ( $columns[0] =~ /^##/) { 		#skip headers
		next;
	} elsif ( $columns[0] =~ /^#/) {		#print header with IDs; skipped unless "filter" or "tab" set
		if ( $setting =~ /tab/ ) {
			print "Chrom\tPos\tID\tRef\tAlt\tGene ID\tGene name\tTranscript\tConsequence\tgnomAD NFE\t1000G EUR\tTOPMED\tESP6500\tCADD (Phred-scaled)\tREVEL\tMPC\tExAC pLI\tReMM\tSutr score\tAUG_event\tAUGgain in frame\tAUGgain with stop\tAUGloss in frame\tAUGloss with stop\tSTOP_event\tSTOPgain with start\tSTOPloss with start\tAcceptor gain score\tAcceptor gain position\tAcceptor loss score\tAcceptor loss position\tDonor gain score\tDonor gain position\tDonor loss score\tDonor loss position\t";
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

		my @callarray;
		my $calls;

		if ( $setting !~ /ukbb/ ) {
			@callarray = splice(@columns,9); #get calls
			$calls=join("\t", @callarray); #get calls as string
			if ( $calls !~ /1/ ) { next; } #skip if there isn't at least one call for this alternative allele
		}

		my $toprint = "$columns[0]\t$columns[1]\t$columns[2]\t$columns[3]\t$columns[4]\t";

		my @annots = split( /\;/, $columns[7]); #get annotations
		my @csqs = split( /\,/, $annots[$#annots]); #get consequences, which are repeated for each transcript	

		my $printthis = 0; #assume not printing this variant
		my $printvcf = $toprint.$columns[5]."\t".$columns[6]."\t"; #get the line ready if outputting a vcf - only want the consequences which pass the filters
		my $lastbutone = $#annots - 1;
		for my $ctr ( 0 .. $lastbutone ) {
			$printvcf = $printvcf.$annots[ $ctr ].";"; #add back the annotations before the consequences
		}
		$printvcf = $printvcf."CSQ=";
	

		my $fred=0;
		while ( $fred <= $#csqs ) { #Now that consequences are being handled separately, I think it's best not to stop early
                        my $input = $csqs[ $fred ];
                        chomp($input);
                        $input = $input."|end"; #added to prevent empty fields at the end of the line causing problems in split
			my @scores = split( /\|/, $input );
			if ( $#scores < 79 ) { print STDERR "Fewer fields than expected for $toprint\n"; } #check the correct number of fields are generated by the split
			my $consq = $scores[1];
			my $gene = $scores[3];
			my $ensg = $scores[4];
			my $transcript = $scores[6];
			my $biotype = $scores[7];
			my $cadd = $scores[28];
			my $exacpli = $scores[30];
			my $revel = $scores[31];
			my $sutr = $scores[42];
			my $remm = $scores[54];
			my $spliceaisnventry = $scores[56];
			my $spliceaiindentry = $scores[58];
			my $topmed = $scores[60];
			my $G1000_eur = $scores[63];
			my $mpc = $scores[68];
			my $gnomad_nfe = $scores[75];
			my $esp6500 = $scores[78];
			my $sutrstring;
			if ( (defined( $sutr )) && ( $sutr ne "" ) ) {
				#NB at present this doesn't check for the rest of these being "", only the score itself. This can result in spaces in the output file if "tab" is selected, although everything stays correctly tabulated.
				$sutrstring = $scores[42]."\t".$scores[33]."\t".$scores[51]."\t".$scores[34]."\t".$scores[48]."\t".$scores[35]."\t".$scores[39]."\t".$scores[40]."\t".$scores[41]."\t";
			} else {
				$sutrstring = ".\t.\t.\t.\t.\t.\t.\t.\t.\t";
			}

################################ Setting up MAF #################################################################################
			my $maf=".";
			if ( $setting !~ /ukbb/ ) {
				#use gnomAD NFE MAF if possible, else 1000 Genomes MAF, else TOPMED, else ESP6500, else assume it's an unknown allele and leave maf as "."
				if ((defined($gnomad_nfe)) && ( $gnomad_nfe =~ /[0-9]/) ) {
					$maf=$gnomad_nfe;
				} elsif ((defined($G1000_eur)) && ( $G1000_eur =~ /[0-9]/) ) {
					$maf=$G1000_eur;
				} elsif ((defined($topmed)) && ( $topmed =~ /[0-9]/) ) {
					$maf=$topmed;
				} elsif ((defined($esp6500)) && ( $esp6500 =~ /[0-9]/) ) {
					$maf=$esp6500;
				}

				#set undefined AFs to "."
				if ( ( !(defined($gnomad_nfe))) || ( $gnomad_nfe !~ m/[0-9]/) ) { $gnomad_nfe="."; }
				if ( ( !(defined($G1000_eur)) ) || ( $G1000_eur !~ m/[0-9]/) ) { $G1000_eur="."; }
				if ( ( !(defined($topmed)) ) || ( $topmed !~ m/[0-9]/) ) { $topmed="."; }
				if ( ( !(defined($esp6500))) || ( $esp6500 !~ m/[0-9]/) ) { $esp6500="."; }		

				$maf =~ s/NA/\./g; #substitute "." for NA
                                
				if ( $maf =~ /\&/ ) { #if there are multiple MAF values, take the higher - don't permit false positives
					if ( $maf !~ /[0-9]/ ) { #catch the case where it's just a string of "."
						print STDERR "Warning; no numerical MAF values in $maf\n";
						$maf = ".";
					} else { #if there are some numbers, allow the . to become a 0, it won't be kept
						$maf =~ s/\&\./\&0/g; #swap out any NAs for "0"
						$maf =~ s/\.\&/0\&/g; #swap out any NAs for "0"
						my @multimafs = split( /\&/, $maf );
						$maf = max @multimafs;
					}
				}

				#final check:
				if ( $maf !~ m/[0-9]/ ) { $maf="."; }
			} else { #for UKBB analysis just want the maximum MAF
				my $mafstring = $scores[60]."_".$scores[62]."_".$scores[63]."_".$scores[64]."_".$scores[65]."_".$scores[66]."_".$scores[70]."_".$scores[71]."_".$scores[72]."_".$scores[73]."_".$scores[74]."_".$scores[75]."_".$scores[76]."_".$scores[78]."_".$scores[79];
				$mafstring =~ s/\&/_/g; #substitute "_" for &
				$mafstring =~ s/NA/\./g; #substitute "." for NA
				my @multimafs = split( /_/, $mafstring );
				my $tmpmaf=0;
				my $bob;
				for $bob ( 0 .. $#multimafs ) {
					if ( $multimafs[ $bob ] =~/[0-9]/ ) { #if there's a numerical entry
						if ( $multimafs[ $bob ] > $tmpmaf ) { $tmpmaf = $multimafs[ $bob]; } #if it's greater than tmpmaf, reset tmpmaf
					}
				}
				$maf = $tmpmaf;
			}
			#finally, if MAF is ., assume it's not a known variant and set it to 0 for filtering purposes
			if ( $maf eq "." ) { $maf=0; }


######################################################## Calculating private frequency #################################################################################
			if ( $setting !~ /ukbb/ ) {
		                my $pip;
        		        my $total = 0;
       	       	  		my $mac = 0;
               	 		for $pip ( 0 .. $#callarray ) {
               	         		if ( $callarray[ $pip ] =~ /1\/1/ ) { #if it's a homozygous call, add 2 to the allele count
                                		$mac = $mac + 2;
	                        	} elsif ( $callarray[ $pip ] =~ /1/ ) { #else if there's an alternative allele, add 1 to the allele count
        	                        	$mac = $mac + 1;
                	        	}
                        		$total = $total + 2; #either way, add 2 to the total allele count
	                	}
        	        	my $private = $mac/$total; #calculate private allele count
				#Skip if the private allele count is more than the MAF + 0.4
               			if ( $private > ($maf + 0.4 )) { 
					$fred=$fred +1;
					next; 
				}
			}

################################ Setting up CADD, ReMM, Sutr and pathogenicity scores #########################################################################
			if ( !(defined($cadd)) ) { $cadd="."; } elsif ( $cadd !~ m/[0-9]/ ) { $cadd="."; }
			if ( !(defined($revel)) ) { $revel="."; } elsif ( $revel !~ m/[0-9]/ ) { $revel="."; } #not currently using REVEL
			if ( !(defined($remm)) ) { $remm=0; } elsif ( $remm !~ m/[0-9]/ ) { $remm=0; } #if there isn't a remm score, assume it is 0 - no regulatory effect
			if ( !(defined($sutr)) ) { $sutr=0; } elsif ( $sutr !~ m/[0-9]/ ) { $sutr=0; } #if there isn't a sutr score, assume it is 0 - no difference in transcript abundance

			if ( $sutr =~ /\&/ ) { #if there are multiple Sutr values, take the lower - don't permit false positives (and remember Sutr can be negative)
				#if there are some numbers, remove the .s, they shouldn't be considered for this ("NA" is not the same as "0" for Sutr if there is another score available)
				$sutr =~ s/\&\.//g; #remove .s next to &
				$sutr =~ s/\.\&//g; #remove .s next to &
				$sutr =~ s/-//g; #convert negative numbers to positives (effectively, take absolute values)
				my @multisutrs = split( /\&/, $sutr );
				$sutr = min @multisutrs;
			}
			$sutr = abs($sutr); #take absolute value of sutr
			#Sutr filtering done in final step since it only applies to 5'UTR variants			

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

			my $path = 100; #Currently assuming anything without a CADD score is an indel and potentially pathogenic
				
			if ( $cadd !~ /^\.$/ ){		#use the cadd score, if defined
				$path = $cadd;
			}

			#Set up SpliceAI scores
			my @spliceaisnvarray = split(/\&/, $spliceaisnventry );
			my @spliceaiindarray = split(/\&/, $spliceaiindentry );
			my $spliceaistringsnv = "";
			my $spliceaistringind = "";
			my $spliceaistring = ".\t.\t.\t.\t.\t.\t.\t.\t";
			my $spliceaisnv = 0;
			my $spliceaiind = 0;
			my $spliceai = 0;

			#go through the splice ai fields. If there's a score defined, set spliceai score to whichever is the highest, and add the details to the splice ai string for later printing (tabulation only)
			my $deb;
			for $deb ( 2 .. 5 ) {
				my $har = $deb + 4;
				if ( ( defined( $spliceaisnvarray[ $deb ] ) ) && ( $spliceaisnvarray[ $deb ] =~/[0-9]/) ) {
					if ( $spliceaisnvarray[ $deb ] > $spliceaisnv ) {
						$spliceaisnv = $spliceaisnvarray[ $deb ];
					}
					$spliceaistringsnv = $spliceaistringsnv.$spliceaisnvarray[ $deb ]."\t".$spliceaisnvarray[ $har ]."\t";
				} else {
					$spliceaistringsnv = $spliceaistringsnv.".\t.\t";
				}

				if ( ( defined( $spliceaiindarray[ $deb ] ) ) && ( $spliceaiindarray[ $deb ] =~/[0-9]/) ) {
					if ( $spliceaiindarray[ $deb ] > $spliceaiind ) {
						$spliceaiind = $spliceaiindarray[ $deb ];
					}
					$spliceaistringind = $spliceaistringind.$spliceaiindarray[ $deb ]."\t".$spliceaiindarray[ $har ]."\t";
				} else {
					$spliceaistringind = $spliceaistringind.".\t.\t";
				}
			} 

			#check there isn't data in both the SNV and Indel splice AI files for this variant 
			if ( ( $spliceaisnv > 0 ) && ( $spliceaiind > 0 ) ) {
				print STDERR "Warning! Data found in both SNV and Indel splice AI files for $toprint\n";
			} 

			#set the final score and string if either string is different from "........" (which is what spliceaistring is to start with)
			if ( $spliceaistringsnv ne $spliceaistring ) {
				$spliceaistring = $spliceaistringsnv; 
				$spliceai = $spliceaisnv; 
			} elsif ( $spliceaistringind ne $spliceaistring ) {
				$spliceaistring = $spliceaistringind;
				$spliceai = $spliceaiind; 
			}
			

################################ Dealing with MPC and ExacPLI (not a filter step) #################################################################################

			if ( !(defined($mpc)) ) { $mpc="."; } elsif ( $mpc !~ m/[0-9]/ ) { $mpc="."; }
			if ( !(defined($exacpli)) ) { $exacpli="."; } elsif ( $exacpli !~ m/[0-9]/ ) { $exacpli="."; }

			$mpc =~ s/NA/\./g; #substitute "." for NA

			if ( $mpc =~ /\&/ ) { #if there are multiple MPC values, take the higher - allow false positives
				if ( $mpc !~ /[0-9]/ ) { #catch the case where it's just a string of "." - don't want to return a 0
					print STDERR "Warning; no numerical MPC values in $mpc\n";
					$mpc = ".";
				} else { #if there are some numbers, allow the . to become a 0, it won't be kept
					$mpc =~ s/\&\./\&0/g; #swap out any NAs for "0"
					$mpc =~ s/\.\&/0\&/g; #swap out any NAs for "0"
					my @multimpcs = split( /\&/, $mpc );
					$mpc = max @multimpcs;
				}
			}

################################ Filter step #################################################################################
			#have to cycle through consequences if there are lots of them
			my @indivcsqs = split(/\&/, $consq );
			my $jill;
			for $jill ( 0 .. $#indivcsqs ) {
				my $includethis = 0; #assume not keeping this consequence

				#filter
				if ( $fnum == 0 ) {
					$includethis = 1;
				} elsif ( ( $fnum == 1 ) && (( ($geneconsq =~ /$indivcsqs[ $jill ]/) || ( $remm > 0.95 ) || ( ($indivcsqs[ $jill ] =~ /intron_variant/ ) && ( $spliceai > 0.5 ) ) ) && ( $maf < $commonmaf )) ) { 
					#1. common variant filter - common AF and either consequence within a gene or high ReMM score (> 0.95)
					$includethis = 1;
				} elsif ( ( $fnum == 2 ) && (( ($geneconsq =~ /$indivcsqs[ $jill ]/) || ( $remm > 0.95 ) || ( ($indivcsqs[ $jill ] =~ /intron_variant/ ) && ( $spliceai > 0.5 ) ) ) && ( $maf < $raremaf )) ) {
					#2. rare MAF and either consequence within a gene or high ReMM score (> 0.95)
					$includethis = 1;
				} elsif ( ( $fnum == 3 ) && ( ( $maf < $commonmaf ) && ( ( ($consequences =~ /$indivcsqs[ $jill ]/) && ( $path > $threshold ) ) || ( ($indivcsqs[ $jill ] =~ /5_prime_UTR_variant/) && ( $sutr > 1 ) )  || ( ( $spliceconsq =~ /$indivcsqs[ $jill ]/ ) && ( $spliceai > 0.5 ) )   ) )) { 
					#3. common AF and either severe consequences and high pathogenicity, or 5'UTR variants and Sutr > 1, or splice-affecting variants and spliceAI score > 0.5
					$includethis = 1;
				} elsif ( ( $fnum == 4) && ( ( $maf < $raremaf ) && ( ( ($consequences =~ /$indivcsqs[ $jill ]/) && ( $path > $threshold ) ) || ( ($indivcsqs[ $jill ] =~ /5_prime_UTR_variant/) && ( $sutr > 1 ) )  || ( ( $spliceconsq =~ /$indivcsqs[ $jill ]/ ) && ( $spliceai > 0.5 ) )   ) )) {
					#4. rare variant filter - rare AF and either severe consequences and high pathogenicity, or 5'UTR variants and Sutr > 1, or splice-affecting variants and spliceAI score > 0.5
					$includethis = 1;
				} elsif ( ( $fnum == 5 ) && (  ( ( ($consequences =~ /$indivcsqs[ $jill ]/) && ( $path > $threshold ) ) || ( ($indivcsqs[ $jill ] =~ /5_prime_UTR_variant/) && ( $sutr > 1 ) )  || ( ( $spliceconsq =~ /$indivcsqs[ $jill ]/ ) && ( $spliceai > 0.5 ) )   ) )) {
					#5. no MAF filter at all, and either severe consequences and high pathogenicity, or 5'UTR variants and Sutr > 1, or splice-affecting variants and spliceAI score > 0.5
					$includethis = 1;
				}

				if ( $includethis == 1 ) {
					if ( $setting =~ /tab/ ) {
						print $toprint.$ensg."\t".$gene."\t".$transcript."\t".$consq."\t".$gnomad_nfe."\t".$G1000_eur."\t".$topmed."\t".$esp6500."\t".$cadd."\t".$revel."\t".$mpc."\t".$exacpli."\t".$remm."\t".$sutrstring.$spliceaistring.$calls."\n";
					} else {
						$printthis=1;
						my $consline = $csqs[ $fred ];
						$consline =~ s/^CSQ\=//; #get rid of the initial CSQ=, if it's there
						$printvcf = $printvcf.$consline.",";
					}
				}
			}
			$fred = $fred+1;
		}
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

