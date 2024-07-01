#!/usr/bin/perl -w
#
#
#------------------------------------------------------------------------------
# TrackHub Generator
#------------------------------------------------------------------------------
# TrackHub Generator allows the automatic generation of track hubs from a collection of BED files.
# To run the TrackHub Generator enter the following command prompt in a unix shell:
# perl TrackHubGenerator.pl PATH/TO/NANE ASSEMBLY FBED UCSC EMAIL
#
# Positional arguments:
# PATH/TO/NAME	Path to location where track hub will be created ending with the name of the track hub
# ASSEMBYL		Assembly of the genome for which the track hub shall be used, e.g. hg38
# FBED			Path to the folder containing all BED files which should be included in the track hub. Each BED file will result in a separate track
# UCSC			Path to the filder containing the bedTOBigBed binary and the fetchChromSizes script downloaded from UCSC
# EMAIL			Email address to be used as contact for people using the track hub
#
#
# Dependencies:
# The TrackHub Generator relies on two tools provided through UCSC Genomes which can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/
# - bedToBigBed			binary tool to create binary BigBed files from input bed files
# - fetchChromSizes		script to download chromosome sizes for indexing BigBed files during generation through bedToBigBed
#
#
# For more information please refere to www.sanger.ac.uk/science/tools/trackhub-generator
#------------------------------------------------------------------------------

use File::Basename;
use File::Find::Rule;

my $hubpath = shift;
my $genome = shift;
my $bedfiles = shift;
my $toolfiles = shift;
my $email = shift;

my $hubname = basename($hubpath);
my $name = $hubname;
$name =~ tr/ //ds;

my $dirname = dirname($hubpath);
$hubpath = join("",$dirname,"/",$hubname,"/");
my $shortlabel = substr($hubname,0,17);
my $longlabel = substr($hubname,0,80);
my $genomepath = join("",$hubpath,$genome);


if (!-d $hubpath) {
	mkdir $hubpath or die "Failed to create path for hub: $hubpath";
	mkdir $genomepath or die "Failed to create path for genome in hub path: $genomepath";
} else {
	if(!-d $genomepath) {
		mkdir $genomepath or die "Failed to create path for genome in hub path: $genomepath";
	}
}


if(-d $hubpath) {

	my $strlength = length $hubname;
	my $hubfile = join("",$hubpath,"hub.txt");
	open(OUT, '>', $hubfile) or die "Could not open file $hubfile";
	print OUT "hub ${name}\nshortLabel ${shortlabel}\nlongLabel ${longlabel}\ngenomesFile genomes.txt\nemail ${email}";
	close OUT;

	my @genomes = File::Find::Rule->new->directory->mindepth(1)->maxdepth(1)->in($hubpath);
	my $temp ="";

	my $genomesfile = join("",$hubpath,"genomes.txt");
	open(OUT, '>', $genomesfile) or die "Could not open file $genomesfile";
	foreach(@genomes) {
		$temp = basename($_);
		print OUT "genome $temp\ntrackDb $temp/trackDb.txt\n\n";
	}
	close OUT;


	my $dir = dirname($toolfiles);
	my $base = basename($toolfiles);
	if(-d $toolfiles) {
    	my $fetchChromSizes = join("",$toolfiles,"fetchChromSizes");
    	my $chromsizes = join("",$toolfiles,$genome,".chrom.sizes");
    	if(!-f $chromsizes) {
        	if(-f $fetchChromSizes) {
            	print("Fetching chromosome sizes for $genome from UCSC\n");
            	`${fetchChromSizes} ${genome} > ${chromsizes}`;
        	} else {
            	die "Tool \"fetchChromSizes\" not found: $fetchChromSizes";
        	}
    	}
	} else {
    	die "There is no directory $toolfiles";
	}




	my $bedtobigbedfile = join("",$toolfiles,"bedToBigBed");
	if(-f $bedtobigbedfile) {
		my $trackDbstring = "";
		$dir = dirname($bedfiles);
		$base = basename($bedfiles);
		$bedfiles = join("",$dir,"/",$base,"/");
		my $chromsizes = join("",$toolfiles,$genome,".chrom.sizes");
		my @bedfiles = glob "$bedfiles*.bed";
		for my $file (@bedfiles) {
			my $newtempbed = join("",$file,".2");
			my $result = `sort -k1,1 -k2,2n "${file}" > "${newtempbed}"`;
			my $filename = basename($file);
			$filename =~ s/\.bed//g;
			my $shortfilename = substr($filename,0,17);
			$filename = substr($filename,0,80);
			my $bbfile = join("",$genomepath,"/",$filename,".bb");
			print "Converting $file to bigBed format\n";
			$result = `${bedtobigbedfile} "${newtempbed}" ${chromsizes} "${bbfile}"`;
			$result = `rm "${newtempbed}"`;
			my $trackname = $filename;
			$trackname =~ tr/ //ds;
			$trackDbstring = $trackDbstring . "track ".$trackname."\nbigDataUrl ".$filename.".bb\nshortLabel ".$shortfilename."\nlongLabel ".$filename."\ntype bigBed 12\nvisibility dense\nitemRgb on\n\n";
		}

		my $trackDbfile = join("",$genomepath,"/","trackDb.txt");
		open(OUT, '>', $trackDbfile) or die "Could not open file $trackDbfile";
		print OUT $trackDbstring;
		close OUT;
	} else {
		die "Tool \"bedToBigBed\" not found: $bedtobigbedfile";
	}
} else {
	die "No path for track hub avaliable: $hubpath";
}