#!/usr/bin/perl
# Oncoming PanCancer cfDNA pipeline has header formatting problems (and some 
# others too!) that are preventing adequate parsing of the VCF files.  Fix the 
# header to avoid the VCF Tools parsing warnings.
#
# 10/4/2017 - D Sims
################################################################################
use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use Data::Dump;

my $scriptname = basename($0);
my $version = "v1.0.111517";
my $description = <<"EOT";
Current Thermo Fisher cfDNA VCF files are not conforming to VCF standard and throw
warnings whenever they are processed by standard VCF processing pipelines, such
as VCF Tools.  Fix the VCF header, and any aberrant entries, so that it will 
not throw warnings and errors when processed with standard tools.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF>
    -o, --output    Send output to custom file.  Default is <VCF>_fixed.vcf
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;

GetOptions( "output|o=s"    => \$outfile,
            "version|v"     => \$ver_info,
            "help|h"        => \$help )
        or die $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version if $ver_info;

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "ERROR: Not enough arguments passed to script!\n\n";
    print "$usage\n";
    exit 1;
}
my $vcf = shift;

# Set up the output file handle.
my $new_vcf;
($outfile) ? ($new_vcf = $outfile) : (($new_vcf = $vcf) =~ s/\.vcf/_fixed.vcf/);
open( my $out_fh, ">", $new_vcf); 

#########--------------------- END ARG Parsing ------------------------#########
sub read_vcf {
    my ($vcf, $out_fh) = @_;
    my (@header, @vars);
    my @bad_fields = qw(MOL_RATIO_TO_WILD_TYPE NORM_COUNT_WITHIN_GENE 
        RATIO_TO_WILD_TYPE NORM_MOL_COUNT_WITHIN_GENE);

    open(my $fh, "<", $vcf);
    while (my $line = <$fh>) {
        if ($line =~ /^#/) {
            if (grep {$line =~ /$_/} @bad_fields) {
                push(@header, fix_line($line));
            } else {
                push(@header, $line);
            }
        } else {
            push(@vars,$line);
        }
    }
    print "Done!\n";
    write_vcf($vcf, $out_fh, \@header, \@vars);
}

sub fix_line {
    chomp(my $line = shift);
    (my $trimmed_line = $line) =~ s/##INFO=<(.*)>/$1/;
    my @pairs = $trimmed_line =~ /(\w+=[^,]+)/g;
    $pairs[-1] =~ s/ "$/."/;
    return "##INFO=<" . join(',',@pairs) . ">\n";
}

sub write_vcf {
    my ($vcf, $out_fh, $header, $vars) = @_;

    print "Writing new VCF file...";
    print {$out_fh} $_ for @$header;
    print {$out_fh} $_ for @$vars;
    print "Done!\n";
}

print "Reading VCF file and extracting header and bad lines...";
read_vcf($vcf, $out_fh);
