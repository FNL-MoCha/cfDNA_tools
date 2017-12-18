#!/usr/bin/perl
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Parallel::ForkManager;
use Data::Dump;
use Sort::Versions;

use constant DEBUG => 0;

my $version = "v0.3.121817";

my $gene;
my $threshold = 2;
my $ref_calls;
my $novel_calls;
my $nocall;
my $raw_output;

my $help;
my $ver_info;
my $outfile;

my $scriptname = basename($0);
my $description = <<"EOT";
Print out results of cfDNA fusion pipeline.  
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -g, --gene      Only output data for a specific driver gene or genes 
                    separated by a comma. 
    -t, --threshold Only report fusions above this threshold (DEFAULT: $threshold).
    -R, --Ref       Include reference variants too (DEFAULT: 'False').
    -n, --novel     Include 'Non-Targeted' fusions in the output (DEFAULT: 'False').
    -N, --NOCALL    Include NOCALL or FAIL Fusions in output (DEFAULT: 'False').
    -r, --raw       Raw output rather that pretty printed file.
    -o, --output    Write output to file.
    -v, --version   Display version information.
    -h, --help      Display this help text.
EOT

GetOptions( 
    "Ref|R"         => \$ref_calls,
    "novel|n"       => \$novel_calls,
    "threshold|t=s" => \$threshold,
    "gene|g=s"      => \$gene,
    "output|o=s"    => \$outfile,
    "raw|r"         => \$raw_output,
    "NOCALL|N"      => \$nocall,
    "help|h"        => \$help,
    "version|v"     => \$ver_info,
);

if (DEBUG) {
    print '='x25 . "  DEBUG: Thresholds as passed  " . '='x25 . "\n";
    print "\tMolecular Counts: $threshold\n";

    # Want to remap values for human readable formatting later.
    my ($nocall_txt, $ref_txt, $novel_txt);
    ($nocall) ? ($nocall_txt = 'True') : ($nocall_txt = 'False');
    ($ref_calls) ? ($ref_txt = 'True') : ($ref_txt = 'False');
    ($novel_calls) ? ($novel_txt = 'True') : ($novel_txt = 'False');

    print "\tOuput Reference Calls: $ref_txt\n";
    print "\tOuput Novel Calls: $novel_txt\n";
    print "\tOutput NOCALLs: $nocall_txt\n";
    print '='x82 . "\n";
}

sub help { 
    printf "%s - %s\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
    exit;
}

sub version_info {
    printf "%s - %s\n", $scriptname, $version;
    exit;
}

help if $help;
version_info if $ver_info;

# Check we have some files to process
if ( @ARGV < 1 ) {
    print "ERROR: You must load at least one VCF file\n";
    print $usage;
    exit 1;
}

# Set up output
my $out_fh;
if ( $outfile ) {
    print "Writing output to '$outfile'.\n";
    open( $out_fh, ">", $outfile ) 
        or die "Can't open the output file '$outfile' for writing: $!";
} else {
    $out_fh = \*STDOUT;
}

my @files = @ARGV;
my @genes_list = map{uc($_)} split(/,/, $gene) if $gene;

#######=================  END ARG Parsing  =================#######
my %results;

my $pm = new Parallel::ForkManager(48);
$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, 
            $data_structure_reference) = @_;
        my $vcf  = $data_structure_reference->{input};
        my $name = $data_structure_reference->{id};
        $name //= basename($vcf);
        $results{$$name} = $data_structure_reference->{result};
    }
);

for my $input_file ( @files ) {
    $pm->start and next;
    my ($return_data, $sample_id) = proc_vcf(\$input_file);
    $pm->finish(0,
        {
            result => $return_data,
            input  => $input_file,
            id     => $sample_id,
        }
    );
}
$pm->wait_all_children;

# Generate and print out the final results table(s)
select $out_fh;
for my $sample ( sort keys %results ) {
    unless ($raw_output) {
        print "::: ";
        print join(',', @genes_list) if @genes_list;
        print " Fusions in $sample :::\n\n";
    }
    
    if ( my $fusion_results = $results{$sample}{'fusions'} ) {
        my @keys = map{ join('.', (split(/\|/, $_))[0..1])  } keys %$fusion_results;
        my $fwidth = field_width(\@keys);
        
        my $format = "%-${fwidth}s %-12s %-12s %-15s %-15s\n";
        my @fusion_header = qw (Fusion ID Read_Count Driver_Gene Partner_Gene);
        printf $format, @fusion_header unless $raw_output;

        for my $entry (sort { versioncmp($a,$b) } keys %$fusion_results) {
            if ($gene) {
                next unless grep{ $$fusion_results{$entry}->{'DRIVER'} eq $_ } @genes_list;
            }

            my ($fusion, $junct, $id ) = split(/\|/, $entry);
            next if $$fusion_results{$entry}->{'COUNT'} < $threshold;

            print_data(\$sample, "$fusion.$junct", \$id, $$fusion_results{$entry}, 
                \$format);
        }
        print "\n";
    } else {
        print "\t\t\t<<< No Fusions Detected >>>\n\n" unless $raw_output;
    }
}
sub field_width {
    my $elems = shift;
    my @lens = sort{ $a <=> $b } map{length($_)} @$elems;
    return $lens[-1] + 2;
}
    
sub proc_vcf {
    my $vcf = shift;
    # Version 1,2, and 3 drivers. Not all exist in the current version, but keep all for backward compatibility.
    my @drivers = qw(ABL1 AKT2 AKT3 ALK AR AXL BRAF BRCA1 BRCA2 CDKN2A EGFR ERBB2 ERBB4 ERG ESR1 ETV1 ETV1a 
                     ETV1b ETV4 ETV4a ETV5 ETV5a ETV5d FGFR1 FGFR2 FGFR3 FGR FLT3 JAK2 KRAS MDM4 MET MYB MYBL1 
                     NF1 NOTCH1 NOTCH4 NRG1 NTRK1 NTRK2 NTRK3 NUTM1 PDGFRA PDGFRB PIK3CA PPARG PRKACA PRKACB 
                     PTEN RAD51B RAF1 RB1 RELA RET ROS1 RSPO2 RSPO3 TERT
    );
    my %results;
    (my $sample_name = $$vcf) =~ s/(:?_Fusion_filtered)?\.vcf$//i;
    $sample_name =~ s/_RNA//;

    open( my $in_fh, "<", $$vcf );
    while (<$in_fh>) {
        next if /^#/;
        my @data = split;
        if ( $data[7] =~ /SVTYPE=(Fusion|RNAExonVariant)/ ) {
            my ($count) = map { /MOL_COUNT=(\d+)/ } $data[7];
            
            # Get RNAExonVariant Controls and add to controls section of hash, 
            # but don't process further.
            if ($data[2] =~ /WT$/) {
                $results{'controls'}->{$data[2]} = $count and next;
            }

            # Filter out ref calls if we don't want to view them. Will now get
            # a 'FAIL' message if read counts too low (though might be able to
            # just rely on counts anyway?). Will do threshold filtering at print
            # step
            next if ($count == 0 or $data[6] eq 'FAIL') and ! $ref_calls;

            # Get rid of FAIL and NOCALL calls to be more compatible with MATCHBox output.
            next if $data[6] eq 'NOCALL' and ! $nocall;

            my ($pair, $junct, $id) = map{s/_\d//; $_} split(/\./, $data[2]);
            $id //= '-';
            my $fid = join('|', $pair, $junct, $id);

            # Get rid of Non-targeted fusions.
            next if (grep{ $id eq $_ } qw(Non-Targeted Novel)) and ! $novel_calls;

            my ($gene1, $gene2) = split(/-/, $pair);
            if ( $gene1 eq $gene2) {
                $results{'fusions'}->{$fid}{'DRIVER'} = $results{'fusions'}->{$fid}{'PARTNER'} = $gene1;
            }
            elsif (grep {$_ eq $gene1} @drivers) {
                $results{'fusions'}->{$fid}{'DRIVER'} = $gene1;
                $results{'fusions'}->{$fid}{'PARTNER'} = $gene2;
            }
            elsif (grep {$_ eq $gene2} @drivers) {
                $results{'fusions'}->{$fid}{'DRIVER'} = $gene2;
                $results{'fusions'}->{$fid}{'PARTNER'} = $gene1;
            }
            else {
                $results{'fusions'}->{$fid}{'DRIVER'} = 'UNKNOWN';
                $results{'fusions'}->{$fid}{'PARTNER'} = "$gene1,$gene2";
            }
            $results{'fusions'}->{$fid}{'COUNT'} = $count;
        }
        elsif ($data[7] =~ /SVTYPE=ProcControl/) {
            my ($gene, $counts) = $data[7] =~ /.*?GENE_NAME=(.*?);.*?MOL_COUNT=(\d+)/;
            $results{'controls'}->{$gene} = $counts;
        }
    }
    return \%results, \$sample_name;
}

sub print_data {
    my ($sample_name, $fusion_name, $id, $data, $format) = @_;

    if ($raw_output) {
        print join(',', $$sample_name, $fusion_name, $$id, $$data{'COUNT'}, $$data{'DRIVER'}), "\n";
    } else {
        printf $$format, $fusion_name, $$id, $$data{'COUNT'}, $$data{'DRIVER'}, $$data{'PARTNER'};
    }
    return;
}
