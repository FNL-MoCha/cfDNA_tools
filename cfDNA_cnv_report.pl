#!/usr/bin/perl
# Get CNV data from cfDNA panel VCF file
#################################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Sort::Versions;
use JSON -support_by_pp;
use Parallel::ForkManager;
use Data::Dump;
use Term::ANSIColor;

use constant DEBUG => 1;

my $scriptname = basename($0);
my $version = "v0.2.112117";

# Remove when in prod.
print "\n";
print colored("*" x 75, 'bold yellow on_black'), "\n";
print colored("      DEVELOPMENT VERSION OF $scriptname (version: $version)\n", 'bold yellow on_black');
print colored("*" x 75, 'bold yellow on_black');
print "\n\n";

my $description = <<"EOT";
cfDNA CNV parser.
EOT

my $copy_amp = 0;
my $copy_loss = 0;
my $fold_amp = 0;
my $fold_loss = 0;

my $outfile;
my $novel;
my $geneid;
my $tiles;
my $nocall;
my $format = 'pp';
my $raw_output;

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF_file(s)>
    Filter Options
    --amp         Only report copy gain above this threshold (DEFAULT: $copy_amp).
    --loss        Only report copy loss below this threshold (DEFAULT: $copy_loss).
    --fold_amp    Only report events with fold amplification above this 
                  threshold (DEFAULT: $fold_amp).
    --fold_loss   Only report events with fold loss below this threshold 
                  (DEFAULT: $fold_loss).
    -g, --gene    Print out results for this gene (or comma separated list of 
                  genes) only.
    -N, --NOCALL  Do not output NOCALL results (Default: OFF)

    Output Options
    -o, --output      Send output to custom file.  Default is STDOUT.
    -f, --format      Format to use for output.  Can choose 'csv', 'tsv', or pretty 
                      print (as 'pp') (DEFAULT: $format).
    -r, --raw         Output as a raw CSV format that can be input into Excel. The 
                      header output is not in columns as a part of the whole.
    -v, --version     Print version information
    -h, --help        Print this help information
EOT

my $help;
my $ver_info;

GetOptions( "novel|n"       => \$novel,
            "copy_amp=f"    => \$copy_amp,
            "copy_loss=f"   => \$copy_loss,
            "fold_amp=f"    => \$fold_amp,
            "fold_loss=f"   => \$fold_loss,
            "gene|g=s"      => \$geneid,
            "tiles|t=i"     => \$tiles,
            "output|o=s"    => \$outfile,
            "format|f=s"    => \$format,
            "raw|r"         => \$raw_output,
            "NOCALL|N"      => \$nocall,
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

my @genelist = split(/,/, $geneid) if $geneid;

# Do not allow both fold diff and copy number values to be used as filters. 
if (($copy_amp and $copy_loss) && ($fold_amp and $fold_loss)) {
    print "ERROR: You can not use both the fold difference and copy number threshold ",
        "for filtering. Please use just one or\nthe other.\n";
    exit 1;
}

my %filters = (
    'copy_amp'   => $copy_amp,
    'copy_loss'  => $copy_loss,
    'fold_amp'   => $fold_amp,
    'fold_loss'  => $fold_loss,
    'gene'       => [@genelist],
    'tiles'      => $tiles,
    'novel'      => $novel,
);

if (DEBUG) {
    print '='x35, '  DEBUG  ', '='x35, "\n";
    print "Filters being employed\n";
    while (my ($keys, $values) = each %filters) {
        $values //= 'undef';
        if ($keys eq 'gene') {
            printf "\t%-7s => %s\n",$keys,join(',',@$values);
        } else {
            printf "\t%-7s => %s\n",$keys,$values;
        }
    }
    print '='x79, "\n";
}

my %formats = (
    'csv'   => ',',
    'tsv'   => "\t",
    'pp'    => '',
);

# Set the output format delimiter
my $delimiter;
die "ERROR: '$format' is not a valid option as a delimiter!\n" unless defined $formats{$format};
($format) ? ($delimiter = $formats{$format}) : ($delimiter = $formats{pp});

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "ERROR: No VCF files passed to script!\n\n"; 
    print "$usage\n";
    exit 1;
}
my @vcfs = @ARGV;

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
    print "Writing data to $outfile.\n";
	open( $out_fh, ">", $outfile ) 
        || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}

#########----------------------- END ARG Parsing ---------------------#########
my %cnv_data;
my $pm = new Parallel::ForkManager(48);

$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my $vcf = $data_structure_reference->{input};
        my $name = $data_structure_reference->{id};
        $name //= basename($vcf);
        $cnv_data{$$name} = $data_structure_reference->{result};
    }
);

for my $input_file (@vcfs) {
    $pm->start and next;
    my ($return_data, $sample_id)  = proc_vcf(\$input_file);

    $pm->finish(0, 
        { 
          result  => $return_data, 
          input   => $input_file, 
          id      => $sample_id,
        }
     );
}
$pm->wait_all_children;

#dd \%cnv_data;
#__exit__(__LINE__,'Post read VCF file.');

my $results = proc_cnv_data(\%cnv_data, \%filters);
#dd $results;
#__exit__(__LINE__,'Post result processing.');

print_results($results, $delimiter);

sub proc_cnv_data {
    my ($cnv_data, $filters) = @_;
    my %results;

    for my $sample ( keys %$cnv_data ) {
        $results{$sample} = [];
        my @outfields = qw( END LEN NUMTILES CN FD HS FUNC PVAL RMMDP MMDP );

        for my $cnv ( sort { versioncmp ( $a, $b ) } keys %{$$cnv_data{$sample}} ) {
            my %mapped_cnv_data;
            last if $cnv eq 'NONE';

            my ($chr, $start, $gene, undef) = split( /:/, $cnv );
            %mapped_cnv_data = map{ $_ => $cnv_data{$sample}->{$cnv}->{$_} } @outfields;
            @mapped_cnv_data{qw(chr start gene undef)} = split(/:/, $cnv);
            $mapped_cnv_data{HS} //= 'No';

            # Get OVAT Annot Data
            # XXX: Remove this for now since there is no OVAT annotation at the moment.  Will bring back 
            # if they do re-implement this.
            #my ($gene_class, $variant_class);
            #my $func = $mapped_cnv_data{FUNC};
            #if ( $func && $func =~ /oncomine/ ) {
                #my $json_annot = JSON->new->allow_singlequote->decode($func);
                #my $parsed_annot = $$json_annot[0];
                #$gene_class = $$parsed_annot{'oncomineGeneClass'};
                #$variant_class = $$parsed_annot{'oncomineVariantClass'};
            #} else {
                #$gene_class = $variant_class = '---';
            #}
            #$mapped_cnv_data{GC} = $gene_class;
            #$mapped_cnv_data{VC} = $variant_class;

            my @filtered_data = filter_results(\%mapped_cnv_data, $filters);
            push(@{$results{$sample}}, \@filtered_data) if @filtered_data;
        }
    }
    return \%results;
}

sub print_results {
    my ($data, $delimiter) = @_;
    my @header = qw( Chr Gene Start End Length Tiles CN FD p-val Med_Mol_Cov 
        Med_Read_Cov );
    #my $pp_format = "%-8s %-8s %-11s %-11s %-11s %-8s %-8s %-8s %-8s %-8s %-8s %-18s\n";

    my %formatter = (
        'Chr'          => '%-8s',
        'Gene'         => '%-8s',
        'Start'        => '%-11s', 
        'End'          => '%-11s',
        'Length'       => '%-11s',
        'Tiles'        => '%-8s',
        'CN'           => '%-8s',
        'FD'           => '%-8s',
        'p-val'        => '%-8s',
        'Med_Mol_Cov'  => '%-14s',
        'Med_Read_Cov' => '%-14s',
    );

    select $out_fh;
    my $string_format = join(' ', @formatter{@header}) . "\n";

    # Print out comma separated dataset for easy import into Excel and whatnot.
    if ($raw_output) {
        raw_output($data, \@header);
    } else {
        for my $sample (keys %$data) {
            my ($id, $gender, $mapd, $cellularity) = split( /:/, $sample );
            print "::: CNV Data For $id (Gender: $gender, Cellularity: $cellularity, MAPD: $mapd) :::\n";
            ($delimiter) ? print join($delimiter, @header),"\n" : printf $string_format, @header;
            if ( ! @{$$data{$sample}} ) {
                print ">>>>  No Reportable CNVs Found in Sample  <<<<\n"; 
            } else {
                for my $cnv (@{$$data{$sample}}) {
                    ($delimiter) ? print join($delimiter, @$cnv), "\n" : printf $string_format, @$cnv;
                }
            }
            print "\n";
        }
    }
}

sub raw_output {
    # No fancy output, but want to have CNV data plus sample the extra metadata
    # output in a CSV format that we can read directly into Excel easily.
    my ($data, $header) = @_;

    select $out_fh;
    print join(',', qw(Sample Gender MAPD Cellularity), @$header), "\n";

    for my $sample (keys %$data) {
        my @elems = split(/:/, $sample);
        for my $cnv (@{$$data{$sample}}) {
            print join(',', @elems, @$cnv), "\n";
        }
    }
}

sub filter_results {
    # Filter out CNV data prior to printing it all out. Each call receives
    # a hash of data, and either None returned or a hash of specific fields
    # as defined in "return_data()"
    my ($data, $filters) = @_;
    my @cn_thresholds = @$filters{qw(copy_amp copy_loss fold_amp fold_loss)};

    # Gene level filter
    return if (@{$filters{gene}}) and ! grep {$$data{gene} eq $_} @{$filters{gene}};

    # We made it the whole way through; check for copy number thresholds
    (copy_number_filter($data, \@cn_thresholds)) ? return return_data($data) : return;
}

sub return_data {
    my $data = shift;
    my @fields = qw(chr gene start END LEN NUMTILES CN FD PVAL RMMDP MMDP);
    return @$data{@fields};
}

sub copy_number_filter {
    # Filter data based on either Fold Diff value or Copy Number value.
    my ($data, $threshold) = @_;
    my ($ca, $cl, $fa, $fl) = @$threshold;

    # XXX: Copy data is dependent on p-val. For now, let's just filter out anything
    # with a p-val > 10e-5 per TF QC doc.  Can tailor this a bit later.
    return 0 if $$data{'PVAL'} > 0.00005;

    # If we're using fold diff values, return 1 if FD > amp threshold or less 
    # than loss threshold.
    if ($fa and $fl) {
        return 1 if ($$data{'FD'} > $fa || $$data{'FD'} < $fl);
    }
    elsif ($ca and $cl) {
        return 1 if ($$data{'CN'} > $ca || $$data{'CN'} < $cl);
    }
    else {
        # Return everything if there are no filters.
        return 1;
    }

    # If we got here, the data does not meet threshold requirements and will be
    # filtered out.
    return 0;
}

sub proc_vcf {
    my $vcf = shift;
    my ($sample_id, $gender, $mapd, $cellularity, $sample_name);
    my %results;

    open( my $vcf_fh, "<", $$vcf);
    while (<$vcf_fh>) {
        if ( /^##/ ) {
            if ( $_ =~ /sampleGender=(\w+)/ ) {
                $gender = $1 and next;
            }
            # Need to add to accomodate the new CNV plugin; may not have the same field as the normal IR data.
            if ($_ =~ /AssumedGender=([mf])/) {
                ($1 eq 'm') ? ($gender='Male') : ($gender='Female');
                next;
            }
            elsif ( $_ =~ /mapd=(\d\.\d+)/ ) {
                $mapd = $1 and next;
            }
            elsif ( $_ =~ /CellularityAsAFractionBetween0-1=(.*)$/ ) {
                $cellularity = $1 and next;
            }
        } 

        my @data = split;
        if ( $data[0] =~ /^#/ ) {
            $sample_name = $data[-1] and next;
        }
        next unless $data[4] eq '<CNV>';
        $sample_id = join( ':', $sample_name, $gender, $mapd, $cellularity );

        # Let's handle NOCALLs for MATCHBox compatibility (prefer to filter on my own though).
        if ($nocall && $data[6] eq 'NOCALL') {
            ${$cnv_data{$sample_id}->{'NONE'}} = '';
            next;
        }

        my $varid = join( ':', @data[0..3] );
        
        # Kludgy, but need to deal with hotspots (HS) field; not like others!
        $data[7] =~ s/HS/HS=Yes/;
        $data[7] =~ s/SD;/SD=NA;/; # sometimes data in this field and sometimes not.  

        my @format = split( /;/, $data[7] );
        my ($cn) = $data[9] =~ /:([^:]+)$/;
        push( @format, "CN=$cn" );

        %{$results{$varid}} = map { split /=/ } @format;
    }
    if (DEBUG) {
        print "="x35, "  DEBUG  ", "="x35, "\n";
        print "\tSample Name:  $sample_name\n";
        print "\tCellularity:  $cellularity\n";
        print "\tGender:       $gender\n";
        print "\tMAPD:         $mapd\n";
        print "="x79, "\n";
    }
    return \%results, \$sample_id;
}

sub __exit__ {
    my ($line, $msg) = @_;
    print "\n\n";
    print colored("Got exit message at line: $line with message: $msg", 'bold white on_green');
    print "\n";
    exit;
}
