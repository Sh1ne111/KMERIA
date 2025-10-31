#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

#Copyright to Chen Shuai (chensss1209@gmail.com) 
#Date:  2025-10-25

our $VERSION = "1.0.0";

my $action = "";
my $sample_list = "";
my $depth_file = "";
my $pheno_file = "";
my $output_depth = "";
my $output_pheno = "";
my $help = 0;

GetOptions(
    "action=s"       => \$action,
    "samples=s"      => \$sample_list,
    "depth=s"        => \$depth_file,
    "pheno=s"        => \$pheno_file,
    "out-depth=s"    => \$output_depth,
    "out-pheno=s"    => \$output_pheno,
    "help|h"         => \$help,
) or pod2usage(2);

pod2usage(-verbose => 2) if $help;

if ($action eq "check") {
    check_sample_order();
} elsif ($action eq "reorder") {
    reorder_files();
} elsif ($action eq "extract") {
    extract_from_depth();
} else {
    print "Error: Unknown action '$action'. Valid actions are: check, reorder, extract\n";
    pod2usage(1);
}

sub check_sample_order {
    unless ($sample_list) {
        die "Error: --samples is required for 'check' action\n";
    }
    
    print "=" x 70 . "\n";
    print "Checking Sample Order Consistency\n";
    print "=" x 70 . "\n\n";
    
    # Read sample list
    my @samples = read_file_samples($sample_list);
    print "Reference sample list: $sample_list\n";
    print "  Total samples: " . scalar(@samples) . "\n";
    print "  First 5: " . join(", ", @samples[0..4]) . "\n\n";
    
    my $all_consistent = 1;
    
    # Check depth file
    if ($depth_file && -f $depth_file) {
        print "Checking depth file: $depth_file\n";
        my @depth_samples = read_file_samples($depth_file);
        my $consistent = compare_sample_order(\@samples, \@depth_samples, "depth file");
        $all_consistent = 0 unless $consistent;
        print "\n";
    } else {
        print "Depth file not provided or not found\n\n";
    }
    
    # Check phenotype file
    if ($pheno_file && -f $pheno_file) {
        print "Checking phenotype file: $pheno_file\n";
        my @pheno_samples = read_file_samples($pheno_file);
        my $consistent = compare_sample_order(\@samples, \@pheno_samples, "phenotype file");
        $all_consistent = 0 unless $consistent;
        print "\n";
    } else {
        print "Phenotype file not provided or not found\n\n";
    }
    
    print "=" x 70 . "\n";
    if ($all_consistent) {
        print "✓ All files have CONSISTENT sample order!\n";
        print "=" x 70 . "\n";
    } else {
        print "✗ WARNING: Sample order INCONSISTENCY detected!\n";
        print "  Use 'reorder' action to fix the order\n";
        print "=" x 70 . "\n";
        exit 1;
    }
}

sub reorder_files {
    unless ($sample_list) {
        die "Error: --samples is required for 'reorder' action\n";
    }
    
    # Read reference sample order
    my @samples = read_file_samples($sample_list);
    print "Reference sample order: $sample_list (" . scalar(@samples) . " samples)\n\n";
    
    # Reorder depth file
    if ($depth_file && -f $depth_file) {
        unless ($output_depth) {
            $output_depth = $depth_file;
            $output_depth =~ s/\.tsv$/_reordered.tsv/;
            $output_depth .= "_reordered.tsv" unless $output_depth =~ /_reordered/;
        }
        
        print "Reordering depth file...\n";
        print "  Input:  $depth_file\n";
        print "  Output: $output_depth\n";
        
        my %depth_data = read_depth_file($depth_file);
        write_depth_file($output_depth, \@samples, \%depth_data);
        
        print "  ✓ Depth file reordered successfully\n\n";
    }
    
    # Reorder phenotype file
    if ($pheno_file && -f $pheno_file) {
        unless ($output_pheno) {
            $output_pheno = $pheno_file;
            $output_pheno =~ s/\.tsv$/_reordered.tsv/;
            $output_pheno =~ s/\.txt$/_reordered.txt/;
            $output_pheno .= "_reordered.txt" unless $output_pheno =~ /_reordered/;
        }
        
        print "Reordering phenotype file...\n";
        print "  Input:  $pheno_file\n";
        print "  Output: $output_pheno\n";
        
        my %pheno_data = read_pheno_file($pheno_file);
        write_pheno_file($output_pheno, \@samples, \%pheno_data);
        
        print "  ✓ Phenotype file reordered successfully\n\n";
    }
    
    print "Reordering complete!\n";
    print "Please verify the output files before using them in the pipeline.\n";
}

sub extract_from_depth {
    unless ($depth_file && -f $depth_file) {
        die "Error: --depth file is required and must exist for 'extract' action\n";
    }
    
    unless ($sample_list) {
        $sample_list = "sample_list_from_depth.txt";
    }
    
    print "Extracting sample list from depth file...\n";
    print "  Input:  $depth_file\n";
    print "  Output: $sample_list\n\n";
    
    my @samples = read_file_samples($depth_file);
    
    open(my $fh, '>', $sample_list) or die "Cannot create $sample_list: $!";
    foreach my $sample (@samples) {
        print $fh "$sample\n";
    }
    close($fh);
    
    print "✓ Extracted " . scalar(@samples) . " samples\n";
    print "✓ Sample list saved to: $sample_list\n";
}


sub read_file_samples {
    my ($file) = @_;
    my @samples;
    
    open(my $fh, '<', $file) or die "Cannot open $file: $!";
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;  # Skip empty lines
        next if $line =~ /^#/;     # Skip comments
        
        # Extract first column (sample name)
        my @fields = split(/\s+/, $line);
        push @samples, $fields[0];
    }
    close($fh);
    
    return @samples;
}

sub compare_sample_order {
    my ($ref_samples, $test_samples, $label) = @_;
    
    my $ref_count = scalar(@$ref_samples);
    my $test_count = scalar(@$test_samples);
    
    print "  Total samples in $label: $test_count\n";
    
    if ($ref_count != $test_count) {
        print "  ✗ ERROR: Sample count mismatch!\n";
        print "    Reference: $ref_count samples\n";
        print "    $label: $test_count samples\n";
        
        my %ref_hash = map { $_ => 1 } @$ref_samples;
        my %test_hash = map { $_ => 1 } @$test_samples;
        
        my @missing_in_test = grep { !exists $test_hash{$_} } @$ref_samples;
        my @extra_in_test = grep { !exists $ref_hash{$_} } @$test_samples;
        
        if (@missing_in_test) {
            print "    Missing in $label: " . join(", ", @missing_in_test[0..4]) . "\n";
        }
        if (@extra_in_test) {
            print "    Extra in $label: " . join(", ", @extra_in_test[0..4]) . "\n";
        }
        
        return 0;
    }
    
    my $order_consistent = 1;
    my @mismatches;
    
    for (my $i = 0; $i < $ref_count; $i++) {
        if ($ref_samples->[$i] ne $test_samples->[$i]) {
            $order_consistent = 0;
            push @mismatches, {
                pos => $i,
                ref => $ref_samples->[$i],
                test => $test_samples->[$i]
            };
            last if scalar(@mismatches) >= 5;  # Only show first 5 mismatches
        }
    }
    
    if ($order_consistent) {
        print "  ✓ Sample order is CONSISTENT\n";
        return 1;
    } else {
        print "  ✗ WARNING: Sample order is INCONSISTENT\n";
        print "    First mismatches:\n";
        foreach my $mm (@mismatches) {
            printf("      Position %d: '%s' (reference) vs '%s' ($label)\n", 
                   $mm->{pos}, $mm->{ref}, $mm->{test});
        }
        return 0;
    }
}

sub read_depth_file {
    my ($file) = @_;
    my %data;
    
    open(my $fh, '<', $file) or die "Cannot open $file: $!";
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        next if $line =~ /^#/;
        
        my @fields = split(/\t/, $line);
        if (scalar(@fields) >= 2) {
            $data{$fields[0]} = $fields[1];
        }
    }
    close($fh);
    
    return %data;
}

sub write_depth_file {
    my ($file, $samples, $data) = @_;
    
    open(my $fh, '>', $file) or die "Cannot create $file: $!";
    foreach my $sample (@$samples) {
        if (exists $data->{$sample}) {
            print $fh "$sample\t$data->{$sample}\n";
        } else {
            warn "WARNING: No depth data found for sample: $sample\n";
            print $fh "$sample\tNA\n";
        }
    }
    close($fh);
}

sub read_pheno_file {
    my ($file) = @_;
    my %data;
    
    open(my $fh, '<', $file) or die "Cannot open $file: $!";
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        next if $line =~ /^#/;
        
        my @fields = split(/\s+/, $line);
        if (scalar(@fields) >= 1) {
            my $sample = shift @fields;
            $data{$sample} = join("\t", @fields);
        }
    }
    close($fh);
    
    return %data;
}

sub write_pheno_file {
    my ($file, $samples, $data) = @_;
    
    open(my $fh, '>', $file) or die "Cannot create $file: $!";
    foreach my $sample (@$samples) {
        if (exists $data->{$sample}) {
            print $fh "$sample\t$data->{$sample}\n";
        } else {
            warn "WARNING: No phenotype data found for sample: $sample\n";
            print $fh "$sample\tNA\n";
        }
    }
    close($fh);
}

__END__

=head1 NAME

sample_order_manager.pl - Manage and verify sample order consistency

=head1 SYNOPSIS

sample_order_manager.pl --action <action> [options]

 Actions:
   check      Check if sample order is consistent across files
   reorder    Reorder depth/phenotype files to match reference sample list
   extract    Extract sample list from depth file

 Options:
   --samples FILE      Reference sample list file (required for most actions)
   --depth FILE        Sample depth file (sample_depth.tsv)
   --pheno FILE        Phenotype file (sample_pheno.tsv)
   --out-depth FILE    Output reordered depth file
   --out-pheno FILE    Output reordered phenotype file
   --help|-h           Show this help message

=head1 DESCRIPTION

This script helps ensure sample order consistency across all input files
used in the KMERIA pipeline. Sample order inconsistency is a common source
of errors in genotype-phenotype association studies.

=head1 EXAMPLES

=head2 Check sample order consistency

perl sample_order_manager.pl --action check \
  --samples sample.list \
  --depth sample_depth.tsv \
  --pheno sample_pheno.tsv

=head2 Reorder files to match reference sample list

perl sample_order_manager.pl --action reorder \
  --samples sample.list \
  --depth sample_depth.tsv \
  --pheno sample_pheno.tsv \
  --out-depth sample_depth_reordered.tsv \
  --out-pheno sample_pheno_reordered.tsv

=head2 Extract sample list from depth file

perl sample_order_manager.pl --action extract \
  --depth sample_depth.tsv \
  --samples output_sample.list

=head1 WORKFLOW RECOMMENDATION

Before running KMERIA pipeline:

1. Decide on your reference sample order (usually from sample.list)
2. Check consistency:
   perl sample_order_manager.pl --action check --samples sample.list \
     --depth sample_depth.tsv --pheno sample_pheno.tsv

3. If inconsistent, reorder files:
   perl sample_order_manager.pl --action reorder --samples sample.list \
     --depth sample_depth.tsv --pheno sample_pheno.tsv

4. Verify the reordered files:
   perl sample_order_manager.pl --action check --samples sample.list \
     --depth sample_depth_reordered.tsv --pheno sample_pheno_reordered.tsv

5. Run KMERIA pipeline with reordered files

=head1 FILE FORMATS

=head2 Sample list (--samples)

Plain text, one sample per line:
  sample1
  sample2
  sample3

=head2 Depth file (--depth)

Tab-separated, no header:
  sample1    45.2
  sample2    52.8
  sample3    38.9

=head2 Phenotype file (--pheno)

Tab or space-separated, no header:
  sample1  1.5  0
  sample2  2.3  1
  sample3  1.8  1

=head1 AUTHOR

Version 1.0.0 - Sample order management tool for KMERIA pipeline

=cut
