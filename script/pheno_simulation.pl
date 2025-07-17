#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $genotype_file = "";
my $heritability = 0.5;      
my $num_qtl = 100;           
my $ploidy = 6;              
my $num_phenotypes = 100;   
my $output_prefix = "simulation";

GetOptions(
    "genotype|g=s"     => \$genotype_file,
    "heritability|h=f" => \$heritability,
    "qtl|q=i"          => \$num_qtl,
    "ploidy|p=i"       => \$ploidy,
    "phenotypes|n=i"   => \$num_phenotypes,
    "output|o=s"       => \$output_prefix
) or die "Error!\nUsage: $0 -g <genotype_file> -h <h2> -q <QTL_number> -p <ploidy> -n <pheno_number> -o <output_prefix>\n";

die "Please define the genotype file (-g)\n" unless $genotype_file;
die "genotype file not found!: $genotype_file\n" unless -e $genotype_file;

print "Loading genotype file...\n";

my @individuals = ();  
my @snp_ids = ();      
my %genotypes = ();    
my %marker_info = ();  

open(my $fh, '<', $genotype_file) or die "file cannot be opened $genotype_file: $!\n";

my $header = <$fh>;
chomp $header;
my @header_parts = split(/\s+/, $header);

# Check the first 5 cols: Marker, Chrom, Position, REF, ALT
if (scalar(@header_parts) < 6) {
    die "Error: File format incorrect,（Marker, Chrom, Position, REF, ALT, Sample ...）\n";
}

my @expected_headers = ("Marker", "Chrom", "Position", "REF", "ALT");
for my $i (0..4) {
    if ($header_parts[$i] ne $expected_headers[$i]) {
        print "Warnning: " . ($i+1) . "Col title '$header_parts[$i]'，Expected '$expected_headers[$i]'\n";
    }
}

@individuals = @header_parts[5..$#header_parts];
print "Checked " . scalar(@individuals) . " 个个体\n";
print "Individual: " . join(", ", @individuals[0..min(4, $#individuals)]);
print "..." if scalar(@individuals) > 5;
print "\n";

    chomp $line;
    next if $line =~ /^\s*$/;  
    
    my @parts = split(/\s+/, $line);
    
    if (scalar(@parts) != scalar(@header_parts)) {
        print "Warning: The" . ($snp_count + 2) . "rows and cols not matched，next line\n";
        next;
    }
    
    my $marker_id = $parts[0];
    my $chrom = $parts[1];
    my $position = $parts[2];
    my $ref_allele = $parts[3];
    my $alt_allele = $parts[4];
    
    $marker_info{$marker_id} = {
        chrom => $chrom,
        position => $position,
        ref => $ref_allele,
        alt => $alt_allele
    };
    
    push @snp_ids, $marker_id;
    
    for my $i (0..$#individuals) {
        my $genotype = $parts[$i + 5];  
        
        
        if (!defined $genotype || $genotype eq 'NA' || $genotype eq '') {
            $genotype = 0;
        } elsif ($genotype !~ /^\d+$/) {
            print "Warning: marker $marker_id individual $individuals[$i] genotype value '$genotype' not number，set 0\n";
            $genotype = 0;
        }
        
        $genotypes{$marker_id}{$individuals[$i]} = $genotype;
    }
    $snp_count++;
}
close($fh);

print "Loading $snp_count  SNPs\n";

my $max_dosage = 0;
my $dosage_distribution = {};  

for my $snp_id (@snp_ids) {
    for my $ind (@individuals) {
        my $dosage = $genotypes{$snp_id}{$ind};
        $max_dosage = $dosage if $dosage > $max_dosage;
        $dosage_distribution->{$dosage}++;
    }
}

print "gene dosage distribution: ";
for my $dosage (sort {$a <=> $b} keys %$dosage_distribution) {
    print "$dosage(" . $dosage_distribution->{$dosage} . ") ";
}
print "\n";

if ($max_dosage > $ploidy) {
    print "警告: 检测到的最大基因剂量 ($max_dosage) 超过设定的倍性 ($ploidy)\n";
    print "将倍性调整为 $max_dosage\n";
    $ploidy = $max_dosage;
}

print "使用倍性: $ploidy\n";

# 随机选择QTL
print "随机选择 $num_qtl 个QTL...\n";
die "QTL数量不能超过SNP总数\n" if $num_qtl > $snp_count;

# Fisher-Yates select QTLs
my @snp_indices = (0..$#snp_ids);
for my $i (reverse 0..$#snp_indices) {
    my $j = int(rand($i + 1));
    ($snp_indices[$i], $snp_indices[$j]) = ($snp_indices[$j], $snp_indices[$i]);
}

my @selected_qtl_indices = @snp_indices[0..($num_qtl-1)];
my @qtl_snps = @snp_ids[@selected_qtl_indices];

print "已选择QTL位点\n";
print "示例QTL: " . join(", ", @qtl_snps[0..min(4, $#qtl_snps)]) . "\n";

print "为QTL分配效应值...\n";
my %qtl_effects = ();

for my $qtl (@qtl_snps) {
    # Box-Muller transform
    my $u1 = rand();
    my $u2 = rand();
    my $z0 = sqrt(-2 * log($u1)) * cos(2 * 3.14159265359 * $u2);
    $qtl_effects{$qtl} = $z0;
}

my $sum_var = 0;
for my $ind (@individuals) {
    my $genetic_value = 0;
    for my $qtl (@qtl_snps) {
        my $dosage = $genotypes{$qtl}{$ind};
        my $normalized_dosage = $dosage / $ploidy;
        $genetic_value += $qtl_effects{$qtl} * $normalized_dosage;
    }
    $sum_var += $genetic_value ** 2;
}

my $genetic_var = $sum_var / scalar(@individuals);
my $scaling_factor = $genetic_var > 0 ? sqrt($heritability / $genetic_var) : 1;

for my $qtl (@qtl_snps) {
    $qtl_effects{$qtl} *= $scaling_factor;
}

print "QTL效应值已标准化\n";

print "生成 $num_phenotypes 组表型数据...\n";
my @all_phenotypes = ();

for my $rep (1..$num_phenotypes) {
    my %phenotypes = ();
    
    for my $ind (@individuals) {
        my $genetic_value = 0;
        for my $qtl (@qtl_snps) {
            my $dosage = $genotypes{$qtl}{$ind};
            my $normalized_dosage = $dosage / $ploidy;
            $genetic_value += $qtl_effects{$qtl} * $normalized_dosage;
        }
         my $error_var = (1 - $heritability);
        # 
        my $u1 = rand();
        my $u2 = rand();
        my $environmental_effect = sqrt(-2 * log($u1)) * cos(2 * 3.14159265359 * $u2) * sqrt($error_var);
        
        $phenotypes{$ind} = $genetic_value + $environmental_effect;
    }
    
    push @all_phenotypes, \%phenotypes;
}


my $qtl_output = "${output_prefix}_QTL_effects.txt";
print "Output QTL effects: $qtl_output\n";

open(my $qtl_fh, '>', $qtl_output) or die "cannot create file $qtl_output: $!\n";
print $qtl_fh "QTL_ID\tChrom\tPosition\tREF\tALT\tEffect\n";
for my $qtl (@qtl_snps) {
    printf $qtl_fh "%s\t%s\t%s\t%s\t%s\t%.6f\n", 
           $qtl, 
           $marker_info{$qtl}{chrom}, 
           $marker_info{$qtl}{position},
           $marker_info{$qtl}{ref},
           $marker_info{$qtl}{alt},
           $qtl_effects{$qtl};
}
close($qtl_fh);

my $phenotype_output = "${output_prefix}_phenotypes.txt";
print "Ouput phenotype: $phenotype_output\n";

open(my $phen_fh, '>', $phenotype_output) or die "cannot create file $phenotype_output: $!\n";

print $phen_fh "Individual";
for my $rep (1..$num_phenotypes) {
    print $phen_fh "\tPhenotype_$rep";
}
print $phen_fh "\n";
for my $ind (@individuals) {
    print $phen_fh $ind;
    for my $rep (0..($num_phenotypes-1)) {
        printf $phen_fh "\t%.6f", $all_phenotypes[$rep]{$ind};
    }
    print $phen_fh "\n";
}
close($phen_fh);

print "Individual number: " . scalar(@individuals) . "\n";
print "SNP number: $snp_count\n";
print "QTL number: $num_qtl\n";
print "h2: $heritability\n";
print "ploidy: $ploidy\n";
print "Pheno number: $num_phenotypes\n";
print "QTL effect file: $qtl_output\n";
print "Pheno output: $phenotype_output\n";

my $total_genetic_var = 0;
my $total_phenotypic_var = 0;

for my $rep (0..($num_phenotypes-1)) {
    my @genetic_values = ();
    my @phenotype_values = ();
    
    for my $ind (@individuals) {
        my $genetic_value = 0;
        for my $qtl (@qtl_snps) {
            my $dosage = $genotypes{$qtl}{$ind};
            my $normalized_dosage = $dosage / $ploidy;
            $genetic_value += $qtl_effects{$qtl} * $normalized_dosage;
        }
        push @genetic_values, $genetic_value;
        push @phenotype_values, $all_phenotypes[$rep]{$ind};
    }
    
    my $genetic_mean = 0;
    my $phenotype_mean = 0;
    $genetic_mean += $_ for @genetic_values;
    $phenotype_mean += $_ for @phenotype_values;
    $genetic_mean /= @genetic_values;
    $phenotype_mean /= @phenotype_values;
    
    my $genetic_var_rep = 0;
    my $phenotype_var_rep = 0;
    for my $i (0..$#genetic_values) {
        $genetic_var_rep += ($genetic_values[$i] - $genetic_mean) ** 2;
        $phenotype_var_rep += ($phenotype_values[$i] - $phenotype_mean) ** 2;
    }
    $genetic_var_rep /= (@genetic_values - 1);
    $phenotype_var_rep /= (@phenotype_values - 1);
    
    $total_genetic_var += $genetic_var_rep;
    $total_phenotypic_var += $phenotype_var_rep;
}

my $avg_genetic_var = $total_genetic_var / $num_phenotypes;
my $avg_phenotypic_var = $total_phenotypic_var / $num_phenotypes;
my $realized_heritability = $avg_phenotypic_var > 0 ? $avg_genetic_var / $avg_phenotypic_var : 0;

printf "h2: %.4f (target: %.4f)\n", $realized_heritability, $heritability;

sub min {
    my ($a, $b) = @_;
    return $a < $b ? $a : $b;
}
