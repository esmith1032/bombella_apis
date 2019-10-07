#!/usr/bin/perl
##Eric Smith

use strict;
use warnings;

my $base_dir = "/Research/postdoc/p_apium/orthoMCL/glox_outgroup/all_genomes_new_paras_added_bombella";
my $proteome_dir = "${base_dir}/proteomes";
my %prot_counts;
my %types_counts;
my $present_list = '';
my $species_present = 0;
my $total_species = 0;
my @proteomes = glob("$proteome_dir/*.fasta");

my $groups_file = "${base_dir}/filtered_groups.txt";

my $singles_file = "${base_dir}/single_copy_orthologs.txt";
my $multi_file = "${base_dir}/multi_copy_orthologs.txt";
my $var_file = "${base_dir}/variable_copy_orthologs.txt";
my $singletons_file = "${base_dir}/singletons.txt";

my $summary_file = "${base_dir}/orthomcl_filtered_summary.txt";

foreach my $prot (@proteomes) {
  next if $prot =~ /glox/;
  $total_species++;
  $prot =~ /\S+\/(\w+)\.fasta/;
  my $taxon = $1;
  $prot_counts{$taxon} = 0;
  $types_counts{sco}{$taxon}{cogs} = 0;
  $types_counts{mco}{$taxon}{cogs} = 0;
  $types_counts{var}{$taxon}{cogs} = 0;
  $types_counts{ton}{$taxon}{cogs} = 0;
  $types_counts{sco}{$taxon}{genes} = 0;
  $types_counts{mco}{$taxon}{genes} = 0;
  $types_counts{var}{$taxon}{genes} = 0;
  $types_counts{ton}{$taxon}{genes} = 0;
}
$types_counts{sco}{total}{cogs} = 0;
$types_counts{sco}{total}{genes} = 0;
$types_counts{mco}{total}{cogs} = 0;
$types_counts{mco}{total}{genes} = 0;
$types_counts{var}{total}{cogs} = 0;
$types_counts{var}{total}{genes} = 0;
$types_counts{ton}{total}{cogs} = 0;
$types_counts{ton}{total}{genes} = 0;

open SCO, ">$singles_file";
open MCO, ">$multi_file";
open VAR, ">$var_file";
open TON, ">$singletons_file";
open GROUPS, "$groups_file";

while (my $line = <GROUPS>) {
  chomp $line;
  
  my ($group, $gene_list) = split /:/, $line;
  my @genes = split /\s+/, $gene_list;
  shift @genes;
  
  foreach my $gene_id (@genes) {
    my ($species, $id) = split /\|/, $gene_id;
    $prot_counts{$species}++;
    unless ($present_list =~ /$species/) {
      $present_list .= $species;
      $species_present++;
    }
  }
  
  if ($species_present == 1) {
    print TON "$line\n";
    $types_counts{ton}{$present_list}{cogs}++;
    $types_counts{ton}{$present_list}{genes} += $prot_counts{$present_list};
    $types_counts{ton}{total}{cogs}++;
    $types_counts{ton}{total}{genes} += $prot_counts{$present_list};
  } elsif ($species_present < $total_species) {
    print VAR "$line\n";
    $types_counts{var}{total}{cogs}++;
    foreach my $test_species (keys %prot_counts) {
      if ($present_list =~ /$test_species/) {
        $types_counts{var}{total}{genes} += $prot_counts{$test_species};
        $types_counts{var}{$test_species}{cogs}++;
        $types_counts{var}{$test_species}{genes} += $prot_counts{$test_species};
      }
    }
  } else {
    my $sco = "yes";
    foreach my $test_species (keys %prot_counts) {
      $sco = "no" if $prot_counts{$test_species} > 1;
    }
    if ($sco eq "yes") {
      print SCO "$line\n";
      $types_counts{sco}{total}{cogs}++;
      foreach my $species (keys %prot_counts) {
        $types_counts{sco}{total}{genes} += $prot_counts{$species};
        $types_counts{sco}{$species}{cogs}++;
        $types_counts{sco}{$species}{genes} += $prot_counts{$species};
      }
    } else {
      print MCO "$line\n";
      $types_counts{mco}{total}{cogs}++;
      foreach my $species (keys %prot_counts) {
        $types_counts{mco}{total}{genes} += $prot_counts{$species};
        $types_counts{mco}{$species}{cogs}++;
        $types_counts{mco}{$species}{genes} += $prot_counts{$species};
      }
    }
  }
  $present_list = '';
  $species_present = 0;
  foreach my $reset_species (keys %prot_counts) {
    $prot_counts{$reset_species} = 0;
  }
}

close GROUPS;
close TON;
close VAR;
close MCO;
close SCO;

open SUM, ">$summary_file";
print SUM "type\ttotal_cogs\ttotal_genes";

foreach my $print_species (sort keys %prot_counts) {
  print SUM "\t${print_species}_cogs\t${print_species}_genes";
}

print SUM "\n";

foreach my $type (keys %types_counts) {
  print SUM "$type\t$types_counts{$type}{total}{cogs}\t$types_counts{$type}{total}{genes}";
  foreach my $count_species (sort keys %prot_counts) {
    print SUM "\t$types_counts{$type}{$count_species}{cogs}\t$types_counts{$type}{$count_species}{genes}";
  }
  print SUM "\n";
}

close SUM;
