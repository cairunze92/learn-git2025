#perl 01Select_ATG2TGA_cds.pl total_transcripts_exon.fasta_orf-1
use strict;
my(@cds,$n,$name,%fasta);
open (IN1, "<", "Tetrahy-v1.mRNA.fa") or die "Can't open Tetrahy-v1.mRNA.fa: $!";
open (OUT1, ">", "Tetrahy-v1.mRNA.fa_ATG2TGA");
open (OUT2, ">", "Tetrahy-v1.mRNA.fa_ATG2TGA_longest");

{@cds=<IN1>;close IN1;}
$n=0;

foreach(@cds){
   if(/^>(\S+)/){
     $name=$1; $n=1;
   }elsif($n==1){
     chomp;
     $fasta{$name}.=$_;
   }
}

my %longest;
foreach my $id (keys %fasta){
   if($fasta{$id}=~/^ATG.+?TGA$/){
     my $len = length($fasta{$id});
     if(!defined $longest{$id} || $len > length($fasta{$longest{$id}})){
       $longest{$id} = $id;
     }
     print OUT1 ">$id\n$fasta{$id}\n";
   }
}

foreach my $id (keys %longest){
   print OUT2 ">$id\n$fasta{$longest{$id}}\n";
}
close OUT1;
close OUT2;
