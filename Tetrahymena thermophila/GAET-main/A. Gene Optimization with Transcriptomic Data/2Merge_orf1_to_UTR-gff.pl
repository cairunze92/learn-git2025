# perl 2Merge_orf1_to_UTR-gff.pl stringtie.gtf transcripts.fa_ATG2TGA_longest
use strict;
use warnings;

my(@GFF, @list, %fa_seq);
my $current_id = "";

open (IN1, "<", "$ARGV[0]") or die "Can't open GTF file $ARGV[0] : $!";
open (IN2, "<", "$ARGV[1]") or die "Can't open FASTA file $ARGV[1] : $!";
open (OUT1, ">", "$ARGV[0]_UTR.gff3");
open (OUT2, ">", "$ARGV[0]_Not_in_orffinder_mRNA.txt");
{
    @GFF = <IN1>;  close IN1;
    @list = <IN2>; close IN2;
}

# 解析FASTA文件，构建序列哈希表
print STDERR "解析FASTA文件...\n";
foreach my $line (@list) {
    chomp $line;
    if ($line =~ /^>(\S+)/) {
        $current_id = $1;
        $fa_seq{$current_id} = "";
    } elsif ($current_id ne "") {
        $fa_seq{$current_id} .= $line;
    }
}

print STDERR "找到 " . scalar(keys %fa_seq) . " 个序列\n";

# 识别每个序列中的ORF
my(%list2, %gene, $ori, $n, $mRNA, %mRNA, $ORF);
my $orf_count = 0;

foreach my $seq_id (keys %fa_seq) {
    my $seq = $fa_seq{$seq_id};
    my ($start_pos, $end_pos, $orf_seq) = find_longest_orf($seq);
    
    if ($start_pos > 0 && $end_pos > $start_pos) {
        $orf_count++;
        $ORF = "ORF1"; # 使用简单的ORF命名
        $ori = "+"; # 假设正向
        $list2{$seq_id} = "$start_pos\t$end_pos\t$ori\t$ORF";
        print STDERR "找到ORF: $seq_id\t$ORF\t$start_pos\t$end_pos\t" . ($end_pos - $start_pos + 1) . "\n" if $orf_count < 5;
    }
}

print STDERR "总共找到 $orf_count 个ORF\n";
$n = keys %list2;
print STDERR "有 $n 个转录本包含ORF\n";

# 处理GTF文件中的转录本和外显子信息
my($transcript_id, $transcrpt_ori, $orf_ori);
$transcript_id = "";

foreach my $line (@GFF) {
    if ($line =~ /transcript.+?transcript_id\s+"(.+?)";.+?FPKM "(.+?)";/) {
        $transcript_id = $1;
    } elsif ($transcript_id ne "" && $line =~ /exon.+?transcript_id\s+"($transcript_id)";/) {
        $mRNA{$transcript_id} .= $line;
    }
}

$n = keys %mRNA;
print STDERR "找到 $n 个转录本有外显子信息\n";

# 整合ORF和转录本信息，生成GFF3
my($chr, $gene_start, $gene_end, $mRNA_id, $gene_id, @exon, @exon_coord, $mark, @CDS_coord, $L1, $L2, $i, @n, @m, $start, $end, $length);

foreach my $line (@GFF) {
    if ($line =~ /(chr\w+|scaffold\w+|STRG\.\d+)\s+\S+\s+transcript\s+(\d+)\s+(\d+)\s+\S+\s+(.)\s+.*?gene_id\s+"(.+?)";\s+transcript_id\s+"(.+?)";/) {
        $chr = $1;
        $gene_start = $2;
        $gene_end = $3;
        $mRNA_id = $6;
        $gene_id = $5;
        $transcrpt_ori = $4;
        
        if (exists $list2{$mRNA_id}) {
            print STDERR "处理转录本: $mRNA_id\n" if $. < 5;
            @exon = split /\n/, $mRNA{$mRNA_id};
            undef @exon_coord;
            undef @CDS_coord;
            $length = 0;
            
            foreach my $exon_line (@exon) {
                if ($exon_line =~ /\S+\s+\S+\s+exon\s+(\d+)\s+(\d+)\s+.*?transcript_id\s+"($mRNA_id)";/) {
                    push(@exon_coord, ($1, $2));
                    $length += $2 - $1 + 1;
                    push(@CDS_coord, ($1, $2));
                }
            }
            
            if ($list2{$mRNA_id} =~ /(\d+)\t(\d+)\t(\+|-)\t(.+)/) {
                $L1 = $1;
                $L2 = $2;
                $orf_ori = $3;
                $ORF = $4;
                
                # 确定最终方向
                if ($transcrpt_ori eq "+" && $orf_ori eq "+") {
                    $ori = "+";
                } elsif ($transcrpt_ori eq "+" && $orf_ori eq "-") {
                    $ori = "-";
                } elsif ($transcrpt_ori eq "-" && $orf_ori eq "+") {
                    $ori = "-";
                    my $a = $L1;
                    my $b = $L2;
                    $L1 = $length - $b - 1;
                    $L2 = $length - $a - 1;
                } elsif ($transcrpt_ori eq "-" && $orf_ori eq "-") {
                    $ori = "+";
                    my $L3 = $L2;
                    $L2 = $length - $L1 - 1;
                    $L1 = $length - $L3 - 1;
                }
                
                if ($transcrpt_ori eq ".") {
                    $ori = $orf_ori;
                    if ($ori eq "-") {
                        my $a = $L1;
                        my $b = $L2;
                        $L1 = $length - $b - 1;
                        $L2 = $length - $a - 1;
                    }
                }
                
                $i = 0;
                undef @n;
                undef @m;
                $mark = 0;
                my @temp_CDS = @CDS_coord;
                @CDS_coord = @temp_CDS;
                
                while ($i < @exon_coord) {
                    $start = $exon_coord[$i] + $L1;
                    
                    if ($i+1 < @exon_coord && $exon_coord[$i+1] > $start) {
                        if ($L1 > 0) {
                            push @m, $exon_coord[$i], $start - 1;
                            shift @CDS_coord if @CDS_coord;
                            unshift @CDS_coord, $start;
                            $L1 = 0;
                        }
                    } else {
                        push @m, $exon_coord[$i], $exon_coord[$i+1] if $i+1 < @exon_coord;
                        shift @CDS_coord if @CDS_coord;
                        shift @CDS_coord if @CDS_coord;
                        $L1 = $L1 - ($exon_coord[$i+1] - $exon_coord[$i]) - 1 if $i+1 < @exon_coord;
                    }
                    
                    $end = $exon_coord[$i] + $L2;
                    
                    if ($i+1 < @exon_coord && $exon_coord[$i+1] > $end) {
                        if ($L2 > 0) {
                            push @n, $end + 1, $exon_coord[$i+1];
                            pop @CDS_coord if @CDS_coord;
                            push @CDS_coord, $end;
                            $L2 = 0;
                        } elsif ($L2 <= 0) {
                            push @n, $exon_coord[$i], $exon_coord[$i+1];
                            pop @CDS_coord if @CDS_coord;
                            pop @CDS_coord if @CDS_coord;
                        }
                    } elsif ($i+1 < @exon_coord && $exon_coord[$i+1] == $end) {
                        $L2 = 0;
                        $mark = 1;
                    } elsif ($i+1 < @exon_coord) {
                        my $len = $exon_coord[$i+1] - $exon_coord[$i] + 1;
                        $L2 = $L2 - ($exon_coord[$i+1] - $exon_coord[$i]) - 1;
                    }
                    
                    $i += 2;
                }
                
                if ($mark == 0 && @n > 0) {
                    pop @CDS_coord if @CDS_coord;
                    push @CDS_coord, $n[0] - 1;
                }
                
                # 输出GFF3
                if ($ori eq "+") {
                    print OUT1 "$chr\tAUGUSTUS\tgene\t$gene_start\t$gene_end\t.\t$ori\t.\tID=${mRNA_id}_$ORF;Name=${mRNA_id}_$ORF;Note=\n";
                    print OUT1 "$chr\tAUGUSTUS\tmRNA\t$gene_start\t$gene_end\t.\t$ori\t.\tID=${mRNA_id}_$ORF;Parent=${mRNA_id}_$ORF\n";
                    
                    $i = 0;
                    while ($i < @m) {
                        print OUT1 "$chr\tAUGUSTUS\tfive_prime_UTR\t$m[$i]\t$m[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.utr5;Parent=${mRNA_id}_$ORF\n";
                        $i += 2;
                    }
                    
                    $i = 0;
                    while ($i < @CDS_coord) {
                        print OUT1 "$chr\tAUGUSTUS\tCDS\t$CDS_coord[$i]\t$CDS_coord[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.cds;Parent=${mRNA_id}_$ORF\n";
                        $i += 2;
                    }
                    
                    $i = 0;
                    while ($i < @n) {
                        print OUT1 "$chr\tAUGUSTUS\tthree_prime_UTR\t$n[$i]\t$n[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.utr3;Parent=${mRNA_id}_$ORF\n";
                        $i += 2;
                    }
                } elsif ($ori eq "-") {
                    print OUT1 "$chr\tAUGUSTUS\tgene\t$gene_start\t$gene_end\t.\t$ori\t.\tID=${mRNA_id}_$ORF;Name=${mRNA_id}_$ORF;Note=\n";
                    print OUT1 "$chr\tAUGUSTUS\tmRNA\t$gene_start\t$gene_end\t.\t$ori\t.\tID=${mRNA_id}_$ORF;Parent=${mRNA_id}_$ORF\n";
                    
                    $i = 0;
                    while ($i < @m) {
                        print OUT1 "$chr\tAUGUSTUS\tthree_prime_UTR\t$m[$i]\t$m[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.utr3;Parent=${mRNA_id}_$ORF\n";
                        $i += 2;
                    }
                    
                    $i = 0;
                    while ($i < @CDS_coord) {
                        print OUT1 "$chr\tAUGUSTUS\tCDS\t$CDS_coord[$i]\t$CDS_coord[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.cds;Parent=${mRNA_id}_$ORF\n";
                        $i += 2;
                    }
                    
                    $i = 0;
                    while ($i < @n) {
                        print OUT1 "$chr\tAUGUSTUS\tfive_prime_UTR\t$n[$i]\t$n[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.utr5;Parent=${mRNA_id}_$ORF\n";
                        $i += 2;
                    }
                }
            }
        } else {
            print OUT2 "$mRNA_id\n";
        }
    }
}

# 找出最长ORF的辅助函数
sub find_longest_orf {
    my ($seq) = @_;
    my $seq_length = length($seq);
    my @start_codons = ("ATG");
    my @stop_codons = ("TAA", "TAG", "TGA");
    
    my $longest_start = -1;
    my $longest_end = -1;
    my $longest_length = 0;
    my $longest_seq = "";
    
    # 对三个读码框分别搜索
    for (my $frame = 0; $frame < 3; $frame++) {
        my $pos = $frame;
        while ($pos <= $seq_length - 3) {
            my $codon = substr($seq, $pos, 3);
            
            if (grep {$_ eq $codon} @start_codons) {
                my $orf_start = $pos;
                my $orf_pos = $pos + 3;
                my $found_stop = 0;
                
                while ($orf_pos <= $seq_length - 3) {
                    my $current_codon = substr($seq, $orf_pos, 3);
                    
                    if (grep {$_ eq $current_codon} @stop_codons) {
                        my $orf_end = $orf_pos + 2;
                        my $orf_length = $orf_end - $orf_start + 1;
                        
                        if ($orf_length > $longest_length) {
                            $longest_length = $orf_length;
                            $longest_start = $orf_start + 1; # 转为1-索引
                            $longest_end = $orf_end + 1;     # 转为1-索引
                            $longest_seq = substr($seq, $orf_start, $orf_length);
                        }
                        
                        $found_stop = 1;
                        last;
                    }
                    
                    $orf_pos += 3;
                }
                
                # 如果没有找到终止密码子，跳到下一个起始密码子
                if (!$found_stop) {
                    $pos += 3;
                    next;
                }
            }
            
            $pos += 3;
        }
    }
    
    return ($longest_start, $longest_end, $longest_seq);
}
