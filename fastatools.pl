#!/usr/bin/perl
# created by:    ejr
# creation date: 2013-06-21
# last modified: 2016-01-26 
# Provide tools for manipulation of FASTA files
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil floor);

## USAGE STATEMENT
if (!defined($ARGV[0]) or $ARGV[0] eq 'help') {
die ("
Program: fastatools (Tools for FASTA files)
Version: 2013-06-21 (EJR)

Usage:   fastatools <command> [options] fasta_in_file

Command: 
         stats       FASTA file statistics

         histogram   Print histogram of gc or lengths
                        -gc         print histogram of gc content
                        -len        print histogram of lengths

         table       Print table of sequence lengths or gc
                        -gc         print table of gc contents
                        -len        print table of lengths
                        -seq        print table of sequences

         filter      Output filtered FASTA file
                        -gc_high n  maximum %GC to retain
                        -gc_low n   minimum %GC to retain
                        -len_high n maximum sequence length to retain
                        -len_low n  minimum sequence length to retain
                        -seq_regex  regex for matching sequence
                        -name_regex regex for matching name

         nr          Combine completely identical sequences for FASTA

         upper       Convert all sequences to uppercase

         help        Print this usage statement

");
}

my $gc = 0;
my $len = 0;
my $seqtab = 0;
my $gc_low = 0;
my $gc_high = 100;
my $len_low = 0;
my $len_high = 1000000000;
my $seq_regex = '.';
my $name_regex = '.';

GetOptions (
            'gc'            => \$gc,
            'len'           => \$len,
            'seq'           => \$seqtab,
            'gc_low=i'      => \$gc_low,
            'gc_high=i'     => \$gc_high,
            'len_low=i'     => \$len_low,
            'len_high=i'    => \$len_high,
            'seq_regex=s'   => \$seq_regex,
            'name_regex=s'  => \$name_regex,
            );

# run program "mode"
my $program = shift(@ARGV);
if ($program eq "stats") {
    fasta_stats();
}elsif ($program eq "histogram") {
    fasta_histogram();
}elsif ($program eq "table") {
    fasta_table();
}elsif ($program eq "filter") {
    fasta_filter();
}elsif ($program eq "nr") {
    fasta_nr();
}elsif ($program eq "upper") {
    fasta_upper();
}else {
    die "Invalid Command\n";
}

exit(0);

## SUBROUTINES

sub fasta_stats {

my $total_at;               # total AT for all sequences
my $total_gc;               # total GC for all sequences
my $total_n;                # total N for all sequences
my $total_length;           # total length of all sequences
my @seq_lengths;            # array of sequence lengths
my $length_10_longest;      # total length of 10 longest seqs
my $length_100_longest;     # total length of 100 longest seqs
my $length_1000_longest;    # total length of 1000 longest seqs
my $length_10000_longest;   # total length of 10000 longest seqs
my $length_gt10k;           # total length of seqs longer than 10k
my $length_gt100k;          # total length of seqs longer than 100k
my $n50;                    # 50% of total length is in seqs this size or longer
my $n90;                    # 90% of total length is in seqs this size or longer

# Read in FASTA FILE
$/ = ">";                   # We do it this way to handle STDIN
<>;

while (my $line = <>) {
        chomp $line;
        my ($header, @seq) = split /\n/, $line;
        my $sequence = join '', @seq;
        my $seq_length = length($sequence);
        $total_gc += $sequence =~ tr/GCgc/GCgc/;
        $total_at += $sequence =~ tr/ATat/ATat/;
        $total_n  += $sequence =~ tr/Nn/Nn/;
        $total_length += $seq_length;

        if ( $seq_length >= 10000) { $length_gt10k  += $seq_length; }
        if ( $seq_length >= 100000){ $length_gt100k += $seq_length; }
        push @seq_lengths, $seq_length;
}

my $average_length = $total_length / scalar(@seq_lengths);
my $gc_content = 100 * $total_gc / ($total_gc + $total_at); 
my $at_content = 100 * $total_at / ($total_gc + $total_at); 
my $n_content  = 100 * $total_n  / $total_length;

my @sorted_seq_lengths = sort {$a <=> $b} @seq_lengths;

# Slice up array for calculating length of longest seqs
if (scalar(@sorted_seq_lengths) >= 10) {
    my @longest_10 = @sorted_seq_lengths[-10..-1]; 
    foreach my $i ( @longest_10) { $length_10_longest += $i; }
} else {
    $length_10_longest = "#N/A";
}

if (scalar(@sorted_seq_lengths) >= 100) {
    my @longest_100 = @sorted_seq_lengths[-100..-1];
    foreach my $i ( @longest_100) { $length_100_longest += $i; }
} else {
    $length_100_longest = "#N/A";
    
}

if (scalar(@sorted_seq_lengths) >= 1000) {
    my @longest_1000 = @sorted_seq_lengths[-1000..-1];
    foreach my $i ( @longest_1000) { $length_1000_longest += $i; }
} else {
    $length_1000_longest = "#N/A";

}
if (scalar(@sorted_seq_lengths) >= 10000) {
    my @longest_10000 = @sorted_seq_lengths[-10000..-1];
    foreach my $i ( @longest_10000) { $length_10000_longest += $i; }
} else {
    $length_10000_longest = "#N/A";

}

my $n50_length = $total_length * .5;
my $n90_length = $total_length * .9;

# Calculate N50
my $running_total = 0;
for (my $i= scalar(@sorted_seq_lengths) -1 ;$i >= 0; $i-- ) {
        $running_total += $sorted_seq_lengths[$i];
        if ($running_total >= $n50_length) { 
            $n50 = $sorted_seq_lengths[$i]; 
            last;
        } 
}

# Calculate N90
$running_total = 0;
for (my $i= scalar(@sorted_seq_lengths) -1 ;$i >= 0; $i-- ) {
        $running_total += $sorted_seq_lengths[$i];
        if ($running_total >= $n90_length) { 
            $n90 = $sorted_seq_lengths[$i]; 
            last;
        } 
}

my $longest = pop(@sorted_seq_lengths);
my $shortest = shift(@sorted_seq_lengths);
unless (defined($length_gt10k)) {$length_gt10k = "#N/A"};
unless (defined($length_gt100k)) {$length_gt100k = "#N/A"};

# Output
print "| | |\n";
print "|---|---:|\n";
print "|Number of Sequences: | ", commify(scalar(@seq_lengths)) ,"|\n";
print "|Total Length: | ", commify($total_length), "|\n";
print "|Average Length: | ", commify(int($average_length)), "|\n";
print "|Longest Sequence: | ", commify($longest) ,"|\n";
print "|Shortest Sequence: | ", commify($shortest) ,"|\n";
printf ("|%%GC: | %d%%|\n", $gc_content);
printf ("|%%N: | %d%%|\n", $n_content);
print "|N50: | ",commify($n50),"|\n";
print "|N90: | ",commify($n90),"|\n";

}

# Format a number with commas
sub commify {
        local $_  = shift;
        1 while s/^(-?\d+)(\d{3})/$1,$2/;
        return $_;
}

sub fasta_filter {

# read in fasta file
$/ = '>';
<>;

while (my $line = <>) {
    chomp $line;
    my ($header, @seq) = split /\n/, $line;
    my $sequence = join '', @seq;
    my $length = length($sequence);

    my $GC = $sequence =~ tr/GCgc/GCgc/;
    my $AT = $sequence =~ tr/ATat/ATat/;
    my $GCcontent;
    if ($GC+$AT != 0) { 
        $GCcontent = $GC / ($GC+$AT);
        $GCcontent = $GCcontent * 100;
    } else {
        $GCcontent = 0;
    }
# file can be filtered on GC and/or length
    if ($name_regex or $seq_regex) {
        if ($GCcontent <= $gc_high and $GCcontent >= $gc_low) {
            if ($length <= $len_high and $length >= $len_low) {
                if ($header =~ /$name_regex/) {
                    if ($sequence =~ /$seq_regex/) {
                        $sequence =~ s/(.{1,80})/$1\n/g;
                        print ">$header\n";
                        print "$sequence";
                    }
                }
            }
        }
    }else{
        if ($GCcontent <= $gc_high and $GCcontent >= $gc_low) {
            if ($length <= $len_high and $length >= $len_low) {
                  $sequence =~ s/(.{1,80})/$1\n/g;
                  print ">$header\n";
                  print "$sequence";
            }
        }
    }

}

}

sub fasta_histogram {

# this is a very simple histogram with no graphical output.
#
# GC HISTOGRAM
if ($gc) {
    my @GC;
    my $max = 1;
    my $min = 0;
    my $bin_size = 1;
    my %histogram;

    $/ = ">";
    <>;
    while (my $line = <>) {
        chomp $line;
        my ($header, @seq) = split /\n/, $line;
        my $sequence = join '', @seq;
        my $length = length($sequence);
        my $gc = $sequence =~ tr/GCgc/GCgc/;
        my $at = $sequence =~ tr/ATat/ATat/;
        if($gc + $at != 0) {
        my $gc_content = $gc / ($gc + $at);
        $gc_content = $gc_content * 100;
        $gc_content = floor($gc_content);
        push @GC, $gc_content;
        }

    }

    my $std = &stdev(\@GC);

    foreach my $content (@GC) {
        $histogram{ceil(($content + 1) / $bin_size) -1}++;
    }   
    print "Bin\tCount\n";

    foreach my $bin (sort { $a <=> $b } keys %histogram){
        print $bin*$bin_size , "\t$histogram{$bin}\n";
    }

}

# LENGTH HISTOGRAM
if ($len) {
    my @LEN;
    my $max = 0;
    my $min = 9999999999;
    my $bin_size = 100;
    my %histogram;

    $/ = ">";
    <>;
    while (my $line = <>) {
        chomp $line;
        my ($header, @seq) = split /\n/, $line;
        my $sequence = join '', @seq;
        my $length = length($sequence);
        push @LEN, $length;
        if ($length > $max) {$max = $length;}
        if ($length < $min) {$min = $length;}

    }

    my $std = &stdev(\@LEN);

    foreach my $content (@LEN) {
        $histogram{ceil(($content + 1) / $bin_size) -1}++;
    }  

    print "Bin\tCount\n";
    foreach my $bin (sort { $a <=> $b } keys %histogram){
        print $bin*$bin_size , "\t$histogram{$bin}\n";
    }

}

}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub fasta_table {
# output table of FASTA file lengths or GC suitable for loading into R or Excel

$/ = '>';
<>;

if ($gc) { print join("\t", 'seq.id', 'perc_gc'), "\n"; }
if ($len){ print join("\t", 'seq.id', 'length' ), "\n"; }
if ($seqtab){ print join("\t", 'seq.id', 'seq' ), "\n"; }
while (my $line = <>) {
    chomp $line;
    my ($header, @seq) = split /\n/, $line;
    my $sequence = join '', @seq;
    my $length = length($sequence);
    if ($length == 0) {$length = .001}
    my $gc_num = $sequence =~ tr/GCgc/GCgc/;
    my $gc_content = $gc_num / $length;

    print "$header";
    if ($gc) {
        printf ("\t%.02f",$gc_content);
        
    }
    if ($len) {
        print "\t$length";
    }
    if ($seqtab) {
        print "\t$sequence";
    }
    print "\n";

}
}

sub fasta_nr {

my %fasta;

$/=">";
<>;
while (my $line = <>) {
    chomp $line;
    my ($header, @seq) = split /\n/, $line;
    $header =~ s/^(.+)\s*/$1/;
    my $sequence = join '', @seq;
    $sequence =~ s/\*//g;
    push @{$fasta{$sequence}}, $header;
}

$/="\n";

foreach my $seq (keys %fasta) {
    print ">";
    my $header =  join " ", @{$fasta{$seq}};
    $header =~ s/\:/ /g;
    print $header, "\n";;
    $seq =~ s/(.{1,60})/$1\n/g;
    print $seq;
}

}

sub fasta_upper {

my %fasta;
my $header;
while (my $line = <>) {
    chomp $line;
    if ($line =~ />/) {
        $header = $line;
    }else{
        $fasta{$header} .= $line;
    }
}

foreach my $header (keys %fasta) {
    print ">", $header, "\n";
    $fasta{$header} = uc($fasta{$header});
    $fasta{$header} =~ s/(.{1,60})/$1\n/g;
    print $fasta{$header};
}
}
