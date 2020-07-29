#!/usr/bin/env perl
#AUTHOR = Adam M. Session

open (IN, "$ARGV[0]") || die "LTRTable OutPrefix Fasta\n";

while (<IN>)
{
    chomp;
    my @fields = split /\t/;
    $P5 = "$fields[10]:$fields[0]:$fields[1]_5p";
    $P3 = "$fields[10]:$fields[0]:$fields[1]_3p";
    $name = "$fields[10]:$fields[0]:$fields[1]";
    `echo ">$P5" >> $ARGV[1].LTR.fa`;
    `blastdbcmd -dbtype nucl -db $ARGV[2] -entry $fields[10] -range $fields[3]-$fields[4] -outfmt %s >> $ARGV[1].LTR.fa`;
    `echo ">$P3" >> $ARGV[1].LTR.fa`;
    `blastdbcmd -dbtype nucl -db $ARGV[2] -entry $fields[10] -range $fields[6]-$fields[7] -outfmt %s >> $ARGV[1].LTR.fa`;
    `echo ">$name" >> $ARGV[1].inner.fa`;
    `blastdbcmd -dbtype nucl -db $ARGV[2] -entry $fields[10] -range $fields[4]-$fields[6] -outfmt %s >> $ARGV[1].inner.fa`;
}
