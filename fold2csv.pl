#!/usr/bin/perl

use strict;
use warnings;

open(FOLD, "<", shift @ARGV) or die "can't open RNAfold output";

while(my $line = <FOLD>){
    if($line =~ /^>(.*?)$/){
        my @ary = split /\s+/, $1;
        my $id = $ary[0];
        my $seq = <FOLD>;
        chomp $seq;
        my $struct = <FOLD>;
        if($struct =~ /^(\S*?) \(\s?(.*?)\)$/){
            my $mfe = $2;
            my $struct = $1;
            print "$id\t$mfe\t$seq\t$struct\n";
        }
    }
}
