#!/usr/bin/perl

=head1 usage

detect miRNA target site in Input sequence

usage:
	perl mirna-target-scan.pl -m=mrna.fa -s=mirna.fa

=cut

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

## global variable=======================================
my($mirna_file, $mrna_file);

my %mirna_record;

my %target_8mer;
my %target_7mera1;
my %target_7merm8;

my %target_6mer;
my %target_os6mer;

## main================================================
GetOptions("m=s" => \$mrna_file,
	"s=s" => \$mirna_file) or die "Commandline error.";

parse_mirna($mirna_file);

parse_mrna($mrna_file);
## functions================================================
sub parse_mirna{
	my($input_file) = @_;
	
	my $mirna_fh = Bio::SeqIO->new(-file => $input_file, -format => "fasta");
	
	while(my $mirna = $mirna_fh->next_seq){
		#print $mirna->id, revcomp($mirna->seq), "\n";
		my $tmseq = $mirna->seq;
        $tmseq =~ tr/T/U/;
		$mirna->seq($tmseq);
		$mirna_record{$mirna->id} = $mirna->seq;
		
		# canonical target site
		$target_8mer{$mirna->id} = revcomp(substr($mirna->seq, 1, 7)) . "A";
		$target_7mera1{$mirna->id} = revcomp(substr($mirna->seq, 1, 6)) . "A";
		$target_7merm8{$mirna->id} = revcomp(substr($mirna->seq, 1, 7));
		
		# marginal target site
		$target_6mer{$mirna->id} = revcomp(substr($mirna->seq, 1, 6));
		$target_os6mer{$mirna->id} = revcomp(substr($mirna->seq, 2, 6));		
	}
}

sub revcomp{
	my($instr) = @_;
	
	my $outstr = reverse $instr;
	$outstr =~ tr/UACG/AUGC/;
	
	return $outstr;
}

sub parse_mrna{
	my($input_file) = @_;
	
	my $mrna_fh = Bio::SeqIO->new(-file => $input_file, -format => "fasta");
	
	while(my $mrna = $mrna_fh->next_seq){
		scan_target($mrna);				
	}
}

sub scan_target{
	my($mrna) = @_;
	
	foreach my $mirna(sort keys %mirna_record){
		my $tmseq = $mrna->seq;
        $tmseq =~ tr/T/U/;
		
		my($regex, $rep, $site_type);
		
		# find 8mer site
		$regex = $target_8mer{$mirna};
		$site_type = "8mer";
		$rep = 'N' x length($regex);
		
		while($tmseq =~ /$regex/g){			
			#print join("\t", $mirna, $mrna->id, $-[0] + 1, $+[0], $site_type, $regex, $mirna_record{$mirna}, $mrna->seq), "\n";
			print join("\t", $mirna, $mrna->id, $-[0] + 1, $+[0], $site_type), "\n";
		}
		$tmseq =~ s/$regex/$rep/g;
		
		# find 7mera1 site
		$regex = $target_7mera1{$mirna};
		$site_type = "7mera1";
		$rep = 'N' x length($regex);
		
		while($tmseq =~ /$regex/g){			
			print join("\t", $mirna, $mrna->id, $-[0] + 1, $+[0], $site_type), "\n";
		}
		$tmseq =~ s/$regex/$rep/g;
		
		# find 7merm8 site
		$regex = $target_7merm8{$mirna};
		$site_type = "7merm8";
		$rep = 'N' x length($regex);
		
		while($tmseq =~ /$regex/g){			
			print join("\t", $mirna, $mrna->id, $-[0] + 1, $+[0], $site_type), "\n";
		}
		$tmseq =~ s/$regex/$rep/g;
		
		# find 6mer site
		$regex = $target_6mer{$mirna};
		$site_type = "6mer";
		$rep = 'N' x length($regex);
		
		while($tmseq =~ /$regex/g){			
			print join("\t", $mirna, $mrna->id, $-[0] + 1, $+[0], $site_type), "\n";
		}
		$tmseq =~ s/$regex/$rep/g;
		
		# find offset 6mer site
		$regex = $target_os6mer{$mirna};
		$site_type = "os6mer";
		$rep = 'N' x length($regex);
		
		while($tmseq =~ /$regex/g){			
			print join("\t", $mirna, $mrna->id, $-[0] + 1, $+[0], $site_type), "\n";
		}
		$tmseq =~ s/$regex/$rep/g;
	}
	return 0;
}
		
