#!/usr/bin/perl -w
 
 
 use strict;
 use Wurst;
 use File::Slurp;
 use File::Path;
 
 
 if ($#ARGV != 1)
 {
  print"\n#############################################################\n";
  print"#          Center for the Study of Systems Biology          #\n";
  print"#              Georgia Institute of Technology              #\n";
  print"#                                                           #\n";
  print"#                     >>> pdb2fasta <<<                     #\n";
  print"#                                                           #\n";
  print"#                    by Michal Brylinski                    #\n";
  print"#               michal.brylinski\@gatech.edu                 #\n";
  print"#                                                           #\n";
  print"#  Usage: pdb2fasta.pl <PDB file> <output FASTA>            #\n";
  print"#############################################################\n";
  die "\n\n";
 }
 
 
 my $f1 = $ARGV[0];
 my $f2 = $ARGV[1];
 
 my $tmpfXX = int(rand(100000)+100000);
 
 my $scratch = $ENV{'MYSCRATCH'}.'/pdb2fasta-'.int(rand(1000000)+1000000).int(rand(1000000)+1000000);
 
 my $tmpf01 = $scratch.'/xp2f01-'.$tmpfXX.int(rand(100001)+100001); # temporary filename
 
 #$tmpf01 = "tmpf01";
 
 mkpath($scratch);
 
 my @pdb1 = read_file( $f1 );
 
 foreach my $wpdb1 (@pdb1)
 {
  if ( length($wpdb1) > 53 )
  {
   substr($wpdb1, 21, 2) = '  ' if ( substr($wpdb1, 0, 6) eq 'ATOM  ' );
  }
 }
 
 my @fasta1 = ();
 
 my @id01 = split(/\//, $f1);
 my $id02 = pop(@id01);
 
 $id02 =~ s/.pdb//;
 
 write_file($tmpf01, "COMPND    MOL_ID: $id02;\nREMARK\n");
 
 append_file($tmpf01, @pdb1);
 
 my $seq01 = pdb_read( $tmpf01, '', '');
 
 my $seq02 = coord_get_seq($seq01);
 
 my $seq03 = uc( seq_print($seq02) );
 
 #push(@fasta1, '>'.$id02.' '.( length($seq03) - 1 )."\n");
 
 push(@fasta1, '>'.$id02."\n");
 
 push(@fasta1, $seq03);
 
 write_file($f2, @fasta1);
 
 unlink("$tmpf01") if ( -e "$tmpf01" );
 
 rmtree($scratch);
 
 exit(0);
 
 
