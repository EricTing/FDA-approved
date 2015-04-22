#!/usr/bin/perl -w
 
 
 use strict;
 use File::Slurp;
 
 
 if ($#ARGV != 3)
 {
  print"\n#############################################################\n";
  print"#          Center for the Study of Systems Biology          #\n";
  print"#              Georgia Institute of Technology              #\n";
  print"#                                                           #\n";
  print"#                    >>> fasta2fasta <<<                    #\n";
  print"#                                                           #\n";
  print"#                    by Michal Brylinski                    #\n";
  print"#               michal.brylinski\@gatech.edu                 #\n";
  print"#                                                           #\n";
  print"#  Usage: fasta2fasta.pl <max characters>                   #\n";
  print"#                        <entry ID>                         #\n";
  print"#                        <input FASTA>                      #\n";
  print"#                        <output FASTA>                     #\n";
  print"#############################################################\n";
  die "\n\n";
 }
 
 
 my $f0 = $ARGV[0];
 my $f1 = $ARGV[1];
 my $f2 = $ARGV[2];
 my $f3 = $ARGV[3];
 
 my @fas1 = read_file("$f2"); chomp(@fas1);
 
 my $seq1 = '';
 
 foreach my $wfas1 (@fas1)
 {
  $seq1 .= $wfas1 if ( length($wfas1) and !( $wfas1 =~ / / ) and !( $wfas1 =~ />/ ) and !( $wfas1 =~ /:/ ) );
 }
 
 my @out1 = ();
 
 push(@out1, sprintf(">%s %d\n", $f1, length($seq1)));
 
 my @seq2 = split(//, $seq1);
 
 my $w1 = 0;
 
 foreach my $wseq2 (@seq2)
 {
  push(@out1, sprintf ("%s", $wseq2));
  
  if ( ++$w1 >= $f0 )
  {
   push(@out1, sprintf("\n"));
   
   $w1 = 0;
  }
 }
 
 push(@out1, sprintf("\n")) if ( $w1 );
 
 write_file("$f3", @out1);
 
 exit(0);
 
 
