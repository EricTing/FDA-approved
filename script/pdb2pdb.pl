#!/usr/bin/perl -w
 
 
 use strict;
 

 if ($#ARGV != 2)
 {
  print"\n#############################################################\n";
  print"#          Center for the Study of Systems Biology          #\n";
  print"#              Georgia Institute of Technology              #\n";
  print"#                                                           #\n";
  print"#                      >>> pdb2pdb <<<                      #\n";
  print"#                                                           #\n";
  print"#                    by Michal Brylinski                    #\n";
  print"#               michal.brylinski\@gatech.edu                 #\n";
  print"#                                                           #\n";
  print"#  Usage: pdb2pdb.pl <start residue> <chain ID> <PDB file>  #\n";
  print"#############################################################\n";
  die "\n\n";
 }
 
 my $sres = $ARGV[0] - 1;
 my $fcha = $ARGV[1];
 my $fpdb = $ARGV[2];
 
 my $satom = 0;
 
 open (PDB, "$fpdb") || die "Cannot open $fpdb for reading.\n";
  my @pdb1=<PDB>;
  chomp(@pdb1);
 close (PDB);
 
 my @pdb2 = grep(/ATOM/, @pdb1);
 
 my $pat1 = '';
 
 foreach my $wpdb2 (@pdb2)
 {
  if ( length($wpdb2) > 53 )
  {
   if ( substr($wpdb2, 0, 6) eq "ATOM  " )
   {
    my $tat = substr($wpdb2, 12,  4);
    my $tam = substr($wpdb2, 17,  3);
    my $tcs = substr($wpdb2, 30, 24);
    
    if ( $tat ne $pat1 )
    {
     $pat1 = $tat;
     
     my $tch = $fcha;
     
     $tch = ' ' if ( $tch eq '_' );
     
     $sres++ if ( $tat eq " N  " );
     
     my $fac1 = 1.0;
     my $fac2 = 1.0;
     
     if ( length($wpdb2) > 65 )
     {
      $fac1 = substr($wpdb2, 54, 6) * 1.0;
      $fac2 = substr($wpdb2, 60, 6) * 1.0;
     }
     
     my $fac3 = '';
     
     if ( length($wpdb2) > 79 )
     {
      $fac3 = substr($wpdb2, 66, 14);
     }
     
     printf("ATOM%7d%5s%4s%2s%4d%28s%6.2f%6.2f%s\n", ++$satom, $tat, $tam, $tch, $sres, $tcs, $fac1, $fac2, $fac3) if ( substr($wpdb2, 13, 1) ne 'H' );
    }
   }
  }
 }
 
 print "TER\nEND\n";
 
 
