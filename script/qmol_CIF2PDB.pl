#!/home/michal/local/bin/perl -w
 
 
 use strict;
 use File::Copy;
 use File::Path;
 use File::Slurp;
 use Compress::Zlib;
 use Chemistry::Mol;
 use Chemistry::File::PDB;
 use Chemistry::File::SDF;
 use Chemistry::Canonicalize ':all';
 use Chemistry::Bond::Find ':all';
 use Cwd;
 
 
 sub three2one {
  
  my $resX3 = shift(@_);
  my $resX1 = "";
  
     if ( $resX3 eq 'ALA' ) { $resX1 = 'A'; }
  elsif ( $resX3 eq 'ASX' ) { $resX1 = 'B'; }
  elsif ( $resX3 eq 'CYS' ) { $resX1 = 'C'; }
  elsif ( $resX3 eq 'ASP' ) { $resX1 = 'D'; }
  elsif ( $resX3 eq 'GLU' ) { $resX1 = 'E'; }
  elsif ( $resX3 eq 'PHE' ) { $resX1 = 'F'; }
  elsif ( $resX3 eq 'GLY' ) { $resX1 = 'G'; }
  elsif ( $resX3 eq 'HIS' ) { $resX1 = 'H'; }
  elsif ( $resX3 eq 'ILE' ) { $resX1 = 'I'; }
  elsif ( $resX3 eq 'LYS' ) { $resX1 = 'K'; }
  elsif ( $resX3 eq 'LEU' ) { $resX1 = 'L'; }
  elsif ( $resX3 eq 'MET' ) { $resX1 = 'M'; }
  elsif ( $resX3 eq 'ASN' ) { $resX1 = 'N'; }
  elsif ( $resX3 eq 'PRO' ) { $resX1 = 'P'; }
  elsif ( $resX3 eq 'GLN' ) { $resX1 = 'Q'; }
  elsif ( $resX3 eq 'ARG' ) { $resX1 = 'R'; }
  elsif ( $resX3 eq 'SER' ) { $resX1 = 'S'; }
  elsif ( $resX3 eq 'THR' ) { $resX1 = 'T'; }
  elsif ( $resX3 eq 'SEC' ) { $resX1 = 'U'; }
  elsif ( $resX3 eq 'VAL' ) { $resX1 = 'V'; }
  elsif ( $resX3 eq 'TRP' ) { $resX1 = 'W'; }
  elsif ( $resX3 eq 'XAA' ) { $resX1 = 'X'; }
  elsif ( $resX3 eq 'TYR' ) { $resX1 = 'Y'; }
  elsif ( $resX3 eq 'GLX' ) { $resX1 = 'Z'; }
  else { die "Unknown residue passed to three2one: $resX3   (in $ARGV[0])\n"; }
  
  return $resX1;
 }
 
 
 sub one2three {
  
  my $resX1 = shift(@_);
  my $resX3 = "";
  
     if ( $resX1 eq 'A' ) { $resX3 = 'ALA'; }
  elsif ( $resX1 eq 'B' ) { $resX3 = 'ASX'; }
  elsif ( $resX1 eq 'C' ) { $resX3 = 'CYS'; }
  elsif ( $resX1 eq 'D' ) { $resX3 = 'ASP'; }
  elsif ( $resX1 eq 'E' ) { $resX3 = 'GLU'; }
  elsif ( $resX1 eq 'F' ) { $resX3 = 'PHE'; }
  elsif ( $resX1 eq 'G' ) { $resX3 = 'GLY'; }
  elsif ( $resX1 eq 'H' ) { $resX3 = 'HIS'; }
  elsif ( $resX1 eq 'I' ) { $resX3 = 'ILE'; }
  elsif ( $resX1 eq 'K' ) { $resX3 = 'LYS'; }
  elsif ( $resX1 eq 'L' ) { $resX3 = 'LEU'; }
  elsif ( $resX1 eq 'M' ) { $resX3 = 'MET'; }
  elsif ( $resX1 eq 'N' ) { $resX3 = 'ASN'; }
  elsif ( $resX1 eq 'P' ) { $resX3 = 'PRO'; }
  elsif ( $resX1 eq 'Q' ) { $resX3 = 'GLN'; }
  elsif ( $resX1 eq 'R' ) { $resX3 = 'ARG'; }
  elsif ( $resX1 eq 'S' ) { $resX3 = 'SER'; }
  elsif ( $resX1 eq 'T' ) { $resX3 = 'THR'; }
  elsif ( $resX1 eq 'U' ) { $resX3 = 'SEC'; }
  elsif ( $resX1 eq 'V' ) { $resX3 = 'VAL'; }
  elsif ( $resX1 eq 'W' ) { $resX3 = 'TRP'; }
  elsif ( $resX1 eq 'X' ) { $resX3 = 'XAA'; }
  elsif ( $resX1 eq 'Y' ) { $resX3 = 'TYR'; }
  elsif ( $resX1 eq 'Z' ) { $resX3 = 'GLX'; }
  else { die "Unknown residue passed to one2three: $resX1   (in $ARGV[0])\n"; }
  
  return $resX3;
 }
 
 
 if ($#ARGV != 1)
 {
  print"\n#############################################################\n";
  print"#          Center for the Study of Systems Biology          #\n";
  print"#              Georgia Institute of Technology              #\n";
  print"#                                                           #\n";
  print"#                   >>> qmol_CIF2PDB <<<                    #\n";
  print"#                                                           #\n";
  print"#                    by Michal Brylinski                    #\n";
  print"#               michal.brylinski\@gatech.edu                 #\n";
  print"#                                                           #\n";
  print"#  Usage: qmol_CIF2PDB.pl <CIF file, gzipped> <PDB-ID>      #\n";
  print"#############################################################\n";
  die "\n\n";
 }
 
 
 # PARAMETERS --------------------------------------------------------------------------------------
 
 my $fcif1 = $ARGV[0];
 my $fpdb1 = $ARGV[1];
 
 my $jctrip = '/home/michal/local/src/jackal_64bit/bin/ctrip';
 my $lpc    = '/home/michal/local/src/LPC/lpcEx';
 
 die "Cannot find $fcif1\n"  if ( !( -e $fcif1 ) );
 die "Cannot find $jctrip\n" if ( !( -e $jctrip ) );
 die "Cannot find $lpc\n"    if ( !( -e $lpc ) );
 
 my $minat =   6;
 my $maxat = 150;
 my $minbr =   5;
 my $minmt =   3;
 my $minfs =   3;
 
 my $tmpfXX = int(rand(100000)+100000);
 
 my $scratch = 'qmol_cif2pdb-'.int(rand(1000000)+1000000).int(rand(1000000)+1000000);
 
 my $ttmpf01 = 'xcif01-'.$tmpfXX.int(rand(100001)+100001);         # temporary filename
 my $ttmpf02 = 'xcif02-'.$tmpfXX.int(rand(200002)+200002);         # temporary filename
 my $ttmpf03 = 'xcif03-'.$tmpfXX.int(rand(300003)+300003);         # temporary filename
 my $ttmpf04 = 'xcif04-'.$tmpfXX.int(rand(400004)+400004);         # temporary filename
 my $ttmpf05 = 'xcif05-'.$tmpfXX.int(rand(500005)+500005);         # temporary filename
 
 #$ttmpf01 = 'xcif01';
 #$ttmpf02 = 'xcif02';
 #$ttmpf03 = 'xcif03';
 #$ttmpf04 = 'xcif04';
 #$ttmpf05 = 'xcif05';
 
 mkpath($scratch);
 
 my $cwd01 = getcwd();
 
 chdir($scratch) or die "Cannot chdir to $scratch $!";
 
 my @files1 = ();
 
 
 # PARSE CIF ---------------------------------------------------------------------------------------
 
 my $w1 = 0;
 my $pat1 = '';
 my @mol1 = ();
 my %mod1 = ();
 my @tt1 = ();
 
 my $gz = gzopen("$cwd01/$fcif1", "rb") or die "Cannot open $cwd01/$fcif1: $gzerrno\n";
 
 while ( $gz->gzreadline($_) > 0)
 {
  my($wcif1) = $_;
  chomp($wcif1);
  
  $w1 = 0 if ( $wcif1 eq '# ' );
  
  if ( $w1 )
  {
   my $wcif2 = $wcif1;
   
   while ( $wcif2 =~ m/\t/ )
   {
    $wcif2 =~ s/\t/ /g;
   }
   
   while ( $wcif2 =~ m/  / )
   {
    $wcif2 =~ s/  / /g;
   }
   
   my @tt2 = split(/ /, $wcif2);
   
   if ( $tt2[6].':'. $tt2[7] ne $pat1 )
   {
    if ( length($pat1) )
    {
     my $pat2 = '';
     
     foreach my $wtt1 (@tt1)
     {
      my @tt4 = split(/ /, $wtt1);
      
      my $tt5 = pop(@tt4);
      
      $pat2 .= $wtt1.'&' if ( $tt5 < 2 );
     }
     
     substr($pat2, -1, 1) = '' if ( substr($pat2, -1, 1) eq '&' );
     
     push(@mol1, $pat2);
    }
    
    @tt1 = ();
    
    $pat1 = $tt2[6].':'. $tt2[7];
   }
   
   push(@tt1, $wcif2);
  }
  
  $w1 = 1 if ( $wcif1 eq '_atom_site.pdbx_PDB_model_num ' );
  
  if ( $wcif1 =~ m/modres/ )
  {
   if ( length($wcif1) > 14 )
   {
    my $wcif3 = $wcif1;
    
    while ( $wcif3 =~ m/  / )
    {
     $wcif3 =~ s/  / /g;
    }
    
    my @tt3 = split(/ /, $wcif3);
    my $ntt3 = @tt3;
    
    if ( $ntt3 > 8 )
    {
     $mod1{$tt3[3].':'.$tt3[4].':'.$tt3[5]} = $tt3[9] if ( $tt3[0] =~ m/modres/ and $tt3[1] =~ m/modres/ );
    }
   }
  }
 }
 
 $gz->gzclose();
 
 my $pat5 = '';
 
 foreach my $wtt1 (@tt1)
 {
  my @tt4 = split(/ /, $wtt1);
  
  my $tt5 = pop(@tt4);
  
  $pat5 .= $wtt1.'&' if ( $tt5 < 2 );
 }
 
 substr($pat5, -1, 1) = '' if ( substr($pat5, -1, 1) eq '&' );
 
 push(@mol1, $pat5);
 
 # split data
 
 my @pro1 = ();
 my @lig1 = ();
 my @wat1 = ();
 
 foreach my $wmol1 (@mol1)
 {
  my @mol2 = split(/\&/, $wmol1);
  
  my $n1 = 0;
  my $n2 = 0;
  my $n3 = 0;
  my $n4 = 0;
  my $n5 = 0;
  my $n6 = 0;
  my $n7 = 0;
  my $n8 = 0;
  
  foreach my $wmol2 (@mol2)
  {
   my @mol3 = split(/\ /, $wmol2);
   
   $n1++ if ( $mol3[0] eq 'ATOM' );
   $n2++ if ( $mol3[3] eq 'CA' );
   $n3++ if ( $mol3[3] eq 'N' );
   $n4++ if ( $mol3[3] eq 'C' );
   $n5++ if ( $mol3[3] eq 'O' );
   $n6++ if ( $mol3[3] eq 'CB' );
   $n7++ if ( $mol3[5] eq 'ALA' or 
              $mol3[5] eq 'CYS' or 
              $mol3[5] eq 'ASP' or 
              $mol3[5] eq 'GLU' or 
              $mol3[5] eq 'PHE' or 
              $mol3[5] eq 'GLY' or 
              $mol3[5] eq 'HIS' or 
              $mol3[5] eq 'ILE' or 
              $mol3[5] eq 'LYS' or 
              $mol3[5] eq 'LEU' or 
              $mol3[5] eq 'MET' or 
              $mol3[5] eq 'ASN' or 
              $mol3[5] eq 'PRO' or 
              $mol3[5] eq 'GLN' or 
              $mol3[5] eq 'ARG' or 
              $mol3[5] eq 'SER' or 
              $mol3[5] eq 'THR' or 
              $mol3[5] eq 'VAL' or 
              $mol3[5] eq 'TRP' or 
              $mol3[5] eq 'TYR' );
   $n8++ if ( $mol3[5] eq 'WAT' or $mol3[5] eq 'HOH' );
  }
  
  if ( $n1 >= 20 and $n2 >= 30 and $n3 >= 20 and $n4 >= 20 and $n5 >= 20 and $n6 >= 20 and $n7 >= 20 )
  {
   push(@pro1, $wmol1);
  }
  elsif ( $n8 >= 5 )
  {
   push(@wat1, $wmol1);
  }
  else
  {
   push(@lig1, $wmol1);
  }
 }
 
 
 # proteins
 
 my @cha1 = ();
 my %cha2 = ();
 my %cha3 = ();
 
 foreach my $wpro1 (@pro1)
 {
  my $chain = '';
  
  my @pro2 = split(/\&/, $wpro1);
  
  open (PDB, ">$ttmpf01.pdb") || die "Cannot open $ttmpf01.pdb for writing\n";
  
  foreach my $wpro2 (@pro2)
  {
   my @pro3 = split(/\ /, $wpro2);
   my $npro3 = @pro3;
   
   if ( $npro3 > 11 )
   {
    $pro3[5] = $mod1{$pro3[6].':'.$pro3[5].':'.$pro3[8]} if ( exists $mod1{$pro3[6].':'.$pro3[5].':'.$pro3[8]} );
    $pro3[0] = 'ATOM' if ( $pro3[0] eq 'HETATM' );
    
    if ( $pro3[2] ne 'H' )
    {
        if ( length($pro3[3]) == 1 )
     {
      $pro3[3] .= '  ';
     }
     elsif ( length($pro3[3]) == 2 )
     {
      $pro3[3] .= ' ';
     }
     elsif ( length($pro3[3]) > 3 )
     {
      $pro3[3] = substr($pro3[3], 0, 3);
     }
     
     printf PDB ("ATOM%7d%5s%4s%2s%4d%12.3f%8.3f%8.3f\n", $pro3[1], $pro3[3], $pro3[5], $pro3[6], $pro3[8], $pro3[10], $pro3[11], $pro3[12]);
     
     $chain = $pro3[6];
    }
   }
  }
  
  print PDB "END\n";
  
  close (PDB);
  
  if ( !grep(/$chain/, @cha1) )
  {
   `$jctrip -prm 2 $ttmpf01.pdb > /dev/null 2>&1`;
   
   `cp $ttmpf01.pdb $ttmpf01\_fix.pdb`;
   
   my @pro4 = read_file("$ttmpf01\_fix.pdb");
   
   my @pro5 = grep(/ATOM/, @pro4);
   
   my $satom = 0;
   my $sres  = 0;
   
   if ( length($chain) == 1 )
   {
    open (PDB, ">$fpdb1$chain.pdb") || die "Cannot open $fpdb1$chain.pdb for writing\n";
    
    push(@cha1, "$chain");
    $cha2{$chain} = 0;
    $cha3{$chain} = '';
    
    foreach my $wpro4 (@pro4)
    {
     if ( length($wpro4) > 53 )
     {
      if ( substr($wpro4, 0, 6) eq "ATOM  " )
      {
       my $tat = substr($wpro4, 12,  4);
       my $tam = substr($wpro4, 17,  3);
       my $tcs = substr($wpro4, 30, 24);
       
       $sres++ if ( $tat eq " N  " );
       
       printf PDB ("ATOM%7d%5s%4s%6d%28s\n", ++$satom, $tat, $tam, $sres, $tcs);
      }
     }
    }
    
    close (PDB);
    
    push(@files1, "$fpdb1$chain.pdb")
   }
  }
  
  unlink("$ttmpf01.pdb") if ( -e "$ttmpf01.pdb" );
  unlink("$ttmpf01\_fix.pdb") if ( -e "$ttmpf01\_fix.pdb" );
 }
 
 
 # remove redundant conformations
 
 foreach my $wlig1 (@lig1)
 {
  my @lig3 = split(/\&/, $wlig1);
  
  my @alt1 = ();
  
  my $alt2 = 0;
  
  foreach my $wlig3 (@lig3)
  {
   my @alt3 = split(/\ /, $wlig3);
   
   my $alt4 = 0;
   
   for ( my $xh = 0; $xh < $alt2; $xh++ )
   {
    if ( $alt1[$xh][0] eq $alt3[4] )
    {
     $alt1[$xh][1]++;
     
     $alt4 = 1;
     
     last;
    }
   }
   
   if ( !$alt4 )
   {
    $alt1[$alt2][0] = $alt3[4];
    $alt1[$alt2][1] = 1;
    
    $alt2++;
   }
  }
  
  my $alt5 = '';
  my $alt6 = 0;
  
  for ( my $xh = 0; $xh < $alt2; $xh++ )
  {
   if ( $alt1[$xh][1] > $alt6 )
   {
    $alt5 = $alt1[$xh][0];
    
    $alt6 = $alt1[$xh][1];
   }
  }
  
  my $alt7 = '';
  
  foreach my $wlig3 (@lig3)
  {
   my @alt8 = split(/\ /, $wlig3);
   
   $alt7 .= $wlig3.'&' if ( $alt8[4] eq $alt5 );
  }
  
  if ( length($alt7) )
  {
   substr($alt7, -1, 1) = '' if ( substr($alt7, -1, 1) eq '&' );
  }
  
  $wlig1 = $alt7;
 }
 
 
 # refine ligands
 
 foreach my $wlig1 (@lig1)
 {
  my @lig3 = split(/\&/, $wlig1);
  
  my @lig4 = ();
  
  foreach my $wlig3 (@lig3)
  {
   my @lig3a = split(/\ /, $wlig3);
   
   my $w3 = 1;
   
   foreach my $wlig4 (@lig4)
   {
    my @lig4a = split(/\ /, $wlig4);
    
    $w3 = 0 if ( sqrt(($lig3a[10]-$lig4a[10])**2 + ($lig3a[11]-$lig4a[11])**2 + ($lig3a[12]-$lig4a[12])**2) < 1.0 );
    
    last if ( !$w3 );
   }
   
   push(@lig4, $wlig3) if ( $w3 );
  }
  
  my $wlig1 = '';
  
  foreach my $wlig4 (@lig4)
  {
   $wlig1 .= $wlig4.'&';
  }
 }
 
 
 # search for metals
 
 my @met1 = ();
 
 foreach my $wlig1 (@lig1)
 {
  my @lig3 = split(/\&/, $wlig1);
  
  foreach my $wlig3 (@lig3)
  {
   my @lig4 = split(/\ /, $wlig3);
   my $nlig4 = @lig4;
   
   if ( $lig4[0] eq 'HETATM' )
   {
    if ( $lig4[2] eq 'ZN' or 
         $lig4[2] eq 'CO' or 
         $lig4[2] eq 'NI' or 
         $lig4[2] eq 'FE' or 
         $lig4[2] eq 'CU' or 
         $lig4[2] eq 'MN' or 
         $lig4[2] eq 'MG' or 
         $lig4[2] eq 'CA' )
    {
     push(@met1, "$lig4[2]:$lig4[10]:$lig4[11]:$lig4[12]");
    }
   }
  }
 }
 
 
 # search for Fe-S clusters
 
 my @fes1 = ();
 
 foreach my $wlig1 (@lig1)
 {
  my @lig3 = split(/\&/, $wlig1);
  
  foreach my $wlig3 (@lig3)
  {
   my @lig4 = split(/\ /, $wlig3);
   my $nlig4 = @lig4;
   
   if ( $lig4[0] eq 'HETATM' )
   {
    push(@fes1, "$lig4[2]:$lig4[10]:$lig4[11]:$lig4[12]") if ( $lig4[2] eq 'FE' or $lig4[2] eq 'S' );
   }
  }
 }
 
 my $fes2 = Chemistry::Mol->new(id => "FES1", name => "Fe-S");
 
 foreach my $wfes1 (@fes1)
 {
  my @fes3 = split(/\:/, $wfes1);
  
  $fes3[0] = 'Fe' if ( $fes3[0] eq 'FE' );
  
  my $fes4 = $fes2->new_atom(symbol => $fes3[0], coords => [$fes3[1],$fes3[2],$fes3[3]]);
 }
 
 find_bonds($fes2);
 
 my @fes5 = $fes2->separate();
 
 my $nfes5 = @fes5;
 
 
 # join ligands
 
 foreach my $wlig1 (@lig1)
 {
  my @lig3 = split(/\&/, $wlig1);
  
  open (PDB, ">$ttmpf02.pdb") || die "Cannot open $ttmpf02.pdb for writing\n";
  
  my $nat1 = 0;
  
  foreach my $wlig3 (@lig3)
  {
   my @lig4 = split(/\ /, $wlig3);
   my $nlig4 = @lig4;
   
   if ( $nlig4 > 5 )
   {
    if ( $lig4[2] ne 'H' )
    {
        if ( length($lig4[2]) == 1 )
     {
      $lig4[2] = ' '.$lig4[2].'  ';
     }
     elsif ( length($lig4[2]) == 2 )
     {
      $lig4[2] = ' '.$lig4[2].' ';
     }
     elsif ( length($lig4[2]) > 3 )
     {
      $lig4[2] = ' '.substr($lig4[2], 0, 3);
     }
     
     $lig4[2] = 'FE  ' if ( $lig4[2] eq ' FE ' );
     $lig4[2] = 'ZN  ' if ( $lig4[2] eq ' ZN ' );
     $lig4[2] = 'CA  ' if ( $lig4[2] eq ' CA ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' MG ' );
     $lig4[2] = 'CD  ' if ( $lig4[2] eq ' CD ' );
     $lig4[2] = 'CU  ' if ( $lig4[2] eq ' CU ' );
     $lig4[2] = 'MN  ' if ( $lig4[2] eq ' MN ' );
     $lig4[2] = 'NI  ' if ( $lig4[2] eq ' NI ' );
     $lig4[2] = 'BE  ' if ( $lig4[2] eq ' BE ' );
     $lig4[2] = 'AU  ' if ( $lig4[2] eq ' AU ' );
     $lig4[2] = 'CL  ' if ( $lig4[2] eq ' CL ' );
     $lig4[2] = 'BR  ' if ( $lig4[2] eq ' BR ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' RH ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' MO ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' RB ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' TB ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' TE ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' TL ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' EU ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' RE ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' GD ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' RU ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' TA ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' GA ' );
     $lig4[2] = 'MG  ' if ( $lig4[2] eq ' ZR ' );
     
     printf PDB ("HETATM%5d%5s%4s%6d%12.3f%8.3f%8.3f\n", ++$nat1, $lig4[2], 'LIG', 1, $lig4[10], $lig4[11], $lig4[12]);
    }
   }
  }
  
  close (PDB);
  
  if ( $nat1 >= $minat and $nat1 <= 200 and $nat1 <= $maxat )
  {
   my %opt = ();
   
   $opt{'sort'} = 1;
   $opt{'invariants'} = 0;
 
   my $tmol1 = Chemistry::Mol->read("$ttmpf02.pdb");
   
   canonicalize($tmol1, %opt);
   
   $tmol1->write("$ttmpf02.pdb");
  }
  
  elsif ( $nat1 >= $minat and $nat1 <= $maxat )
  {
   my @tlig1 = read_file("$ttmpf02.pdb");
   my @tlig2 = grep(/HETATM/, @tlig1);
   
   my @tlig3 = ();
   
   foreach my $wtlig2 (@tlig2)
   {
    if ( length($wtlig2) > 53 )
    {
     if ( substr($wtlig2, 0, 6) eq 'ATOM  ' or substr($wtlig2, 0, 6) eq 'HETATM' )
     {
      my $pat2 = substr($wtlig2, 12, 4);
      
      substr($pat2, 0, 1) = 'X' if ( substr($pat2, 0, 1) eq ' ' );
      
      push(@tlig3, $pat2.':'.$wtlig2);
     }
    }
   }
   
   @tlig2 = sort { lc($a) cmp lc($b) } @tlig3;
   
   open (PDB, ">$ttmpf02.pdb") || die "Cannot open $ttmpf02.pdb for writing\n";
   
    my $nat2 = 1;
    
    foreach my $wtlig2 (@tlig2)
    {
     my @tlig4 = split(/\:/, $wtlig2);
     
     if ( length($tlig4[1]) > 53 )
     {
      printf PDB ("HETATM%5d%5s%4s%6d%s\n", $nat2++, substr($tlig4[1], 12, 4), 'LIG', 1, substr($tlig4[1], 26, 28)) if ( substr($tlig4[1], 0, 6) eq 'ATOM  ' or substr($tlig4[1], 0, 6) eq 'HETATM' );
     }
    }
    
   close (PDB);
  }
  
  my @xlig1 = read_file("$ttmpf02.pdb");
  my @xlig2 = grep(/HETATM/, @xlig1);
  
  unlink("$ttmpf02.pdb") if ( -e "$ttmpf02.pdb" );
  
  my @ttyp=();
  my @ltyp=();
  my $ntyp=0;
  
  my @xyzok=();
  
  open (PDB, ">$ttmpf02.pdb") || die "Cannot open $ttmpf02.pdb for writing\n";
  
  foreach my $wxlig2 (@xlig2)
  {
   if ( length($wxlig2) > 52 )
   {
    if ( substr($wxlig2, 0, 6) eq 'HETATM' or substr($wxlig2, 0, 6) eq 'ATOM  ' )
    {
     if ( substr($wxlig2, 13, 1) ne "H"  and 
          substr($wxlig2, 13, 1) ne " "  and 
          substr($wxlig2, 12, 1) ne "H"  and 
          substr($wxlig2, 12, 1) ne "N"  and 
          substr($wxlig2, 12, 1) ne "P"  and 
          substr($wxlig2, 12, 1) ne "R"  and 
          substr($wxlig2, 12, 1) ne "I"  and 
          substr($wxlig2, 12, 1) ne "U"  and 
          substr($wxlig2, 12, 2) ne " K" and 
          substr($wxlig2, 12, 2) ne "NA" and 
          substr($wxlig2, 12, 2) ne "LI" and 
          substr($wxlig2, 12, 2) ne " W" and 
          substr($wxlig2, 12, 2) ne " U" and 
          substr($wxlig2, 12, 2) ne "IR" and 
          substr($wxlig2, 12, 2) ne "XE" and 
          substr($wxlig2, 12, 2) ne "KR" and 
          substr($wxlig2, 12, 2) ne "AR" and 
          substr($wxlig2, 12, 2) ne "SR" and 
          substr($wxlig2, 12, 2) ne "MO" and 
          substr($wxlig2, 12, 2) ne "BA" and 
          substr($wxlig2, 12, 2) ne "CS" and 
          substr($wxlig2, 12, 2) ne "AL" and 
          substr($wxlig2, 12, 2) ne " L" and 
          substr($wxlig2, 12, 2) ne "UN" and 
          substr($wxlig2, 12, 2) ne "TE" and 
          substr($wxlig2, 12, 2) ne "SM" and 
          substr($wxlig2, 12, 2) ne "SB" and 
          substr($wxlig2, 12, 2) ne "YB" and 
          substr($wxlig2, 12, 2) ne "GD" and 
          substr($wxlig2, 12, 2) ne " A" and 
          substr($wxlig2, 12, 2) ne "SC" and 
          substr($wxlig2, 12, 2) ne "TL" and 
          substr($wxlig2, 12, 2) ne " Y" and 
          substr($wxlig2, 12, 2) ne "LU" and 
          substr($wxlig2, 12, 2) ne "CB" and 
          substr($wxlig2, 12, 2) ne "CG" and 
          substr($wxlig2, 12, 2) ne "CC" and 
          substr($wxlig2, 12, 2) ne "CO" and 
          substr($wxlig2, 12, 2) ne "CP" and 
          substr($wxlig2, 12, 2) ne "CE" and 
          substr($wxlig2, 12, 2) ne "MC" and 
          substr($wxlig2, 12, 2) ne "GO" and 
          substr($wxlig2, 12, 2) ne "OS" and 
          substr($wxlig2, 12, 2) ne "CN" and 
          substr($wxlig2, 12, 2) ne "SG" and 
          substr($wxlig2, 12, 2) ne "TA" and 
          substr($wxlig2, 12, 2) ne " X" and 
          substr($wxlig2, 12, 2) ne "CX" and 
          substr($wxlig2, 12, 2) ne "CZ" and 
          substr($wxlig2, 12, 2) ne " D" and 
          substr($wxlig2, 12, 2) ne "GC" and 
          substr($wxlig2, 12, 2) ne "GN" and 
          substr($wxlig2, 12, 2) ne "**" and 
          substr($wxlig2, 12, 2) ne "_U" and 
          substr($wxlig2, 13, 3) ne "***" and 
          substr($wxlig2, 12, 2) ne "CO" )
     {
      my $at2 = substr($wxlig2, 12, 2);
      my $at1 = 0;
      
      substr($wxlig2, 17, 3) = 'LIG';
      
         if ($at2 eq " C") {$at1 = 1;}
      elsif ($at2 eq "AC") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "BC") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "ZC") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "CK") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "SI") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "C1") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "C2") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "C3") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "C4") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "C5") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "C6") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "C7") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "C8") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "C9") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "CR") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "CM") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "CF") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "0C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "1C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "2C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "3C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "4C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "5C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "6C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "7C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "8C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "9C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq "'C") {$at1 = 1; $at2=" C"; substr($wxlig2, 12, 2)=" C";}
      elsif ($at2 eq " O") {$at1 = 1;}
      elsif ($at2 eq "AO") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "BO") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "ZO") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "O1") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "O2") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "O3") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "O4") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "O5") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "O6") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "O7") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "O8") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "O9") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "1O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "2O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "3O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "4O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "5O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "6O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "7O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "8O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "9O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "OX") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq "'O") {$at1 = 1; $at2=" O"; substr($wxlig2, 12, 2)=" O";}
      elsif ($at2 eq " N") {$at1 = 1;}
      elsif ($at2 eq "AN") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "BN") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "1N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "2N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "3N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "4N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "5N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "6N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "7N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "8N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "9N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq "'N") {$at1 = 1; $at2=" N"; substr($wxlig2, 12, 2)=" N";}
      elsif ($at2 eq " P") {$at1 = 1;}
      elsif ($at2 eq "AP") {$at1 = 1; $at2=" P"; substr($wxlig2, 12, 2)=" P";}
      elsif ($at2 eq "OP") {$at1 = 1; $at2=" P"; substr($wxlig2, 12, 2)=" P";}
      elsif ($at2 eq "'P") {$at1 = 1; $at2=" P"; substr($wxlig2, 12, 2)=" P";}
      elsif ($at2 eq " S") {$at1 = 1;}
      elsif ($at2 eq "SE") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "S1") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "S2") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "S3") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "S4") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "S5") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "S6") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "S7") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "S8") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "S9") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "AS") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq "'9") {$at1 = 1; $at2=" S"; substr($wxlig2, 12, 2)=" S";}
      elsif ($at2 eq " F") {$at1 = 1;}
      elsif ($at2 eq "F2") {$at1 = 1; $at2=" F"; substr($wxlig2, 12, 2)=" F";}
      elsif ($at2 eq "F3") {$at1 = 1; $at2=" F"; substr($wxlig2, 12, 2)=" F";}
      elsif ($at2 eq "'F") {$at1 = 1; $at2=" F"; substr($wxlig2, 12, 2)=" F";}
      elsif ($at2 eq " I") {$at1 = 1;}
      elsif ($at2 eq "CL") {$at1 = 1;}
      elsif ($at2 eq " c") {$at1 = 1; $at2="CL"; substr($wxlig2, 12, 2)="CL";}
      elsif ($at2 eq "BR") {$at1 = 1;}
      elsif ($at2 eq " b") {$at1 = 1; $at2="BR"; substr($wxlig2, 12, 2)="BR";}
      elsif ($at2 eq "FE") {$at1 = 0;}
      elsif ($at2 eq "ZN") {$at1 = 0;}
      elsif ($at2 eq "CA") {$at1 = 0;}
      elsif ($at2 eq "MG") {$at1 = 0;}
      elsif ($at2 eq "CD") {$at1 = 0;}
      elsif ($at2 eq "CU") {$at1 = 0;}
      elsif ($at2 eq "MN") {$at1 = 0;}
      elsif ($at2 eq "NI") {$at1 = 0;}
      elsif ($at2 eq "BE") {$at1 = 0;}
      elsif ($at2 eq "AU") {$at1 = 0;}
      elsif ($at2 eq " V") {$at1 = 0;}
      elsif ($at2 eq " B") {$at1 = 0;}
      elsif ($at2 eq "RH") {$at1 = 0;}
      elsif ($at2 eq "MO") {$at1 = 0;}
      elsif ($at2 eq "RB") {$at1 = 0;}
      elsif ($at2 eq "TB") {$at1 = 0;}
      
      else
      {
       die "In $fpdb1 unknown atom: -$at2-\n$wxlig2\n\n";
      }

      my $nw = 1;
      
      my $xa = 0;
      
      for ( $xa = 0; $xa < $ntyp; $xa++ )
      {
       if ( $ttyp[$xa] eq $at2 )
       {
        $ltyp[$xa]++;
        $nw = 0;
        last;
       }
      }

      if ( $nw )
      {
       push(@ttyp, $at2);
       $ltyp[$ntyp] = 1;

       if ( $at1 )
       {
        substr($wxlig2, 12, 4) = $at2.$ltyp[$ntyp].' ' if ( $ltyp[$ntyp] <  10 );
        substr($wxlig2, 12, 4) = $at2.$ltyp[$ntyp]     if ( $ltyp[$ntyp] >= 10 );
       }
       else
       {
        substr($wxlig2, 12, 4) = $at2.'  ';
       }
       
       $ntyp++;
      }
      
      else
      {
       if ( $at1 )
       {
        $ltyp[$xa] -= 99 if ( $ltyp[$xa] >= 100 );
        
        substr($wxlig2, 12, 4) = $at2.$ltyp[$xa].' ' if ( $ltyp[$xa] <  10 );
        substr($wxlig2, 12, 4) = $at2.$ltyp[$xa]     if ( $ltyp[$xa] >= 10 );
       }
       else
       {
        substr($wxlig2, 12, 4) = $at2.'  ';
       }
      }

      my $w66 = 1;

      my $nx = substr($wxlig2, 30, 8) * 1.0;
      my $ny = substr($wxlig2, 38, 8) * 1.0;
      my $nz = substr($wxlig2, 46, 8) * 1.0;
      
      foreach my $wxyz ( @xyzok )
      {
       my $ox = substr($wxyz, 30, 8) * 1.0;
       my $oy = substr($wxyz, 38, 8) * 1.0;
       my $oz = substr($wxyz, 46, 8) * 1.0;
      
       my $rr = sqrt((($nx-$ox)**2.0)+(($ny-$oy)**2.0)+(($nz-$oz)**2.0));
      
       $w66 = 0 if ( $rr < 1.0 );
      }
      
      if ( $w66 )
      {
       push(@xyzok, $wxlig2);
       
       substr($wxlig2, 0, 6) = "HETATM" if ( substr($wxlig2, 0, 6) eq "ATOM  " );
       
       my $desc = substr($wxlig2, 0, 30);
       
       my $kx = substr($wxlig2, 30, 8) * 1.0;
       my $ky = substr($wxlig2, 38, 8) * 1.0;
       my $kz = substr($wxlig2, 46, 8) * 1.0;
       
       printf PDB ("%s%8.3f%8.3f%8.3f\n", $desc, $kx, $ky, $kz);
      }
     }
    }
   }
  }
  
  close (PDB);
  
  my @xlig3 = read_file("$ttmpf02.pdb");
  my @xlig4 = grep(/HETATM/, @xlig3);
  my $nxlig4 = @xlig4;
  
  unlink("$ttmpf02.pdb") if ( -e "$ttmpf02.pdb" );
  
  if ( $nxlig4 <= $maxat and $nxlig4 >= $minat )
  {
   # LPC
   
   foreach my $wcha1 (@cha1)
   {
    my @xpro1 = read_file("$fpdb1$wcha1.pdb");
    my @xpro2 = grep(/ATOM/, @xpro1);
    
    write_file("$ttmpf03.pdb", @xpro2);
    append_file("$ttmpf03.pdb", "TER\n");
    append_file("$ttmpf03.pdb", @xlig4);
    append_file("$ttmpf03.pdb", "END\n");
    
    unlink("RES1") if ( -e "RES1" );
    
    `$lpc 1 $ttmpf03.pdb`;
    
    open (LPC, "RES1") || die "Cannot open RES1 for reading.\n";
     my @lpc1=<LPC>;
     chomp(@lpc1);
    close (LPC);
    
    unlink "RES1" if ( -e "RES1");
    unlink "$ttmpf03.pdb" if ( -e "$ttmpf03.pdb");
    
    my $ww1 = 0;
    my $nn1 = 0;
    
    foreach my $wlpc1 (@lpc1)
    {
     $ww1++ if ( $wlpc1 eq '----------------------------------------------------------' or 
                 $wlpc1 eq '                                  Specific contacts' or 
                 $wlpc1 eq '     Residue      Dist    Surf    HB    Arom    Phob    DC' );
     
     $ww1 = 0 if ( !length($wlpc1) );
     
     $nn1++ if ( $ww1 == 4 and $wlpc1 ne '----------------------------------------------------------' );
    }
    
    if ( $nn1 >= $minbr )
    {
     # print ligand
     
     my $pat3 = $fpdb1.$wcha1;
     
       if ( $cha2{$wcha1} < 10 )
     {
      $pat3 .= '.LG'.$cha2{$wcha1};
     }
     elsif ( $cha2{$wcha1} < 100 )
     {
      $pat3 .= '.L'.$cha2{$wcha1};
     }
     else
     {
      $pat3 .= '.'.$cha2{$wcha1};
     }
     
     # check if ligand overlaps with previously found
     
     my $w4 = 1;
     
     my @olig1 = split(/\&/, $cha3{$wcha1});
     
     foreach my $wolig1 (@olig1)
     {
      my @olig2 = split(/\:/, $wolig1);
      
      foreach my $wolig2 (@olig2)
      {
       if ( length($wolig2) > 53 )
       {
        if ( substr($wolig2, 0, 6) eq 'HETATM' )
        {
         my $x1 = substr($wolig2, 30, 8) * 1.0;
         my $y1 = substr($wolig2, 38, 8) * 1.0;
         my $z1 = substr($wolig2, 46, 8) * 1.0;
         
         foreach my $wxlig4 (@xlig4)
         {
          if ( length($wxlig4) > 53 )
          {
           if ( substr($wxlig4, 0, 6) eq 'HETATM' )
           {
            my $x2 = substr($wxlig4, 30, 8) * 1.0;
            my $y2 = substr($wxlig4, 38, 8) * 1.0;
            my $z2 = substr($wxlig4, 46, 8) * 1.0;
            
            $w4 = 0 if ( sqrt(($x1-$x2)**2 + ($y1-$y2)**2 + ($z1-$z2)**2) < 1.0 );
           }
          }
          
          last if ( $w4 );
         }
        }
       }
       
       last if ( $w4 );
      }
      
      last if ( $w4 );
     }
     
     if ( $w4 )
     {
      write_file("$pat3.pdb", @xlig4);
      
      push(@files1, "$pat3.pdb");
      
      $cha2{$wcha1}++;
      
      my $pat4 = '';
      
      foreach my $wxlig4 (@xlig4)
      {
       $pat4 .= $wxlig4.':';
      }
      
      substr($pat4, -1, 1) = '' if ( substr($pat4, -1, 1) eq ':' );
      
      if ( length($cha3{$wcha1}) )
      {
       $cha3{$wcha1} .= '&'.$pat4;
      }
      else
      {
       $cha3{$wcha1} = $pat4;
      }
     }
    }
   }
  }
 }
 
 
 # metals
 
 foreach my $wcha1 (@cha1)
 {
  my @mpro1 = read_file("$fpdb1$wcha1.pdb");
  my @mpro2 = grep(/ATOM/, @mpro1);
  
  my $met3 = "$fpdb1$wcha1.metals";
  
  my @met4 = ();
  
  foreach my $wmet1 (@met1)
  {
   # LPC
   
   my @met2 = split(/\:/, $wmet1);
   
   my @mlig1 = sprintf("HETATM%5d%3s%6s%6d%12.3f%8.3f%8.3f\n", 1, $met2[0], $met2[0], 1, $met2[1], $met2[2], $met2[3]);
   
   write_file("$ttmpf04.pdb", @mpro2);
   append_file("$ttmpf04.pdb", "TER\n");
   append_file("$ttmpf04.pdb", @mlig1);
   append_file("$ttmpf04.pdb", "END\n");
   
   unlink("RES1") if ( -e "RES1" );
   
   `$lpc 1 $ttmpf04.pdb`;
   
   open (LPC, "RES1") || die "Cannot open RES1 for reading.\n";
    my @mlpc1=<LPC>;
    chomp(@mlpc1);
   close (LPC);
   
   unlink "RES1" if ( -e "RES1");
   unlink "$ttmpf04.pdb" if ( -e "$ttmpf04.pdb");
   
   my $mww1 = 0;
   
   my @met5 = ();
   
   foreach my $wmlpc1 (@mlpc1)
   {
    $mww1++ if ( $wmlpc1 eq '----------------------------------------------------------' or 
                 $wmlpc1 eq '                                  Specific contacts' or 
                 $wmlpc1 eq '     Residue      Dist    Surf    HB    Arom    Phob    DC' );
    
    $mww1 = 0 if ( !length($wmlpc1) );
    
    push(@met5, (substr($wmlpc1, 0, 7) * 1).':'.three2one(substr($wmlpc1, 10, 3))) if ( $mww1 == 4 and $wmlpc1 ne '----------------------------------------------------------' );
   }
   
   my $met6 = @met5;
   
   
   if ( $met6 >= $minmt )
   {
    my $met7 = sprintf("%2s%9.3f%9.3f%9.3f%3d", $met2[0], $met2[1], $met2[2], $met2[3], $met6);
    
    foreach my $wmet5 (@met5)
    {
     $met7 .= " $wmet5";
    }
    
    push(@met4, "$met7\n");
   }
  }
  
  if ( @met4 )
  {
   write_file($met3, @met4);
   
   push(@files1, $met3);
  }
 }
 
 
 # Fe-S
 
 foreach my $wcha1 (@cha1)
 {
  my @spro1 = read_file("$fpdb1$wcha1.pdb");
  my @spro2 = grep(/ATOM/, @spro1);
  
  my $nspro2 = @spro2;
  
  my $spro3 = substr($spro2[$nspro2-1],  6, 5) * 1;
  my $spro4 = substr($spro2[$nspro2-1], 22, 4) * 1;
  
  my $fes6 = "$fpdb1$wcha1.fes";
  
  my @fes7 = ();
  
  foreach my $wfes5 (@fes5)
  {
   # LPC
   
   my @fes8 = $wfes5->atoms();
   my @fes15 = $wfes5->bonds();
   
   my @fes9 = ();
   
   my $fes10 = $spro3;
   
   foreach my $wfes8 (@fes8)
   {
    my $fes11 = ' FE ';
    
    $fes11 = ' S  ' if ( $wfes8->symbol() eq 'S' );
    
    push(@fes9, sprintf("HETATM%5d%5s%4s%6d%12.3f%8.3f%8.3f\n", ++$fes10, $fes11, 'FES', $spro4+1, $wfes8->x3(), $wfes8->y3(), $wfes8->z3()));
   }
   
   write_file("$ttmpf05.pdb", @spro2);
   append_file("$ttmpf05.pdb", "TER\n");
   append_file("$ttmpf05.pdb", @fes9);
   append_file("$ttmpf05.pdb", "END\n");
   
   unlink("RES1") if ( -e "RES1" );
   
   `$lpc 1 $ttmpf05.pdb`;
   
   open (LPC, "RES1") || die "Cannot open RES1 for reading.\n";
    my @mlpc1=<LPC>;
    chomp(@mlpc1);
   close (LPC);
   
   unlink "RES1" if ( -e "RES1");
   unlink "$ttmpf05.pdb" if ( -e "$ttmpf05.pdb");
   
   my $mww1 = 0;
   
   my @fes12 = ();
   
   foreach my $wmlpc1 (@mlpc1)
   {
    $mww1++ if ( $wmlpc1 eq '----------------------------------------------------------' or 
                 $wmlpc1 eq '                                  Specific contacts' or 
                 $wmlpc1 eq '     Residue      Dist    Surf    HB    Arom    Phob    DC' );
    
    $mww1 = 0 if ( !length($wmlpc1) );
    
    push(@fes12, (substr($wmlpc1, 0, 7) * 1).':'.three2one(substr($wmlpc1, 10, 3)).':'.(substr($wmlpc1, 14, 8) * 1.0).':'.(substr($wmlpc1, 22, 8) * 1.0)) if ( $mww1 == 4 and $wmlpc1 ne '----------------------------------------------------------' );
   }
   
   my $fes13 = @fes12;
   
   if ( $fes13 >= $minfs )
   {
    push(@fes7, "Fs-S start\n");
    
    foreach my $wfes8 (@fes8)
    {
     my $fes14 = $wfes8->id();
     
     $fes14 =~ s/a//g;
     
     push(@fes7, sprintf("ATM %3d %5s %8.3f %8.3f %8.3f\n", $fes14, $wfes8->symbol(), $wfes8->x3(), $wfes8->y3(), $wfes8->z3()));
    }
    
    my $fes16 = 0;
    
    foreach my $wfes15 (@fes15)
    {
     my @fes14 = sort { $a <=> $b } $wfes15->atoms();
     
     my $fes15 = '';
     
     foreach my $wfes14 (@fes14)
     {
      $wfes14 =~ s/a//g;
      
      $fes15 .= sprintf("%3d", $wfes14);
     }
     
     push(@fes7, sprintf("BND %3d %s\n", ++$fes16, $fes15));
    }
    
    my $fes18 = 0;
    
    foreach my $wfes12 (@fes12)
    {
     my @fes17 = split(/\:/, $wfes12);
     
     push(@fes7, sprintf("RES %3d %4d %s %6.1f %6.1f\n", ++$fes18, $fes17[0], $fes17[1], $fes17[2], $fes17[3]));
    }
    
    push(@fes7, "Fs-S end\n");
   }
  }
  
  if ( @fes7 )
  {
   write_file($fes6, @fes7);
   
   push(@files1, $fes6);
  }
 }
 
 
 # water
 
 foreach my $wcha1 (@cha1)
 {
  my @mpro1 = read_file("$fpdb1$wcha1.pdb");
  my @mpro2 = grep(/ATOM/, @mpro1);
  
  my $wat2 = "$fpdb1$wcha1.water";
  
  my @wat3 = ();
  
  foreach my $wwat1 (@wat1)
  {
   my @wat4 = split(/\&/, $wwat1);
   
   foreach my $wwat4 (@wat4)
   {
    my @wat5 = split(/\ /, $wwat4);
    my $nwat5 = @wat5;
    
    if ( $wat5[2] eq 'O' )
    {
     my $wat6 = 0;
     
     foreach my $wmpro2 (@mpro2)
     {
      my $mx1 = substr($wmpro2, 30, 8) * 1.0;
      my $my1 = substr($wmpro2, 38, 8) * 1.0;
      my $mz1 = substr($wmpro2, 46, 8) * 1.0;
      
      my $mr1 = sqrt(($wat5[10]-$mx1)**2 + ($wat5[11]-$my1)**2 + ($wat5[12]-$mz1)**2);
      
      if ( $mr1 <= 4.5 )
      {
       $wat6 = 1;
       
       last;
      }
     }
     
     push(@wat3, sprintf("WAT%9.3f%9.3f%9.3f\n", $wat5[10], $wat5[11], $wat5[12]));
    }
   }
  }
  
  if ( @wat3 )
  {
   write_file($wat2, @wat3);
   
   push(@files1, $wat2);
  }
 }
 
 foreach my $wfiles1 (@files1)
 {
  move($wfiles1, "$cwd01/$wfiles1");
 }
 
 chdir($cwd01) or die "Cannot chdir to $cwd01 $!";
 
 rmtree($scratch);
 
 exit(0);
 
 
