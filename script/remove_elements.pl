use strict;
use Chemistry::Mol;
use Chemistry::Atom;
use Chemistry::File::SDF;


my @forbidden = ('B'
                 ,'Ga'
                 ,'Pt'
                 ,'Na'
                 ,'K'
                 ,'Mg'
                 ,'Si'
                 ,'Ra'
                 ,'Sb'
                 ,'Li'
                 ,'Se'
                 ,'Bi'
                 ,'Zn'
                 ,'Co'
                 ,'Ag'
                 ,'Ca'
                 ,'Al'
                 ,'As'
                 ,'Gd'
                 ,'Au'
                 ,'Fe'
                 ,'Hg');

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];

my $drug = Chemistry::Mol->read($ifn);
my @all_atoms = $drug->atoms;

foreach my $atom (@all_atoms) {
    my $symbol = $atom->symbol;
    if ($symbol ~~ @forbidden)
    {
        printf ("delete %s in %s\n", $symbol, $ifn);
        $drug->delete_atom($atom);
    }
    $drug->write($ofn);
}
