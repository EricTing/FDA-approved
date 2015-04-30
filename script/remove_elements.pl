use strict;
use Chemistry::Mol;
use Chemistry::Atom;
use Chemistry::File::SDF;

my @to_fix = ('DB08928'
              ,'DB01482'
              ,'DB01294'
              ,'DB06715'
              ,'DB00534'
              ,'DB01377'
              ,'DB00971'
              ,'DB00995'
              ,'DB00325'
              ,'DB08990'
              ,'DB05245'
              ,'DB00364'
              ,'DB00638'
              ,'DB05630'
              ,'DB01375'
              ,'DB00200'
              ,'DB00188'
              ,'DB00526'
              ,'DB01169'
              ,'DB00483'
              ,'DB00115'
              ,'DB00516'
              ,'DB00515'
              ,'DB08913'
              ,'DB06723'
              ,'DB05389'
              ,'DB06402'
              ,'DB00958');

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

foreach (@to_fix) {
    my $ifn = "/work/jaydy/working/obgen_todo/" . $_ . ".sdf";
    my $ofn = "/work/jaydy/working/obgen_fix/" . $_ . ".sdf";

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
}
