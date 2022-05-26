mol load gro Microfibril-15X20.gro
set sel [atomselect top "same residue as (within 100 of resname C0810)"]
$sel writegro Microfibril.gro
exit
