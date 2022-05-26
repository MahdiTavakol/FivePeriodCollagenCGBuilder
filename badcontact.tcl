mol load gro 4-MT-Mineralized.gro
set sel [atomselect top all]
$sel set beta 0
set bad [atomselect top {resid > 155 and ( within 12.46 of resid < 156)}]
$bad set beta 1
set out [atomselect top {beta 0}]
$out writepdb 4-MT-Mineralized-12.46.pdb
mol delete top
mol load pdb 4-MT-Mineralized-12.46.pdb
set sel [atomselect top all]
$sel writegro 4-MT-Mineralized-12.46.gro
mol delete top
exit
