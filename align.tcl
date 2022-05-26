lappend auto_path ~/Downloads/la1.0
lappend auto_path ~/Downloads/orient
mol load pdb 3hr2.pdb
  package require Orient
namespace import Orient::orient

set sel [atomselect top all]
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
set I [draw principalaxes $sel]
set sel [atomselect top all]
$sel writepdb  0-3hr2-aligned.pdb
exit
