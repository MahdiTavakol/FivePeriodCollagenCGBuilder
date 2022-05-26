mol load pdb 0-3hr2-aligned.pdb
set numRes 1054
set File [open "0-3hr2-AA.csv" w]
for {set i 0} {$i < $numRes} {incr i} {
	set j [expr $i+1]
	set atoms [atomselect top "resid $j"]
	set x [lindex [measure center $atoms] 0]
	set y [lindex [measure center $atoms] 1]
	set z [lindex [measure center $atoms] 2]
	puts $File "$j,$x,$y,$z"
}
close $File
exit
