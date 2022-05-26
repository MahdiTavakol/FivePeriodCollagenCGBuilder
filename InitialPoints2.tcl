mol load pdb 0-3hr2-aligned.pdb
set numAtom 3134
set File [open "0-3hr2-AA.csv" w]
for {set i 0} {$i < $numAtom} {incr i} {
	set j [expr $i+1]
	set atoms [atomselect top "serial $j"]
	set x [lindex [measure center $atoms] 0]
	set y [lindex [measure center $atoms] 1]
	set z [lindex [measure center $atoms] 2]
	puts $File "$j,$x,$y,$z"
}
close $File
exit
