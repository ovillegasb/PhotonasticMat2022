
# RUN:
# vmd mol_traj.xyz -dispdev text -e path/save_polymerInfo.tcl

label add Dihedrals 0/11 0/12 0/13 0/14
set angles [label graph Dihedrals 0]
set N [llength $angles]

set ofile [open "dihedrals.dat" w]

for {set i 0} {$i < $N} {incr i} {
    # angle i
    set phi [lindex $angles $i]
    puts $ofile "$i    $phi"
}

puts "Finish!"
close $ofile
exit
