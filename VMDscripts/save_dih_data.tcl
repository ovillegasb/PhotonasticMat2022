
# RUN:
# vmd mol_traj.xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_dih_data.tcl

label add Dihedrals 0/4 0/10 0/11 0/12
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
