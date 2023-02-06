
# RUN:
# vmd *.xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_dih_data.tcl

# PBC in xyz files
pbc set {65.002 60. 60.} -all
pbc box -center origin

# RDF analysis
# all atoms
set RDFs [measure gofr [atomselect 0 "index 0 to 25"] [atomselect 0 "not index 0 to 25"] delta 0.1 rmax 20. usepbc True first 0 last -1 step 1]
set r [lindex $RDFs 0]
set g_r [lindex $RDFs 1]
set nbins [llength $r]

set ofile [open "rdf_all_at_pc_env.csv" w]
puts $ofile "r,g_r"
for {set i 0} {$i < $nbins} {incr i} {
    puts $ofile "[lindex $r $i],[lindex $g_r $i]"
}
close $ofile


puts "Finish!"
exit
