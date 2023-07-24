
# RUN:
# vmd mol_traj.xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_dih_data.tcl

proc save_RDF {} {
    # RDF analysis
    # all atoms
    # only for GRO files
    set RDFs [measure gofr [atomselect 0 "resid 1"] [atomselect 0 "not resid 1"] delta 0.1 rmax 20. usepbc True first 1 last -1 step 1]
    set r [lindex $RDFs 0]
    set g_r [lindex $RDFs 1]
    set nbins [llength $r]
    set ofile [open "rdf_cm_pc_env.csv" w]
    puts $ofile "r,g_r"
    for {set i 0} {$i < $nbins} {incr i} {
        puts $ofile "[lindex $r $i],[lindex $g_r $i]"
    }
    close $ofile
}

# call function
save_RDF

puts "Finish!"
exit