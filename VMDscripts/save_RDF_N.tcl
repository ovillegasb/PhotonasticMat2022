

set current_directory [pwd]
regexp {/(cis|trans)/(\d+)} $current_directory match type number

# RDF analysis
# all atoms
# only for GRO files
set RDFs [measure gofr [atomselect 0 "resid 1 to $number"] [atomselect 0 "not resid 1 to $number"] delta 0.1 rmax 20. usepbc True first 1 last -1 step 1]
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