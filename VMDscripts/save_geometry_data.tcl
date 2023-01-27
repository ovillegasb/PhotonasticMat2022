
# RUN:
# vmd mol_traj.xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_dih_data.tcl

# GEOMETRY STAMP
# r_cc
label add Bonds 0/11 0/14
# r_ar_1
label add Bonds 0/0 0/24
# r_ar_2
label add Bonds 0/9 0/15
# d_cnnc
label add Dihedrals 0/11 0/12 0/13 0/14
# d_ccnn_1,2,3,4
label add Dihedrals 0/0 0/11 0/12 0/13
label add Dihedrals 0/9 0/11 0/12 0/13
label add Dihedrals 0/12 0/13 0/14 0/24
label add Dihedrals 0/12 0/13 0/14 0/15


# GROMACS label add Dihedrals 0/4 0/10 0/11 0/12

set r_cc [label graph Bonds 0]
set r_ar_1 [label graph Bonds 1]
set r_ar_2 [label graph Bonds 2]

set d_cnnc [label graph Dihedrals 0]
set d_ccnn_1 [label graph Dihedrals 1]
set d_ccnn_2 [label graph Dihedrals 2]
set d_ccnn_3 [label graph Dihedrals 3]
set d_ccnn_4 [label graph Dihedrals 4]

set Nframes [molinfo 0 get numframes]
set ofile [open "Geometry_mol.csv" w]
puts $ofile "frame,r_cc,r_ar_1,r_ar_2,d_cnnc,d_ccnn_1,d_ccnn_2,d_ccnn_3,d_ccnn_4"
for {set i 0} {$i < $Nframes} {incr i} {
    puts $ofile "$i,[lindex $r_cc $i],[lindex $r_ar_1 $i],[lindex $r_ar_2 $i],[lindex $d_cnnc $i],[lindex $d_ccnn_1 $i],[lindex $d_ccnn_2 $i],[lindex $d_ccnn_3 $i],[lindex $d_ccnn_4 $i]"
}
 
puts "Finish!"
close $ofile
exit
