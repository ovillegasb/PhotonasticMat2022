
# RUN:
# vmd mol_traj.xyz -dispdev text -e save_geometry_N_pc.tcl


# RUN:
# vmd mol_traj.xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_dih_data.tcl

# GEOMETRY GROMACS
# r_c--c
label add Bonds 0/4 0/12
# r_ar_1
label add Bonds 0/3 0/17
# r_ar_2
label add Bonds 0/5 0/13
# r_co_1
label add Bonds 0/1 0/24
# r_co_2
label add Bonds 0/15 0/22
# r_oh_1
label add Bonds 0/5 0/6
# r_oh_2
label add Bonds 0/22 0/23
# r_nn
label add Bonds 0/10 0/12
# r_cn_1
label add Bonds 0/4 0/10
# r_cn_2
label add Bonds 0/12 0/11
# r_ar_3
label add Bonds 0/1 0/15

# a_cnn_1,2
label add Angles 0/4 0/10 0/11
label add Angles 0/13 0/12 0/11
# a_coh_1,2
label add Angles 0/1 0/24 0/25
label add Angles 0/15 0/22 0/23

# d_cnnc
label add Dihedrals 0/4 0/10 0/11 0/12
# d_ccnn_1,2,3,4
# 0.0
label add Dihedrals 0/3 0/4 0/10 0/11
label add Dihedrals 0/13 0/12 0/11 0/10
# 180.
label add Dihedrals 0/5 0/4 0/10 0/11
label add Dihedrals 0/17 0/12 0/11 0/10


set r_cc [label graph Bonds 0]
set r_ar_1 [label graph Bonds 1]
set r_ar_2 [label graph Bonds 2]
set r_co_1 [label graph Bonds 3]
set r_co_2 [label graph Bonds 4]
set r_oh_1 [label graph Bonds 5]
set r_oh_2 [label graph Bonds 6]
set r_nn [label graph Bonds 7]
set r_cn_1 [label graph Bonds 8]
set r_cn_2 [label graph Bonds 9]
set r_ar_3 [label graph Bonds 10]

set a_cnn_1 [label graph Angles 0]
set a_cnn_2 [label graph Angles 1]
set a_coh_1 [label graph Angles 2]
set a_coh_2 [label graph Angles 3]

set d_cnnc [label graph Dihedrals 0]
set d_ccnn_1 [label graph Dihedrals 1]
set d_ccnn_2 [label graph Dihedrals 2]
set d_ccnn_3 [label graph Dihedrals 3]
set d_ccnn_4 [label graph Dihedrals 4]

set input [molinfo top get filename]
regexp {mol_traj_(\d+)} $input match resid
set output "Geometry_mol_$resid.csv"

set Nframes [molinfo 0 get numframes]
set ofile [open $output w]
puts $ofile "frame,r_cc,r_ar_1,r_ar_2,r_co_1,r_co_2,r_oh_1,r_oh_2,r_nn,r_cn_1,r_cn_2,r_ar_3,a_cnn_1,a_cnn_2,a_coh_1,a_coh_2,d_cnnc,d_ccnn_1,d_ccnn_2,d_ccnn_3,d_ccnn_4"
for {set i 0} {$i < $Nframes} {incr i} {
    puts $ofile "$i,[lindex $r_cc $i],[lindex $r_ar_1 $i],[lindex $r_ar_2 $i],[lindex $r_co_1 $i],[lindex $r_co_2 $i],[lindex $r_oh_1 $i],[lindex $r_oh_2 $i],[lindex $r_nn $i],[lindex $r_cn_1 $i],[lindex $r_cn_2 $i],[lindex $r_ar_3 $i],[lindex $a_cnn_1 $i],[lindex $a_cnn_2 $i],[lindex $a_coh_1 $i],[lindex $a_coh_2 $i],[lindex $d_cnnc $i],[lindex $d_ccnn_1 $i],[lindex $d_ccnn_2 $i],[lindex $d_ccnn_3 $i],[lindex $d_ccnn_4 $i]"
}
close $ofile

puts "Finish!"
exit
