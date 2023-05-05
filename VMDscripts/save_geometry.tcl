# Execution
# vmd mol_0_traj.xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_geometry.tcl


# load functions
source /home/ovillegas/.gitproyects/PhotonasticMat/VMDscripts/FunctionsPhotoN.tcl

# load the labels and arrays
distances_STAMP
dihedrals_STAMP

save_GEOMETRY_STAMP

puts "Finish!"
exit
