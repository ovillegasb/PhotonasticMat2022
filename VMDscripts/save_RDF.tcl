
# RUN:
# vmd *.gro -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_RDF.tcl

# load functions
source ~/.gitproyects/PhotonasticMat/VMDscripts/FunctionsPhotoN.tcl


# call function
save_RDF

puts "Finish!"
exit
