
# RUN:
#  vmd xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/countTypes.tcl


# COLOR CODE:
#------------
# Donneur                 Kr
# AccepteurLibre          Xe
# Graine                  Br
# AncienneGraine          F
# AncienAlibre            Ca
# AncienA                 Fe
# AncienD                 B


# Number of frames
set nframes [molinfo top get numframes]
# puts "N frames $nframes"


# Current frame number
set frame [molinfo top get frame]
# puts "current frames $frame"


# Graine : name Br
# for show only
# mol modselect 0 0 name Br
set grain [atomselect top "name Br" frame $frame]
puts "N Graine [$grain num]"


# Donneur : name Kr
# for show only
# mol modselect 0 0 name Kr
set donor [atomselect 0 "name Kr" frame $frame]
puts "N Donneur [$donor num]"

# AccepteurLibre : name Xe
# for show only
# mol modselect 0 0 name Xe
set acceptorL [atomselect top "name Xe" frame $frame]
puts "N AccepteurLibre [$acceptorL num]"


# AncienneGraine : name F
# for show only
# mol modselect 0 0 name F
# set ancgrain [atomselect top "name F"]
# puts "N AncienneGraine [$ancgrain num]"

# exit