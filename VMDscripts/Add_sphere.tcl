
# draw material Transparent
# draw color gray
# draw sphere {0.0 0.0 0} radius 10.0 resolution 30
# 
# draw color black
# draw text {10.0 10.0 0.0} " 1 nm"

# proc vmd_draw_sphere {mol start end} {
#     # an sphere is made using box center.
#     # set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
#     graphics $mol cylinder $start $middle radius 0.15
#     graphics $mol cone $middle $end radius 0.25
# }

# get the number of frames in the movie
###set num [molinfo top get numframes]
###puts $num

# loop through the frames
# for {set i 0} {$i < $num} {incr i} {
#     # go to the given frame
#     #set filename Frame[format "%04d" $i].gif
#     animate goto $i
#     #set sel [atomselect top "segid P2"]
#     # set com [measure center $sel]
#     # $sel moveby [vecscale -1.0 $com]
#     # set I [draw principalaxes $sel]
#     # set A [orient $sel [lindex $I 0] {0 1 0}]
#     # $sel move $A 
#     # set I [draw principalaxes $sel]
#     # render snapshot $filename 
#     draw material Transparent
#     draw color gray
#     draw sphere {0.0 0.0 0} radius 10.0 resolution 30
# }

graphics top materials on
graphics top material Transparent
graphics top color silver
graphics top sphere {0.0 0.0 0} radius 5.0 resolution 30