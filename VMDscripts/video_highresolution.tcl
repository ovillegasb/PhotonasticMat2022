

# script example:
# /usr/local/lib/vmd/plugins/noarch/tcl/vmdmovie1.9/vmdmovie.tcl

# Size of windows resolution
# example: 1680 1013
display get size

# Command ffmpeg
# vcode libx265 personal, libx264 for shared
# 
# "ffmpeg -an -i filename.%05d.ppm -vcodec libx264 -r 24 filename.mp4"