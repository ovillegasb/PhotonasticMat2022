

# Style
mol modstyle 0 0 VDW 0.800000 50.000000
mol modmaterial 0 0 AOEdgy
color Name C gray
color Name N blue2

# Background
color Display Background white

# Display options
display depthcue on
display projection Orthographic
display shadows on
display ambientocclusion on

pbc box off

# Meterials
#material add copy AOEdgy
#material rename Material23 matOV

# render TachyonInternal test3.tga display %s
# In command render
# "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 1600 1200 -o %s.tga
# or in terminal tcl
# ender TachyonInternal test3.tga "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 1600 1200 -o %s.tga


# Others
# Select solvent
# mol modselect 1 0 resname SOL
# mol modmaterial 1 0 Material23
# material add Material23 copy AOEdgy
# material change diffuse Material23 0.300000
# material change specular Material23 0.500000
# material change shininess Material23 0.000000
# material change opacity Material23 0.050000
# material change outline Material23 0.000000
# material change outlinewidth Material23 0.000000