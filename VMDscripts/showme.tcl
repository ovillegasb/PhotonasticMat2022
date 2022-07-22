
# Show BOX
pbc box -color black

# Style
mol modstyle 0 0 VDW 1.000000 12.000000
mol modmaterial 0 0 AOEdgy

# Background
color Display Background white

# Animation
animate speed 0.78
animate goto 0

# Color Labels
color Labels Bonds black
color Labels Angles orange
color Labels Dihedrals red3
