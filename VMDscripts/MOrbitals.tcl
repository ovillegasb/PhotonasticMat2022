
package require pbctools
pbc box -off
color Display Background white


# From cube loaded
mol modstyle 0 0 VDW 0.200000 32.000000
# new representation
mol color Name
mol representation DynamicBonds 1.600000 0.100000 12.000000
mol selection all
mol material AOEdgy
mol addrep 0

# Colors
color Name C gray
color Name N blue2

# Display options
display depthcue on
display projection Orthographic
display shadows on
display ambientocclusion on

# Isosurface
mol color ColorID 0
mol representation Isosurface -0.020000 0 0 0 1 1
mol selection all
mol material Transparent
mol addrep 0

mol color ColorID 1
mol representation Isosurface 0.020000 0 0 0 1 1
mol selection all
mol material Transparent
mol addrep 0

puts "Its ok!"


# mol new {/home/ovillegas/.bettyboop/.jeanzay/AZOB_TD/isomer-cis/48.cube} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 }
# # mol modstyle 0 0 CPK 1.000000 0.300000 12.000000 12.000000
# mol modstyle 0 0 DynamicBonds 1.600000 0.100000 12.000000
# 
# mol modmaterial 0 0 AOEdgy
# 
# mol modstyle 0 0 DynamicBonds 1.600000 0.100000 12.000000
# mol color Name
# mol representation DynamicBonds 1.600000 0.100000 12.000000
# mol selection all
# mol material AOEdgy
# mol addrep 0
# mol modstyle 1 0 VDW 0.200000 12.000000