
# Start loading all xyz files

# Transparent material
mol modmaterial 0 0 Transparent
mol modstyle 0 0 DynamicBonds 1.600000 0.300000 12.000000

# Adding new representation
mol color Name
mol representation DynamicBonds 1.600000 0.300000 12.000000
mol selection index 0 to 9 90 to 99 1130 to 1139 1940 to 1949 2310 to 2319 2720 to 2729 3650 to 3659 4570 to 4579 4870 to 4879 6050 to 6059 6700 to 6709 8060 to 8069 8080 to 8089 9170 to 9179 10420 to 10429 11770 to 11779 12170 to 12179 12520 to 12529 14640 to 14649 16390 to 16399 18020 to 18029
mol material AOEdgy
mol addrep 0

## Adding new representation
#mol color Name
#mol representation DynamicBonds 1.600000 0.300000 12.000000
#mol selection index 1220 to 1229 740 to 749 11800 to 11809 18540 to 18549 18970 to 18979 20710 to 20719 17440 to 17449 10130 to 10139 4970 to 4979 6390 to 6399 18470 to 18479 1860 to 1869 2130 to 2139 7540 to 7549 6690 to 6699 6470 to 6479 4150 to 4159
#mol material AOEdgy
#mol addrep 0
#
#mol color Name
#mol representation DynamicBonds 1.600000 0.300000 12.000000
## mol representation VDW 0.300000 14.000000
#mol selection index index 5840 to 5849 820 to 829 8050 to 8059 2830 to 2839 620 to 629 1340 to 1349 10160 to 10169 8770 to 8779 17680 to 17689
#mol material AOEdgy
#mol addrep 0
#
#
#mol color Name
#mol representation DynamicBonds 1.600000 0.300000 12.000000
## mol representation VDW 0.300000 14.000000
#mol selection index index 14740 to 14749 15300 to 15309 17480 to 17489 14420 to 14429 9870 to 9879 17330 to 17339 8250 to 8259 20680 to 20689 15040 to 15049 3980 to 3989
#mol material AOEdgy
#mol addrep 0

# Display options
display depthcue on
display projection Orthographic
display shadows on
display ambientocclusion on

pbc box off

# Show only %i and %a in labels
# label add Atoms 0/12093
# label textformat Atoms 1 { %i(%a)  }

# For show
# label show Atoms all

# Or for hide
# label hide Atoms all

# Selecting a polimer
# mol modselect 1 0 index 1220 to 1229 740 to 749 11800 to 11809 18540 to 18549 18970 to 18979 20710 to 20719 17440 to 17449 10130 to 10139 4970 to 4979 6390 to 6399 18470 to 18479 1860 to 1869 2130 to 2139 7540 to 7549 6690 to 6699 6470 to 6479 4150 to 4159