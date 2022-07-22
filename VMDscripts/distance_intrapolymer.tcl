

# # Material of bone background
# material add Mat copy Transparent
# material change ambient Mat 0.000000
# material change diffuse Mat 0.650000
# material change specular Mat 0.500000
# material change shininess Mat 0.530000
# material change mirror Mat 0.000000
# material change opacity Mat 0.080000
# material change outline Mat 0.000000
# material change outlinewidth Mat 0.000000

# Data monomro 1
mol color ColorID 0
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol selection index 120 to 129
mol material AOEdgy
mol addrep 0

# select coordinates
# # or 0 --> top
set mono1 [atomselect 0 "index 120 to 129"]
set cmono1 [measure center $mono1]
graphics top color green
graphics top sphere $cmono1 radius 0.3

# Data monomro 2
mol color ColorID 1
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol selection index 30 to 39
mol material AOEdgy
mol addrep 0

# select coordinates
set mono2 [atomselect 0 "index 30 to 39"]
set cmono2 [measure center $mono2]
graphics top sphere $cmono2 radius 0.3


# Distance between both monomers
set d12 [ veclength [ vecsub $cmono2 $cmono1 ] ]

# Data monomer 3
mol color ColorID 1
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol selection index 290 to 299
mol material AOEdgy
mol addrep 0

# select coordinates
set mono3 [atomselect 0 "index 290 to 299"]
set cmono3 [measure center $mono3]
graphics top sphere $cmono3 radius 0.3

# Distance between both monomers
set d13 [ veclength [ vecsub $cmono3 $cmono1 ] ]

# Data monomer 4
mol color ColorID 3
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol selection index 220 to 229
mol material Transparent
mol addrep 0

# select coordinates
set mono4 [atomselect 0 "index 220 to 229"]
set cmono4 [measure center $mono4]
graphics top sphere $cmono4 radius 0.3

# Distance between both monomers
set d14 [ veclength [ vecsub $cmono4 $cmono1 ] ]

# Data monomer 5
mol color ColorID 3
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol selection index 310 to 319
mol material Transparent
mol addrep 0

# select coordinates
set mono5 [atomselect 0 "index 310 to 319"]
set cmono5 [measure center $mono5]
graphics top sphere $cmono5 radius 0.3

# Distance between both monomers
set d15 [ veclength [ vecsub $cmono5 $cmono1 ] ]

# Data monomer 6
mol color ColorID 2
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol selection index 50 to 59
mol material Transparent
mol addrep 0

# select coordinates
set mono6 [atomselect 0 "index 50 to 59"]
set cmono6 [measure center $mono6]
graphics top sphere $cmono6 radius 0.3

# Distance between both monomers
set d16 [ veclength [ vecsub $cmono6 $cmono1 ] ]


# For add lines
graphics top color black

graphics top line $cmono1 $cmono3 width 3 style dashed
graphics top line $cmono1 $cmono5 width 3 style dashed
graphics top line $cmono1 $cmono6 width 3 style dashed

# info
graphics top list
graphics top info 42
