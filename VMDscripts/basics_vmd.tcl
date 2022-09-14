# Documentation
# http://www.ks.uiuc.edu/Research/vmd/current/docs.html
# https://www.ks.uiuc.edu/Research/vmd/current/ug/node117.html
# https://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/node4.html
# https://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node182.html

# The Basics of Tcl Scripting

# Print information
puts "Hello world"

# Equations
expr -3 * 10

# Variable definition
set x 3

# Variable from equation
set x [expr -3 * 10]

# print a variable
puts $x
puts "the value of x is: $x"

# Equation from variable
expr -3 * $x

# [expr.] - represents the results of the expression inside the brackets

# for {initialization} {test} {increment} {commands}

# Calculate the values of -3*x for integers x from 0 to 10 and output the results into a
# file named myoutput.dat.

set file [open "out.dat" w]
for {set x 0} {$x <= 10} {incr x} {
	puts $file [ expr -3 * $x ]
}

# Load new molecule
mol new structure.pdb


# Style mol
# VMW:               Scale    Resolution
mol modstyle 0 0 VDW 1.000000 12.000000

# Color mol
mol modcolor 0 0 ResType

# Delete mol representation
mol delrep 0 top

# Adding representation
mol representation NewCartoon 0.3 10.0 4.1 0
mol color Structure
mol selection {protein}
mol addrep top

# open menus
menu graphics on
menu graphics off

# Global parameters
display projection orthographic
color Display Background white
axes location off

# UNIX: $HOME/.vmdrc
# create a Tcl script that will be "played" when VMD starts.

# Display options
# Show main window
menu main on
menu tkcon on

# Kayboard shortcuts
user add key a "axes location off"
user add key A "axes location lowerleft"

# Colour definitions
# blue 
color change rgb 0 0.00 0.47 0.72

# Selecting atoms
# creates a new atom selection
# atomselect molid selection

set Hatoms [atomselect 0 "name H"] # or 0 --> top
$Hatoms list
$Hatoms keywords
$Hatoms moveby {10 0 0}
$Hatoms move [transaxis x 40 deg]

$Hatoms get resid
$Hatoms get resname
$Hatoms get {resname resid}
$Hatoms get {x y z}
$Hatoms get mass

#  If you want to obtain some of the structural properties
measure center $Hatoms

# Determine the symmetry
set result [measure symmetry $Hatoms]
# Create array ’symm’ containing the results
array set symm $result
# Print selected elements of the array
puts $symm(pointgroup)
puts $symm(order)
puts $symm(elements)
puts $symm(axes)
# you can use commands to learn about the properties
# (number of atoms, coordinates, total charge, etc) 

# Once you are done with a selection, it is always a
# good ideal to delete it to save memory
$symm delete

# Labels
label add Atoms 0/0
label add Atoms 0/2
label add Bonds 0/0 0/2

# Drawing shapes

# VMD offers a way to display user-defined objects built
# from graphics primitives such as points, lines, cylinders,
# cones, spheres, triangles, and text. The command that can
# realize those functions is graphics, the syntax of which is

# graphics molid command

graphics top point {0 0 10}
graphics top line {-10 0 0} {0 0 0} width 5 style solid
graphics top line {10 0 0} {0 0 0} width 5 style dashed
graphics top cylinder {15 0 0} {15 0 10} radius 10 resolution 60 filled no
graphics top text {0 0 1} "Hola"
graphics top list
graphics top info 0
graphics top delete ID

# Changing a variable from atomselect
set atoms [atomselect 0 "all"]
set q "-1.2 -1.2 0.4 0.4 0.4 0.4 0.4 0.4"
$atoms set charge $q

# Changing coordinates for an atom
$atom set {x y z} {{1.6 0 0}}

# Change color resname
color Resname A14 tan
color change rgb 5 tan 0.390000 0.500000 0.200000