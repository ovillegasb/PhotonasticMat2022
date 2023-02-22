
pbc box -off

# mol modstyle 0 0 Lines 1.000000
mol modstyle 0 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modmaterial 0 0 Transparent

set atoms [atomselect 0 "all"]
set coords [$atoms get {x y z}]
set names [$atoms get name]
set natoms [$atoms num]

# move by
set mv_to {0.0 0.25 0.0}

# Read charges
proc read_Charges {file} {
    set fp [open $file r]
    set file_data [read $fp]
    close $fp

    #  Process data file
    set data [split $file_data "\n"]

    return $data
}

# draw color black
graphics top color black

set charges [read_Charges "charges.dat"]
puts $charges

for {set i 0} {$i < $natoms} {incr i} {
    set at [lindex $coords $i]
    graphics top text [vecadd $mv_to $at] "[lindex $charges $i]" size 1.2 thickness 2
}

#graphics top delete all
