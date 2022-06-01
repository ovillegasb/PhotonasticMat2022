
# Script recording the evolution of the reactivity of a
# system during polymerization in STAMP.

# RUN in folder ./XYZ
# vmd -dispdev text -e path/save_polymerInfo.tcl

# list all xyz
set xyzlist [glob *.xyz]

# numbers of files
set Nxyz [llength $xyzlist]
puts "XYZ files number: $Nxyz"

# Freq save xyz
set step 100


set ofile [open "polystate.dat" w]
puts $ofile "# i_xyz graine donneur acceptor_l"

for {set i 0} {$i < $Nxyz} {incr i} {

    # Select from list
    set ifile [lindex $xyzlist $i]

    # Load new molecule
    puts "File: $ifile"
    mol new $ifile

    # Graine : name Br
    # for show only
    # mol modselect 0 0 name Br
    set grain [atomselect top "name Br"]
    puts "N Graine [$grain num]"

    # Donneur : name Kr
    # for show only
    # mol modselect 0 0 name Kr
    set donor [atomselect top "name Kr"]
    puts "N Donneur [$donor num]"
    
    # AccepteurLibre : name Xe
    # for show only
    # mol modselect 0 0 name Xe
    set acceptorL [atomselect top "name Xe"]
    puts "N AccepteurLibre [$acceptorL num]"


    puts $ofile "[expr $i * $step] [$grain num] [$donor num] [$acceptorL num]"
    mol delete all
}

close $ofile
exit