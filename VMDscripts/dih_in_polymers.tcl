

# Style
mol modstyle 0 0 CPK 1.000000 0.300000 52.000000 12.000000
mol modmaterial 0 0 Transparent
mol showrep 0 0 0

mol color Name
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol selection index 0 to 25
mol material AOEdgy
mol addrep 0

set atoms [atomselect 0 "index 11 12 13 14"]
$atoms list

# get coordinates
set coords [$atoms get {x y z}]

# length
puts [llength $coords]

# accese to per index
# lindex listname index

set C1 [lindex $coords 0]
set N1 [lindex $coords 1]
set N2 [lindex $coords 2]
set C2 [lindex $coords 3]


# draw color black
##graphics top color black

# draw text [vecadd {-0.1 0.5 0.} $C5] "C5"
##graphics top text [vecadd {-0.1 0.5 0.} $C1] "C" size 1.2 thickness 2
#draw text [vecadd {-0.1 0.5 0.} $N1] "N1"
##graphics top text [vecadd {-0.1 0.5 0.} $N1] "N" size 1.2 thickness 2
#draw text [vecadd {-0.1 0.5 0.} $N2] "N2"
##graphics top text [vecadd {-0.1 0.5 0.} $N2] "N" size 1.2 thickness 2
#draw text [vecadd {-0.1 0.5 0.} $C7] "C7"
##graphics top text [vecadd {-0.1 0.5 0.} $C2] "C" size 1.2 thickness 2

mol addrep 0
mol modstyle 1 0 Bonds 0.100000 12.000000
mol modselect 1 0 "index 11 12 13 14"
mol modcolor 1 0 ColorID 25

label add Dihedrals 0/11 0/12 0/13 0/14

