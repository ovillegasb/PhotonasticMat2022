
package require hbonds


set current_directory [pwd]
regexp {/(cis|trans)/(\d+)} $current_directory match type number

# # Lee la entrada del usuario y almacénala en la variable input
# gets stdin input

hbonds -sel1 [atomselect 0 "resid 1 to $number"] -sel2  [atomselect 0 "not resid 1 to $number"] -writefile yes -upsel yes -dist 3.0 -ang 20 -plot no -outfile "hbonds_all.dat" -frames 1:-1
# OH -- O PC
# Observar si agrego o no H
##hbonds -sel1 [atomselect 0 "resid 1 to $number and name O"] -sel2  [atomselect 0 "not resid 1 to $number"] -writefile yes -upsel yes -dist 3.0 -ang 20 -plot no -outfile "hbonds_to_O_pc.dat" -frames 1:-1
# PC OH -- O PC?
##hbonds -sel1 [atomselect 0 "resid 1 to $number"] -sel2  [atomselect 0 "not resid 1 to $number  and name O"] -writefile yes -upsel yes -dist 3.0 -ang 20 -plot no -outfile "hbonds_to_O_pol.dat" -frames 1:-1
# OH -- N PC
##hbonds -sel1 [atomselect 0 "resid 1 to $number and name N"] -sel2  [atomselect 0 "not resid 1 to $number"] -writefile yes -upsel yes -dist 3.0 -ang 20 -plot no -outfile "hbonds_to_N.dat" -frames 1:-1

puts "Finish!"
exit

## set hbonds [measure hbonds 3.0 20 [atomselect 0 "resid 1 to 50"] [atomselect 0 "not resid 1 to 50"]]
## set donnors [lindex $hbonds 0]
## set acceptors [lindex $hbonds 1]
## set hydrogens [lindex $hbonds 2]