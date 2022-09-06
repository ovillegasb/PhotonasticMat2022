#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Adapt the FAtomes file.

This script is made to add two things:

1) activate the photochrome parameters,
2) add the new parameters of the latest Stamp version.
"""

import sys
from datetime import date

try:
    file = sys.argv[1]
except IndexError:
    print("File not selected")
    print("Usage : UpdateFAtomes.py FAtomes.in\n")
    exit()

print("Entry file:", file)

print("Indicate the conformer:")
iso = input("cis or trans: ")

polymerCoor = ""
polymerConnectivity = ""

with open(file, "r") as FATOME:
    for line in FATOME:
        if "PositionDesAtomesCart" in line:
            N = int(FATOME.readline())
            print("N PositionDesAtomesCart:", N)
            for i in range(N):
                coord = FATOME.readline()
                if i >= 26:
                    polymerCoor += coord

        elif "Zmatrice" in line:
            N = int(FATOME.readline())
            print("N Connectivity:", N)
            for i in range(N):
                conn = FATOME.readline()
                if i >= 26:
                    polymerConnectivity += conn

out = "FAtomes.in"

lines = ""
lines += """*#####################################################################*
*                       File FAtomes Stamp4 (22.07.21)                *
*                       Created on : ({})                     *
*#####################################################################*
""".format(date.today())

lines += """* #####################
* ## Types atomiques ##
* #####################
NbTypesAtomes 11
*
* Atome 0
nom             ca
nomXYZ          C
nomFF           ca
type            Atome
masse           12.e-3 kg/mol
structure       FICHIER
maille_long     65.00200000 60.00000000 60.00000000 ang
maille_angle    90.0 90.0 90.0 degre
maille_orient   0                               ** 0={a~x / b~xy} - 1={a~z / b~yz} - 2={a~z / c~yz} - 3={b~z / c~yz}
maille_ref
* Atome 1
nom             ha
nomXYZ          H
nomFF           ha
type            Atome
masse           1.e-3  kg/mol
* Atome 2
nom             oh
nomXYZ          O
nomFF           oh
type            Atome
masse           16.e-3 kg/mol
* Atome 3
nom             ho
nomXYZ          H
nomFF           ho
type            Atome
masse           1.e-3  kg/mol
* Atome 4
nom             ne
nomXYZ          N
nomFF           ne
type            Atome
masse           14.e-3 kg/mol
* Atome 5
nom             nf
nomXYZ          N
nomFF           nf
type            Atome
masse           14.e-3 kg/mol
*  Atome 6 = 0bis
nom             cb 
nomXYZ          C
nomFF           ca
type            Atome
masse           12.e-3 kg/mol
* Atome BTDN_1
nom             CT
nomXYZ          C
nomFF           CT
type            Atome
masse           12.011e-03 kg/mol
charge          -0.12 e-
* Atome BTDN_2
nom             CM
nomXYZ          C
nomFF           CM
type            Atome
masse           12.011e-03 kg/mol
charge          -0.115 e-
* Atome BTDN_3
nom             HT
nomXYZ          H
nomFF           HT
type            Atome
masse           1.008e-03 kg/mol
charge          0.06 e-
* Atome BTDN_4
nom             HM
nomXYZ          H
nomFF           HM
type            Atome
masse           1.008e-03 kg/mol
charge          0.115 e-
"""

lines += """* ################################
* ## Potentiel intermoleculaire ##
* ################################
Potentiel  0  0 LJ sigma   3.3996695084e+00 ang epsilon   3.7318093892e-03 eV rc   2.5000000000e+00 -
Potentiel  1  1 LJ sigma   2.5996424595e+00 ang epsilon   6.5089698650e-04 eV rc   2.5000000000e+00 -
Potentiel  2  2 LJ sigma   3.0664733878e+00 ang epsilon   9.1299150639e-03 eV rc   2.5000000000e+00 -
Potentiel  3  3 LJ sigma   0.0000000000e+00 ang epsilon   0.0000000000e+00 eV rc   2.5000000000e+00 -
Potentiel  4  4 LJ sigma   3.2499985238e+00 ang epsilon   7.3768325136e-03 eV rc   2.5000000000e+00 -
Potentiel  5  5 LJ sigma   3.2499985238e+00 ang epsilon   7.3768325136e-03 eV rc   2.5000000000e+00 -
Potentiel  6  6 LJ sigma   3.3996695084e+00 ang epsilon   3.7318093892e-03 eV rc   2.5000000000e+00 -
*
Potentiel  7  7   LJ sigma 3.5    ang epsilon 0.066  kcal/mol rc 12.5 ang
Potentiel  7  8   LJ sigma 3.525  ang epsilon 0.0708 kcal/mol rc 12.5 ang
Potentiel  7  9   LJ sigma 2.958  ang epsilon 0.0445 kcal/mol rc 12.5 ang
Potentiel  7  10  LJ sigma 2.91   ang epsilon 0.0445 kcal/mol rc 12.5 ang
Potentiel  8  8   LJ sigma 3.55   ang epsilon 0.076  kcal/mol rc 12.5 ang
Potentiel  8  9   LJ sigma 2.979  ang epsilon 0.0477 kcal/mol rc 12.5 ang
Potentiel  8  10  LJ sigma 2.931  ang epsilon 0.0477 kcal/mol rc 12.5 ang
Potentiel  9  9   LJ sigma 2.5    ang epsilon 0.030  kcal/mol rc 12.5 ang
Potentiel  9  10  LJ sigma 2.46   ang epsilon 0.030  kcal/mol rc 12.5 ang
Potentiel  10 10  LJ sigma 2.42   ang epsilon 0.030  kcal/mol rc 12.5 ang
Regle_melange Lorentz-Berthelot
"""

lines += """* #####################
* ## Champ de forces ##
* #####################
ChampDeForces
55
bond_gaff ca ca       1.3870000000 ang     478.4000000000 kcal/mol/ang2
bond_gaff ca ha       1.0870000000 ang     344.3000000000 kcal/mol/ang2
bond_gaff ca oh       1.3620000000 ang     386.1000000000 kcal/mol/ang2
bond_gaff ca ne       1.4310000000 ang     361.8000000000 kcal/mol/ang2
bond_gaff ca nf       1.4310000000 ang     361.8000000000 kcal/mol/ang2
bond_gaff oh ho       0.9740000000 ang     369.6000000000 kcal/mol/ang2
bond_gaff ne nf       1.2475 ang     797.8847 kcal/mol/ang2
angle_gaff ca ca ca     119.9700000000 deg      67.1800000000 kcal/mol
angle_gaff ca ca ha     120.0100000000 deg      48.4600000000 kcal/mol
angle_gaff ca ca oh     119.9400000000 deg      69.8500000000 kcal/mol
angle_gaff ca ca ne     103.4 deg      53.2182 kcal/mol
angle_gaff ca ca nf     103.4 deg      53.2182 kcal/mol
angle_gaff ca oh ho     109.4700000000 deg      48.8500000000 kcal/mol
angle_gaff ca ne nf     110.1 deg      76.4997 kcal/mol
angle_gaff ca nf ne     110.1 deg      76.4997 kcal/mol
torsion_gaff ca ca ca ca      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg
torsion_gaff ca ca ca ha      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg
torsion_gaff ca ca ca oh      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg
torsion_gaff ca ca ca ne      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg
torsion_gaff ca ca ca nf      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg
torsion_gaff ca ca oh ho       1.8000000000 kcal/mol   2.00 -   2.00 -     180.0000000000 deg
torsion_gaff ca ne nf ca       17.6364 kcal/mol   1.00 -  2.00 -     180.0 deg
torsion_gaff ca ne nf ca        4.7014 kcal/mol   1.00 -  1.00 -       0.0 deg
torsion_gaff ha ca ca ha      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg
torsion_gaff ha ca ca oh      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg
torsion_gaff ha ca ca ne      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg
torsion_gaff ha ca ca nf      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg
torsion_gaff ca ca nf ne       2.39005736   kcal/mol   1.00 -   2.00 -     180.0000000000 deg
torsion_gaff ca ca ne nf       2.39005736   kcal/mol   1.00 -   2.00 -     180.0000000000 deg
impropre_gaff    ca ca ca ha          1.1000000000 kcal/mol
impropre_gaff    ca ca ca oh          1.1000000000 kcal/mol
impropre_gaff    ca ca ca ne          1.1000000000 kcal/mol
impropre_gaff    ca ca ca nf          1.1000000000 kcal/mol
*
bond_opls              CT      CT                      1.529 ang       268.0 kcal/mol/ang2     0.0 -   0.0 -
bond_opls              CT      HT                      1.090 ang       340.0 kcal/mol/ang2     0.0 -   0.0 -
bond_opls              CT      CM                      1.510 ang       317.0 kcal/mol/ang2     0.0 -   0.0 -
bond_opls              CM      HM                      1.080 ang       340.0 kcal/mol/ang2     0.0 -   0.0 -
bond_opls              CM      CM                      1.340 ang       549.0 kcal/mol/ang2     0.0 -   0.0 -
angle_opls             HT      CT      HT              107.8 degre     33.0  kcal/mol          0.0 -   0.0 -
angle_opls             CT      CT      HT              110.7 degre     37.5  kcal/mol          0.0 -   0.0 -
angle_opls             HT      CT      CM              109.5 degre     35.0  kcal/mol          0.0 -   0.0 -
angle_opls             CT      CM      HM              117.0 degre     35.0  kcal/mol          0.0 -   0.0 -
angle_opls             CT      CM      CM              124.0 degre     70.0  kcal/mol          0.0 -   0.0 -
angle_opls             HM      CM      CM              120.0 degre     35.0  kcal/mol          0.0 -   0.0 -
angle_opls             CT      CT      CM              111.1 degre     63.0  kcal/mol          0.0 -   0.0 -
torsion_opls           HT      CT      CT      HT      0.000  kcal/mol 0.000   kcal/mol         0.300  kcal/mol
torsion_opls           HT      CT      CT      CM      0.000  kcal/mol 0.000   kcal/mol         0.366  kcal/mol
torsion_opls           CT      CT      CM      HM      0.000  kcal/mol 0.000   kcal/mol         0.000  kcal/mol
torsion_opls           HT      CT      CM      HM      0.000  kcal/mol 0.000   kcal/mol         0.318  kcal/mol
torsion_opls           CM      CM      CT      CT      0.346  kcal/mol 0.405   kcal/mol        -0.940  kcal/mol
torsion_opls           CT      CM      CM      CT      0.000  kcal/mol 1400.00 kcal/mol         0.000  kcal/mol
torsion_opls           HM      CM      CM      HM      0.000  kcal/mol 1400.00 kcal/mol         0.000  kcal/mol
torsion_opls           CT      CM      CM      HM      0.000  kcal/mol 1400.00 kcal/mol         0.000  kcal/mol
torsion_opls           HT      CT      CM      CM      0.000  kcal/mol 0.000   kcal/mol        -0.372  kcal/mol
torsion_opls           CM      CT      CT      CM      1.300  kcal/mol -0.05   kcal/mol         0.200  kcal/mol
"""

IsoCoor = {
"cis": """* #########################
* ## Positions atomiques ##
* #########################
PositionDesAtomesCart angstrom
21626
  ca       2.1892260000       4.9700910000       5.6843710000
  ha       1.4830670000       5.7595370000       5.4729160000  
  ca       1.7376040000       3.7811770000       6.2099020000
  ha       0.6896710000       3.6252110000       6.4175470000
  ca       2.6313610000       2.7513100000       6.4849510000
  oh       2.1388090000       1.6206740000       7.0139280000
  ho       2.8460740000       0.9908180000       7.1435740000
  ca       3.9782580000       2.9192500000       6.1894930000
  ha       4.6774020000       2.1134560000       6.3808940000
  ca       4.4193880000       4.0977250000       5.6247200000
  ha       5.4577330000       4.2236580000       5.3519570000
  cb       3.5431600000       5.1522680000       5.4005570000
  ne       4.0641630000       6.2734010000       4.7218930000
  nf       3.7535430000       7.4464410000       4.9490330000
  cb       3.0101660000       7.8853630000       6.0642480000
  ca       3.2412890000       7.4616980000       7.3734380000
  ha       3.9594330000       6.6787010000       7.5678870000
  ca       2.5804470000       8.0504000000       8.4262600000
  ha       2.7650410000       7.7401300000       9.4452270000
  ca       1.6545680000       9.0629750000       8.1967740000
  oh       1.0254040000       9.5883750000       9.2592600000
  ho       0.4292200000      10.2781280000       8.9720190000
  ca       1.4315680000       9.5014190000       6.8976740000
  ha       0.7282590000      10.3039090000       6.7122240000
  ca       2.1265970000       8.9363590000       5.8487770000
  ha       1.9929700000       9.2981550000       4.8396740000
""",
"trans":"""* #########################
* ## Positions atomiques ##
* #########################
PositionDesAtomesCart angstrom
21626
  ca       0.0010000000       4.0370000000       3.3910000000
  ha       0.0010000000       4.6900000000       4.2540000000
  ca       0.0010000000       2.6620000000       3.5480000000
  ha       0.0010000000       2.2090000000       4.5330000000
  ca       0.0010000000       1.8280000000       2.4220000000
  oh       0.0020000000       0.4810000000       2.6470000000
  ho       0.0010000000       0.0000000000       1.8120000000
  ca       0.0010000000       2.3800000000       1.1400000000
  ha       0.0010000000       1.7350000000       0.2650000000
  ca       0.0010000000       3.7630000000       0.9880000000
  ha       0.0000000000       4.2100000000       0.0000000000
  cb       0.0010000000       4.6050000000       2.1040000000
  ne       0.0010000000       5.9910000000       1.8310000000
  nf       0.0010000000       6.7420000000       2.8360000000
  cb       0.0010000000       8.1270000000       2.5630000000
  ca       0.0010000000       8.6950000000       1.2760000000
  ha       0.0010000000       8.0420000000       0.4130000000
  ca       0.0010000000      10.0700000000       1.1190000000
  ha       0.0010000000      10.5230000000       0.1340000000
  ca       0.0010000000      10.9040000000       2.2450000000
  oh       0.0020000000      12.2520000000       2.0190000000
  ho       0.0020000000      12.7320000000       2.8550000000
  ca       0.0010000000      10.3520000000       3.5270000000
  ha       0.0010000000      10.9970000000       4.4020000000
  ca       0.0010000000       8.9700000000       3.6790000000
  ha       0.0000000000       8.5220000000       4.6660000000
"""}

lines += IsoCoor[iso]

lines += polymerCoor

lines += """* ###############
* ## Z-matrice ##
* ###############
Zmatrice
21626
 0         2         11          1 
 1         0 
 2         4          0          3 
 3         2 
 4         5          7          2 
 5         6          4 
 6         5 
 7         4          8          9 
 8         7 
 9         7         10         11 
10         9 
11         9          0         12 
12        11         13 
13        12         14 
14        13         15         24 
15        14         16         17 
16        15 
17        15         18         19 
18        17 
19        17         22         20 
20        19         21 
21        20 
22        24         23         19 
23        22 
24        14         22         25 
25        24 
"""

lines += polymerConnectivity

lines += """* #######################
* ## Charges atomiques ##
* #######################
ModificationChargeDesAtomes e-
26
0  -9.66410768e-02
1   1.09659923e-01
2  -3.11919077e-01
3   1.95152923e-01
4   3.36819923e-01
5  -5.22373077e-01
6   3.76906923e-01
7  -3.11919077e-01
8   1.95152923e-01
9  -9.66410768e-02
10   1.09659923e-01
11   1.51418923e-01
12  -1.35278077e-01
13  -1.35278077e-01
14   1.51418923e-01
15  -9.66410768e-02
16   1.09659923e-01
17  -3.11919077e-01
18   1.95152923e-01
19   3.36819923e-01
20  -5.22373077e-01
21   3.76906923e-01
22  -3.11919077e-01
23   1.95152923e-01
24  -9.66410768e-02
25   1.09659923e-01
*
* ==============================================================================
*           Contribution de dispersion repulsion en intramoleculaire
* ==============================================================================
ContribDispRepIntra
28
cb cb 0.0 0.0 0.0 0.0 0.0 0.0
ca ca 0.0 0.0 0.50 0.0 0.0 0.8333
ca ha 0.0 0.0 0.50 0.0 0.0 0.8333
ca oh 0.0 0.0 0.50 0.0 0.0 0.8333
ca ho 0.0 0.0 0.50 0.0 0.0 0.8333
ca ne 0.0 0.0 0.50 0.0 0.0 0.8333
ca nf 0.0 0.0 0.50 0.0 0.0 0.8333
ca cb 0.0 0.0 0.50 0.0 0.0 0.8333
ha ha 0.0 0.0 0.50 0.0 0.0 0.8333
ha oh 0.0 0.0 0.50 0.0 0.0 0.8333
ha ho 0.0 0.0 0.50 0.0 0.0 0.8333
ha ne 0.0 0.0 0.50 0.0 0.0 0.8333
ha nf 0.0 0.0 0.50 0.0 0.0 0.8333
ha cb 0.0 0.0 0.50 0.0 0.0 0.8333
oh oh 0.0 0.0 0.50 0.0 0.0 0.8333
oh ho 0.0 0.0 0.50 0.0 0.0 0.8333
oh ne 0.0 0.0 0.50 0.0 0.0 0.8333
oh nf 0.0 0.0 0.50 0.0 0.0 0.8333
oh cb 0.0 0.0 0.50 0.0 0.0 0.8333
ho ho 0.0 0.0 0.50 0.0 0.0 0.8333
ho ne 0.0 0.0 0.50 0.0 0.0 0.8333
ho nf 0.0 0.0 0.50 0.0 0.0 0.8333
ho cb 0.0 0.0 0.50 0.0 0.0 0.8333
ne ne 0.0 0.0 0.50 0.0 0.0 0.8333
ne nf 0.0 0.0 0.50 0.0 0.0 0.8333
ne cb 0.0 0.0 0.50 0.0 0.0 0.8333
nf nf 0.0 0.0 0.50 0.0 0.0 0.8333
nf cb 0.0 0.0 0.50 0.0 0.0 0.8333
"""

with open(out, "w") as OUT:
    OUT.write(lines)
