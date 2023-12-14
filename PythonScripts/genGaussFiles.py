"""Script dedicated to the generation of gaussian files from an xyz structure."""

from molcraft.structure import load_xyz
import sys


def GaussInputFile_method2(name, coord, title="", print_info=False):
    lines = ""
    lines += "%nprocshared=20\n"
    lines += "%mem=40GB\n"
    lines += "%chk={}.chk\n".format(name)
    lines += "#p oniom(td=(nstates=6) cam-b3lyp/6-311+g(d,p):hf/6-31g(d))=ptembed\n"
    lines += "\n"
    lines += "{}\n".format(title)
    lines += "\n"
    lines += "0 1 0 1 0 1\n"
    
    for i in coord.index:
        lines += " %s%17d%14.8f%14.8f%14.8f%2s\n" % (
            coord.loc[i, "atsb"],
            0,
            coord.loc[i, "x"],
            coord.loc[i, "y"],
            coord.loc[i, "z"],
            coord.loc[i, "layer"]
        )
        
    lines += "\n"
        
    if print_info:
        print(lines)
    else:
        with open(f"{name}.com", "w") as INPUT:
            INPUT.write(lines)


def GaussInputFile_method3(name, coord, title="", print_info=False):
    lines = ""
    lines += "%nprocshared=12\n"
    lines += "%mem=2GB\n"
    lines += "%chk={}.chk\n".format(name)
    lines += "#p td=(nstates=6) cam-b3lyp/6-311+g(d,p)\n"
    lines += "\n"
    lines += "{}\n".format(title)
    lines += "\n"
    lines += "0 1\n"
    
    for i in coord.index:
        lines += " %s%14.8f%14.8f%14.8f\n" % (
            coord.loc[i, "atsb"],
            coord.loc[i, "x"],
            coord.loc[i, "y"],
            coord.loc[i, "z"]
        )
        
    lines += "\n"
        
    if print_info:
        print(lines)
    else:
        with open(f"{name}.com", "w") as INPUT:
            INPUT.write(lines)


def GaussInputFile_method4(name, coord, title="", print_info=False):
    lines = ""
    lines += "%nprocshared=12\n"
    lines += "%mem=2GB\n"
    lines += "%chk={}.chk\n".format(name)
    lines += "#p td=(nstates=6) cam-b3lyp/6-311+g(d,p) scrf=(Read)\n"
    lines += "\n"
    lines += "{}\n".format(title)
    lines += "\n"
    lines += "0 1\n"
    
    for i in coord.index:
        lines += " %s%14.8f%14.8f%14.8f\n" % (
            coord.loc[i, "atsb"],
            coord.loc[i, "x"],
            coord.loc[i, "y"],
            coord.loc[i, "z"]
        )
        
    lines += "\nEps=3.500\n\n"
        
    if print_info:
        print(lines)
    else:
        with open(f"{name}.com", "w") as INPUT:
            INPUT.write(lines)
    

try:
    file = sys.argv[1]
    iso = sys.argv[2]
except IndexError:
    print("No file has been selected")
    print("\tUsage: python genGaussFiles.py file.xyz iso")
    exit()

print("Selected file:", file)


azo_pol = load_xyz(file)
azo_pol["layer"] = "L"
azo_pol.loc[range(26), "layer"] = "H"
print(azo_pol)

azo = azo_pol.loc[0:25, :]

GaussInputFile_method2("azoO" + iso[0].upper() + "_method2", azo_pol, title="azoO PBOH - Method 2 - QM/QM", print_info=False)
GaussInputFile_method3("azoO" + iso[0].upper() + "_method3", azo, title="azoO PBOH - Method 3 - QM", print_info=False)
GaussInputFile_method4("azoO" + iso[0].upper() + "_method4", azo, title="azoO PBOH - Method 4 - QM SCRF", print_info=False)
