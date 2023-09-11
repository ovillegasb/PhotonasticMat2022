
import os
import glob
import sys


folder = sys.argv[1].replace("/", "")
print("Input folder:", folder)

files = glob.glob(f"{folder}/*")
print("Number of files:", len(files))

folder_out = folder + "_notSolvent"
print("Output folder:", folder_out)

if not os.path.exists(folder_out):
    os.mkdir(folder_out)

for f in files:
    lines = ""
    f_out = f.split("/")[-1]
    with open(f, "r") as INP:
        for line in INP:
            if line.startswith("#"):
                line = " ".join([a for a in line.split() if "SCRF" not in a])

            lines += line


    with open(f"{folder_out}/{f_out}", "w") as OUT:
        OUT.write(lines)

print("Number of files processed:", len(glob.glob(f"{folder_out}/*")))
print("Finish!")
