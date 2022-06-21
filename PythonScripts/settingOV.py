"""
Module used to configure a jupyter session.

Records output folders or useful functions.

"""

import platform

mySystem = platform.uname()

folderImg = ""
if mySystem.node == "freya":
    # LAPTOP office work
    folderImg += "/home/ovillegas/Pictures"
else:
    print("Unrecognized system")
    print("-------------------")
    print(f"System: {mySystem.system}")
    print(f"Node Name: {mySystem.node}")
    print(f"Release: {mySystem.release}")
    print(f"Version: {mySystem.version}")
    print(f"Machine: {mySystem.machine}")
    print(f"Processor: {mySystem.processor}")
