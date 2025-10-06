import os
import glob


index = 0;
for name in glob.glob('Data/jpsitcs_00*.hipo'):
    print(name)

    name_parts = name.split("_");
    name_parts2 = name_parts[1].split(".")
    run = int(name_parts2[0])
    
    cmd = "./AnaData.exe %d" %( run )
    print("Processing the run %d" %(run))
    os.system(cmd);
    
    index = index + 1
