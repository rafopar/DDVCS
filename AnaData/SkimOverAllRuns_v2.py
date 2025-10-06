import os
import glob
import time
import subprocess

processes = set()

index = 0;
#for name in glob.glob('Data_F18In/jpsitcs_00*.hipo'):
#for name in glob.glob('Data_F18Out/jpsitcs_00*.hipo'):
for name in glob.glob('Data_F18InEarly/jpsitcs_00*.hipo'):
    print(name)

    name_parts = name.split("_");
    name_parts2 = name_parts[2].split(".")
    run = int(name_parts2[0])
    
    cmd = "./AnaData.exe %d" %( run )
    print("Processing the run %d" %(run))

    subprocess.Popen([cmd], shell = True)
    #os.system(cmd);
    
    index = index + 1

    if (index % 5) == 0:
        print("Index is %d Waiting for the 2nd bunch of runss"%(index))
        time.sleep(400)
