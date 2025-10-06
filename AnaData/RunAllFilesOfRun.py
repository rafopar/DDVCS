import os
import glob


index = 0;
for name in glob.glob('Data_DST/006638/rec_clas_*'):
    print(name)

    cmd = "./AnaData.exe %s %d" %( name, index )
    
    os.system(cmd);
    
    index = index + 1
