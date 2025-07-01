import subprocess
import shlex
import os
import time
import sys

def check_NumberofProcesses(proc_dict):
    nRunningProcess = 0
    for ind in proc_dict:
        if proc_dict[ind].poll() is None:
            nRunningProcess = nRunningProcess + 1
    return nRunningProcess;


runs = [5117,5125,5126,5128,5130,5163,5165,5166,5167,5169,5181,5191,5196,5197,5198,5199,5202,5203,5204,5212,5215,5219,5220,5221,5222,5230,5231,5232,5233,5234,5247,5248,5249,5257,5258,5261,5303,5304,5306,5317,5318,5319,5346,5356,5357,5358,5359,5360,5366,5367,5374,5375,5379,5381,5391,5393,5407]

print(runs)

run_counter = 0
proc_ana = {}

for curRun in runs:

    print("Cur run is %d"%(curRun) )

    cmd_Run = "./AnaGEMC.exe %d"%(curRun)

    proc_ana[curRun] = subprocess.Popen([cmd_Run], shell = True)
    run_counter = run_counter + 1

    if run_counter % 18 == 0:

        time.sleep(2)

        stillRunning = True

        while stillRunning:

            nProc = check_NumberofProcesses(proc_ana)
            
            if nProc > 8:
                print("* Still %d Analysis are running for this batch" %(nProc))
                print("* Sleeping...")
                time.sleep(10)
            else:
                print( "Going to start next batch of Analysis" )
                stillRunning = False

