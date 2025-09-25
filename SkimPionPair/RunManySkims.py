import subprocess
import shlex
import os
import time
import glob
import sys

def check_NumberofProcesses(proc_dict):
    nRunningProcess = 0
    for file_ind in proc_dict:
        if proc_dict[file_ind].poll() is None:
            nRunningProcess = nRunningProcess + 1
    return nRunningProcess;

def WaitWhileRunning(proc, procDescription):
    StillRunning = True
    while StillRunning:
        if proc.poll() is None:
            print("* Waiting 5 seconds, while %s is running" % (procDescription))
            time.sleep(5);
        else:
            return 1



if __name__ == "__main__":

    if len(sys.argv) != 2 :
        print( "Wrong syntax, exiting" )
        print( "The command should look like RunManySkims.py $Run" )
        exit(1)


    processes = set()
    run = int(sys.argv[1])

    nPerFile = 4500
    

    files = glob.glob("Data/%d/rec_clas*.hipo" %(run));

    
    proc_Skim = {}

    file_counter = 0

    for curFile in files:
        splited_fname = curFile.split(".")

        file_ind = splited_fname[2]
        print("file_ind = %s"%(file_ind))

        skim_cmd = "./SkipPionPair.exe %d %s %d"%(run, file_ind, nPerFile)
        #skim_cmd = "./SkimSinglePin.exe %d %s %d"%(run, file_ind, nPerFile)
        
        proc_Skim[file_ind] = subprocess.Popen([skim_cmd], shell = True)

        file_counter = file_counter + 1        

        if file_counter % 18 == 0:

            time.sleep(2)

            stillRunning = True

            while stillRunning:

                nProc = check_NumberofProcesses(proc_Skim)

                if nProc > 8:
                    print("* Still %d Skims are running for this batch" %(nProc))
                    print("* Sleeping...")
                    time.sleep(10)
                else:
                    print( "Going to start next batch of Skims" )
                    stillRunning = False

    time.sleep(5)

    stillRunning = True

    while stillRunning:
        nRunningProcess = 0

        print("* Checking if Skims of all files is finished...")
        for file_ind in proc_Skim:

            print( "Process Status for %d process is %s" %( proc_Skim[file_ind].pid, proc_Skim[file_ind].poll() ) )

            if proc_Skim[file_ind].poll() is None:
                nRunningProcess = nRunningProcess + 1
            print("The Skim for the file %s is still running" %(file_ind))

        if nRunningProcess == 0:
            stillRunning = False;
            break
        else:
            print("There are still %d Skims are running." %(nRunningProcess) )
            wait_sec = 30;
            print("Waiting for %d seconds" %(wait_sec))
            time.sleep(wait_sec)



    print(" \n\n\n\n * Adding all root files together")

    cmd_hadd = "hadd -f Skim_PionPair_Run_%d_All.root Skim_PionPair_%d_*.root"%(run, run)
    #    Skim_PionPair_5117_00000-00004.root
    proc_Hadd = subprocess.Popen([cmd_hadd], shell = True)
    
    WaitWhileRunning(proc_Hadd, "Hadd")

    cmd_rmRoot = "rm -f Skim_PionPair_%d_*.root"%( run )
    proc_rm = subprocess.Popen([cmd_rmRoot], shell = True)

    WaitWhileRunning(proc_rm, "Cleaning root files")

    
    print("All done.")
