import os
import sys
import subprocess


if __name__ == "__main__":

    if len(sys.argv) != 3 :
        print( "Wrong syntax, exiting" )
        print( "The command should look like python linkPars.py NewRun RunItShouldPointTo" )
        exit(1)

    processes = set()

    newRun = int(sys.argv[1])
    PointToRun = int(sys.argv[2])
    
    cmd_mupEloss = "ln -s mup_ElossFunc_Pars_Run_%d.dat mup_ElossFunc_Pars_Run_%d.dat" %(PointToRun, newRun)
    cmd_mumEloss = "ln -s mum_ElossFunc_Pars_Run_%d.dat mum_ElossFunc_Pars_Run_%d.dat" %(PointToRun, newRun)
    cmd_mupThCorr = "ln -s mup_DeltaTheta_Corr_Run_%d.dat mup_DeltaTheta_Corr_Run_%d.dat" %(PointToRun, newRun)
    cmd_mumThCorr = "ln -s mum_DeltaTheta_Corr_Run_%d.dat mum_DeltaTheta_Corr_Run_%d.dat" %(PointToRun, newRun)
    cmd_mupPhiCorr = "ln -s mup_DeltaPhi_Corr_Run_%d.dat mup_DeltaPhi_Corr_Run_%d.dat" %(PointToRun, newRun)
    cmd_mumPhiCorr = "ln -s mum_DeltaPhi_Corr_Run_%d.dat mum_DeltaPhi_Corr_Run_%d.dat" %(PointToRun, newRun)
    
    
    subprocess.Popen([cmd_mupEloss], shell = True)
    subprocess.Popen([cmd_mumEloss], shell = True)
    subprocess.Popen([cmd_mupThCorr], shell = True)
    subprocess.Popen([cmd_mumThCorr], shell = True)
    subprocess.Popen([cmd_mupPhiCorr], shell = True)
    subprocess.Popen([cmd_mumPhiCorr], shell = True)

    

    
