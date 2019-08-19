import os
import subprocess
import matplotlib.pyplot as plt
import photon_correlation as pc
import numpy
import newmakefigpulse as npf
import gc
def isbf(file):
    rval = 'r'
    for i in range(len(file)-1):
        if file[i] == 'b' and file[i+1] == 'f':
            rval = 'b'
            break

    return rval

def analyze(filepath, filedir, fullfilename, numlines, order, mode, gnpwr, 
            numbins, pulsebins, channels, makefig, makeafig, pulsed, picyzoom =[-1,-1], 
            reprate = 1, deadtime = 70000, dopic =1, normalize = 0, binwidth = 10**12, foutedit = ""):
    
    dfilepath = filepath + "RawData/" + filedir
    ffilepath = filepath + "Figures/" + filedir
    os.chdir(dfilepath+fullfilename)
    suffix = ".txt"
    t2time = "-"+str(2**gnpwr)+","+str(numbins)+","+str(2**gnpwr) #min, numbins, max in ps
    time = t2time
    fileoutname = fullfilename + foutedit
    if pulsed == 1:
        pulse = "-"+str((pulsebins)/2)+","+str(pulsebins)+","+str(pulsebins/2)
        taurep = (10**12)/(reprate*10**6) #ps
        gn = str(2**(int(numpy.log2(taurep/2))+1))
        time = "-"+gn+","+str(numbins)+","+gn
       
    
        
    file = fullfilename
    cmd = "photon_gn,"
    if pulsed == 1:
        cmd1 = ("photon_synced_t2,--sync-channel," + str(channels) +",--file-in," + dfilepath +file+ "/" + file + suffix)
        cmd1 = cmd1.split(",")
        photons = subprocess.Popen(cmd1, stdout = subprocess.PIPE)
        mode = "t3"
    else:
        cmd = cmd + ("--file-in," + dfilepath +file+ "/" + file + suffix + ",")
        mode = "t2"

    
    cmd = (cmd + "--order," + str(order) + ",--mode," + mode +",--channels," + str(channels) + ",--file-out," + fileoutname)
    bashCommand = cmd.split(",")
    bashCommand.append("--time")
    bashCommand.append(time)
    if pulsed == 1:
        bashCommand.append("--pulse")
        bashCommand.append(pulse)
    
    bashCommand.append("--print-every")
    bashCommand.append(str(int(numlines/20)))#every 5%


    print(bashCommand)
    print(os.getcwd())
    if pulsed != 1:           
        process = subprocess.Popen(bashCommand)
        output, error = process.communicate()
    else:
        process = subprocess.Popen(bashCommand, stdin = photons.stdout)
        output, error = process.communicate()
    process = subprocess.Popen("pwd")
    output, error = process.communicate()

    cmd = ("photon_intensity,")
    if pulsed == 1:
        cmd1 = ("photon_synced_t2,--sync-channel," + str(channels) +",--file-in," + dfilepath +file+ "/" + file + suffix)
        cmd1 = cmd1.split(",")
        photons = subprocess.Popen(cmd1, stdout = subprocess.PIPE)
        mode = "t3"
    else:
        cmd = cmd + ("--file-in," + dfilepath +file+ "/" + file + suffix + ",")
        mode = "t2"
    cmd = cmd + ("--mode," + mode + ",--channels," + str(channels) + ",--bin-width," + str(binwidth)+ ",--file-out," + fileoutname + str(binwidth)+ "bins.txt")
    bashCommand = cmd.split(",")
    print("Binning started")
    if pulsed == 1:
        proc = subprocess.Popen(bashCommand, stdin = photons.stdout)
    proc = subprocess.Popen(bashCommand)
    print("Process done")
    print(fileoutname + str(binwidth)+ "bins.txt")
    print("Working directory = ")
    process = subprocess.Popen("pwd")
    output, error = process.communicate()
    

    if makefig == 1:
        c = isbf(file)

        if pulsed == 0:
            filename = fileoutname + "fig"
            print(dfilepath+file+"/"+fileoutname+".g2.run")
            makeafig("g2", filename, [-1,-1], [-1,-1], 0, pulsed, normalize = normalize, filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
        else:
            npf.makepulsedfig(dfilepath+file+"/"+fileoutname+".g2.run/", "g2", ffilepath+file+"/", fileoutname, deadtime = deadtime, reprate = reprate, timespace = 1000)
            npf.makepulsedfig(dfilepath+file+"/"+fileoutname+".g2.run/", "g2", ffilepath+file+"/", fileoutname, deadtime = deadtime, reprate = reprate, timespace = 1000, xzoom = 2000)
            #npf.makepulsedfig(dfilepath, "g2", sfilepath, fileoutname+"250nszoom", reprate = reprate, timespace = 1000, xzoom = 100) 
        
        


    os.chdir(dfilepath+file)
    

    if dopic == 1:
        print(os.getcwd())
        fileout = file + "-PIC"
        cmd = "photon_intensity_correlate,"
        if pulsed == 1:
            cmd1 = ("photon_synced_t2,--sync-channel," + str(channels) +",--file-in," + dfilepath +file+ "/" + file + suffix)
            cmd1 = cmd1.split(",")
            photons = subprocess.Popen(cmd1, stdout = subprocess.PIPE)
        else:
            cmd = cmd + ("--file-in," + dfilepath +file+ "/" + file + suffix + ",")
        cmd = cmd + ("--order," + str(order) + ",--mode," + mode +
                ",--channels," + str(channels) + ",--file-out," + fileout)
        bashCommand = cmd.split(",")
        bashCommand.append("--time-scale")
        bashCommand.append("log-zero")
        print(bashCommand)

        if pulsed != 1:           
            process = subprocess.Popen(bashCommand)
            output, error = process.communicate()
        else:
            process = subprocess.Popen(bashCommand, stdin = photons.stdout)
            output, error = process.communicate()

        
        if makefig == 1:
            print(os.getcwd())
            filename = "PIC"
            makeafig(fileout, filename,[-1,-1], [-1,-1], 1, 0,filepath = filepath, filedir = filedir + file+"/", color = c)
        

