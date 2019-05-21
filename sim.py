import os
import subprocess
import math
import matplotlib.pyplot as plt
import photon_correlation as pc

import makefig2 as m

import write2 as w
import writens as wns
import analyze2 as a

def isbf(file):
    rval = 'r'
    for i in range(len(file)-1):
        if file[i] == 'b' and file[i+1] == 'f':
            rval = 'b'
            break

    return rval

def simulate(filepath, filedir, fullfilename, write = 1, analyze = 1, makefig = 1, diffuse = 1, antibunch = 1,
             pulsed = 0, numlines = 10**7, maxlines = 10**8, endtime = 10**12,
             temp = 298, concentration = 2*10**-8, absXsec = 3.6*10**-10,
             k_emission = 10000, emwavelength = 815, bfrate = 0, r = 10,
             eta = 8.9 * 10**(-13), n = 1.3, accConc = 0, racc = 0, 
             crashtransferprob = 1 , tripem = 1, dftime = float("inf"),
             reprate = 1, wavelength = 532, laserpwr = 0.5, pulselength = 80, foclen = 310000,
             NA = 1.4, darkcounts = 1, sensitivity = 0.1, deadtime = 70000, afterpulse = 0, order = 2,
             mode = "t2", gnpwr = 20, numbins = 4096, pulsebins = 99, channels = 3,
             picyzoom = 100, timestep = 200, ac = -1, compare = 0, dopic = 1, normalize = 0, units = 'ps'):
    
    if units == 'ns':
        writefn = wns.write
    
    else:
        writefn = w.write
    
    suffix = ".txt"
    if not os.path.isdir(filepath):
        os.mkdir(filepath)

    if not os.path.isdir(filepath + "RawData/"):
         os.mkdir(filepath + "RawData/")
         os.mkdir(filepath + "Figures/")

    if not os.path.isdir(filepath + "RawData/"+ filedir):
        os.mkdir(filepath + "RawData/"+ filedir)
        os.mkdir(filepath + "Figures/"+ filedir)
        os.mkdir(filepath + "RawData/"+ filedir+fullfilename)
        os.mkdir(filepath + "Figures/"+ filedir+fullfilename)

    elif not os.path.isdir(filepath + "RawData/"+ filedir+fullfilename):
        os.mkdir(filepath + "RawData/"+ filedir+fullfilename)
        os.mkdir(filepath + "Figures/"+ filedir+fullfilename)

    os.chdir(filepath + "RawData/"+ filedir+fullfilename)

    print(os.getcwd())
    if write == 1:
        writefn(filepath, filedir, fullfilename, antibunch, diffuse, pulsed, numlines, maxlines, endtime,
             temp, concentration, absXsec,k_emission, emwavelength, bfrate, r,
             eta, n, accConc, racc, crashtransferprob, tripem, dftime, 
             reprate,wavelength, laserpwr, pulselength, foclen,
             NA, darkcounts, sensitivity, deadtime, afterpulse, timestep, channels, ac, compare)

    if analyze == 1:
        a.analyze(filepath, filedir, fullfilename, numlines, order, mode, gnpwr, 
                    numbins, pulsebins, channels, makefig, m.makeafig, pulsed, picyzoom, reprate, deadtime, dopic, normalize)


    elif makefig == 1: #only really used if it crashed part way through analysis
        
        file = fullfilename
        fileoutname = file
        c = isbf(fullfilename)
        filename = "fig"
  
        m.makeafig("g2", filename, [-1,-1], [-1,-1], 0, pulsed, normalize = normalize, filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
        filename = "750nszoom"
        m.makeafig("g2", filename, [-75,75], [-1,-1], 0, pulsed, normalize = normalize,filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
        filename = "250nszoom"
        m.makeafig("g2", filename, [-250,250], [-1,-1], 0, pulsed, normalize = normalize,filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
        filename = "log"
        m.makeafig("g2", filename,[-1,-1],[-1,-1], 1, pulsed, normalize = normalize, filepath = filepath, filedir = filedir + file+"/", fileoutdir = fileoutname+".g2.run/", color = c)
'''
        print(filepath+"/"+file)
        fileout = file + "-PIC"
        filename = "PIC"
        m.makeafig(fileout, filename,[-1,-1], [-1,-1], 1, pulsed,filepath+"/"+file, color = c)
        filename = filename + "xzoom"
        m.makeafig(fileout, filename, [0,10**6], picyzoom, 1, pulsed,filepath+"/"+file, color = c)
        filename = "PICtight"
        m.makeafig(fileout,filename,[0.0001,10], picyzoom, 1, pulsed,filepath+"/"+file, color = c)
        filename = "PICmed"
        m.makeafig(fileout,filename,[0,10**10], picyzoom, 1, pulsed,filepath+"/"+file, color = c)'''




