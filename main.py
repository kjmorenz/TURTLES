import sim as sim
import math

#control parameters - 1 means on, 0 means off
write = 1
analyze = 1
makefig = 1
normalize = 1

diffuse = 1
antibunch = 1
pulsed = 1

compare = 0

#General saving folder parameters
filepath = "C:/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/Simulation-cleanedMar92018/"
filedir = "2019-05-30-redotento9longtaufewerems/"
filedir = "2019-06-06-longonetry/"
filedir = "2019-06-13-newexseries/"
#filedir = "2019-04-15-diffmanyems/"
#filename = "Longerbutshortg2Lowerpwr" #fewerems
filedir = "2019-07-09/"
filename = "lowapyesdcnotperfsensmedlpshorttaulongerALLDARKS"
filedir = "New pulsed4-long/"
filename = "1p4ustau-short-lowerac"


#file length - goes to the max of these two
numlines = 10**8
maxlines = 5*10**8
endtime = 10*10**17#10**14#17# ms 1*10**12 #ps

#Sample parameters
temp = 298 #K
k_emission = 1.4*10**6# 10**4#1# emission lifetime in ms #10**9 #emission lifetime in ps
#k_emission = 9#,8,7,6,5]
#k_emission = 10000
emwavelength = 815#1000#
bfrate = 1 #fraction of successful emitting bright fission events - make sure to consider laser wavelength and emission wavelength
bfrates = [0,1]#0.3,0.5,0.7,,0.01,0.1,1


r = 10 #nm - hydrodynamic radius of particles

eta = 8.9 * 10**(-13) # kg/nm s - dynamic viscosity of solvent (water)
n = 1.3 # index of refraction - sample in water

concentration = 5*10**(-4) 
concentration = 0.00002#2*10**(-5)
concentration = 8 #1.81*10**8 emitters per focal vol according to other calc for Yb perov
concentration = 8*10**(-7)/18 #1 emitter at 405 nm focal vol
#concentrations = [concentration * 10, concentration, concentration / 10, concentration / 100, concentration / 1000]
concentration = 2*10**(-8) # 4.55 e-8 for 1 em avg in 405 focal vol

absXsec = 7.180616773321853*10**(-9)# per emitter numabs'd = phperpulse*absXsec*numEms - this is reasonable based on absXsec for CdSe is 550000*r^3/cm (from Bawendi paper Ruvim sent me)

#Laser parameters
reprate = 0.05#1 #MHz 

wavelength = 532 #nm405 #
beamdiam = 5000000 #5 mm in nm
laserpwr = 2.6133154984659886*10**(-8)#5*10**(-5)#2.6*10**(-8)#0.00015 # mW (0.52 mWinto back of objective 23 #mW)
laserpwr = 3.7740302137966845*10**(-9) # according to calculated num ems and ac 10**4
laserpwr = 3.7740302137966845*10**(-2) # according to calculated num ems and ac 10**4
laserpwr = 3.7740302137966845*10**(-8)
laserpwr = 3.7740302137966845*10**(-3)
#ac = 2*10**4 #overrides laser power
acs = [10,100,1000,10000,100000]
acs = [5000,15000,25000,30000,40000,50000]#10000,20000,
pulselength = 80 #ps - not used

accConc = 0#concentration/100
racc = 0#r
crashtransferprob = 1
tripem = 1
dftime = float("inf")# 10*k_emission

#Objective parameters
foclen = 310000 #310 microns in nm (working distance + coverslip thickness)
NA = 1.4

#Detector parameters
darkcounts = 100 #Hz per detector, assume 2 detectors
sensitivity = 0.01#001#001
deadtime = 70000# #70 ns in ps
afterpulse = 0.0001#001 #percent of time a photon is emitted a deadtime after one is detected

#Correlation parameters
gnname = ""
order = 2 #g2
mode = "t2"
gnpwr = 23#19#35#5#ms 
numbins = 256#2**10 #
pulsebins = (8)-1#should always be an odd number
channels = 2

#Miscellaneous simulation parameters
picyzoom = [-1,-1]
timestep = 1000#**12 #average number of photons per "round" of calculations - maybe broken
'''
write = 0
concentration = concentration * 10
write = 1
for concf in range(2):
    for bfrate in bfrates:
        f = filename + "bf" + str(bfrate)
        if concf == 0:
            f = f+ "-10ems"
        elif concf == 1:
            f = f + "-100ems"
        else:
            f = f + "1em"
        sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                        pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                        k_emission, emwavelength, bfrate, r,eta, n, 
                        accConc, racc, crashtransferprob, tripem, dftime,
                        reprate, wavelength, laserpwr, pulselength, foclen,
                        NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                        mode, gnpwr, numbins, pulsebins, channels, 
                        dopic = 0, normalize = normalize, 
                        units = 'new', timestep=timestep, foutedit = gnname)
    concentration = concentration * 10
    write = 1


'''

for diffuse in [0,1]:
    for ac in acs:
        f = filename + "avgexc" + str(ac/((reprate*10**6))) + "perpulse"
        if diffuse == 1:
            f = f +"diff"
        print()
        print(f)
        print()
        sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                                k_emission, emwavelength, bfrate, r,eta, n, 
                                accConc, racc, crashtransferprob, tripem, dftime,
                                reprate, wavelength, laserpwr, pulselength, foclen,
                                NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                                mode, gnpwr, numbins, pulsebins, channels, ac = ac*sensitivity, 
                                dopic = 0, normalize = normalize, 
                                units = 'new', timestep=timestep)
'''
    sim.simulate(filepath, filedir, f + "old", write, analyze, makefig, diffuse, antibunch,
                            pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                            k_emission, emwavelength, bfrate, r,eta, n, 
                            accConc, racc, crashtransferprob, tripem, dftime,
                            reprate, wavelength, laserpwr, pulselength, foclen,
                            NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                            mode, gnpwr, numbins, pulsebins, channels, ac = ac*sensitivity, dopic = 0, normalize = normalize, 
                            units = 'old', timestep=timestep)

for bfrate in bfrates:
    f = filename + "pmqy" + str(bfrate)
    sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                            pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                            k_emission, emwavelength, bfrate, r,eta, n, 
                            accConc, racc, crashtransferprob, tripem, dftime,
                            reprate, wavelength, laserpwr, pulselength, foclen,
                            NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                            mode, gnpwr, numbins, pulsebins, channels, dopic = 0, normalize = normalize, 
                            units = 'new', timestep=timestep)


                
    

for bfrate in bfrates:
    
#   
    #for sensitivity in [0.001, 0.0001, 0.000001]:
    f = filename + "bf" + str(int(100*bfrate)) + "sens" + str(sensitivity)

    sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                                k_emission, emwavelength, bfrate, r,eta, n, 
                                accConc, racc, crashtransferprob, tripem, dftime,
                                reprate, wavelength, lp, pulselength, foclen,
                                NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                                mode, gnpwr, numbins, pulsebins, channels, dopic = 0, normalize = normalize, timestep=timestep)

concentration = concentration / 100
for bfrate in bfrates:
    if bfrate == 0:
        write = 1
#   f = filename
    for sensitivity in [0.001, 0.0001, 0.000001]:
        f = "0.01conc" + filename + "bf" + str(int(100*bfrate)) + "sens" + str(sensitivity)

        sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                                k_emission, emwavelength, bfrate, r,eta, n, 
                                accConc, racc, crashtransferprob, tripem, dftime,
                                reprate, wavelength, lp, pulselength, foclen,
                                NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                                mode, gnpwr, numbins, pulsebins, channels, dopic = 0, normalize = normalize, timestep=timestep)
concentration = concentration * 100
laserpwr = laserpwr / 100
for bfrate in bfrates:
    
#   f = filename
    for sensitivity in [0.001, 0.0001, 0.000001]:
        f = "0.01lp" + filename + "bf" + str(int(100*bfrate)) + "sens" + str(sensitivity)

        sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                                k_emission, emwavelength, bfrate, r,eta, n, 
                                accConc, racc, crashtransferprob, tripem, dftime,
                                reprate, wavelength, lp, pulselength, foclen,
                                NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                                mode, gnpwr, numbins, pulsebins, channels, dopic = 0, normalize = normalize, timestep=timestep)

---

    concentration = concentration/(10**6)
    for i in range(3):
        concentration = concentration*10**2
        f = filename + "bf" + str(int(bfrate*100)) + "conc" + str(concentration)
        
            

    

    accConc = concentration/1000
    racc = 10
    concentration = concentration/1000
    filenames = ["NoDF-", "DF-"]
    for i in range(2):
        concentration = concentration * 10 * (10**i)
        for bf in bfrate:
            for j in range(2):
                tripem = 10**(7+2*j)
                for k in range(3):
                    accConc = accConc * 10**(2*k)
                    for filename in filenames:
                        if filename == "DF-":
                            for i in range(3):
                                dftime = k_emission * 10 ** (2*i)
                                f = filename + str(int(dftime/100)) + "-AC-" + str(k) + "-Tem-" + str(j)
                                if pulsed == 1:
                                    f= f + "-1MHz"

                                sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                        pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                                        k_emission, emwavelength, bf, r,eta, n, 
                                        accConc, racc, crashtransferprob, tripem, dftime,
                                        reprate, wavelength, laserpwr, pulselength, foclen,
                                        NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                                        mode, gnpwr, numbins, pulsebins, channels)
                        else:
                            dftime = float("inf")
                            f = filename + "-AC-" + str(k) + "-Tem-" + str(j) 
                            if pulsed == 1:
                                f= f + "-1MHz"

                            sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                        pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                                        k_emission, emwavelength, bf, r,eta, n, 
                                        accConc, racc, crashtransferprob, tripem, dftime,
                                        reprate, wavelength, laserpwr, pulselength, foclen,
                                        NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                                        mode, gnpwr, numbins, pulsebins, channels) '''

    

    
    


                


