import sim as sim
import math

#control parameters - 1 means on, 0 means off
write = 1
analyze = 1
makefig = 1
normalize = 1

diffuse = 0
antibunch = 1
pulsed = 0

compare = 0

#General saving folder parameters
filepath = "C:/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/Simulation-cleanedMar92018/"
filedir = "2019-05-07-perov-nofit/"
#filedir = "2019-04-15-diffmanyems/"
#filename = "Longerbutshortg2Lowerpwr" #fewerems
filename = "nslowsens2"

#file length - goes to the max of these two
numlines = 10**4
maxlines = 10**5
endtime = 1*10**9 #ns

#Sample parameters
temp = 298 #K
k_emission = 10**6 #emission lifetime in ns
#k_emission = 9#,8,7,6,5]
#k_emission = 10000
emwavelength = 1000#815
#bfrate = 0.1 #fraction of successful emitting bright fission events - make sure to consider laser wavelength and emission wavelength
bfrates = [0.01,1]#0.1,0.3,0.5,0.7,

r = 10 #nm - hydrodynamic radius of particles

eta = 8.9 * 10**(-13) # kg/nm s - dynamic viscosity of solvent (water)
n = 1.3 # index of refraction - sample in water

concentration = 5*10**(-4) 
concentration = 2*10**(-2)

absXsec = 7.180616773321853*10**(-9)# per emitter numabs'd = phperpulse*absXsec*numEms - this is reasonable based on absXsec for CdSe is 550000*r^3/cm (from Bawendi paper Ruvim sent me)

#Laser parameters
reprate = 1 #MHz

wavelength = 405#532 #nm
beamdiam = 5000000 #5 mm in nm
laserpwr = 2.6*10**(-8)#0.00015 # mW (0.52 mWinto back of objective 23 #mW)
#laserpwr = [0.05,0.1,0.2,0.3,0.5,0.75,1,10]#[0.001,0.01,0.025,

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
darkcounts = 200#200 #s^-1
sensitivity = 0.0001
sens = [0.0001, 0.000001, 0.001] 
deadtime = 80 #ns
afterpulse = 0.0001 #percent of time a photon is emitted a deadtime after one is detected

#Correlation parameters
order = 2 #g2
mode = "t2"
gnpwr = 25#20
numbins = 2**10 #256
pulsebins = (8)-1#should always be an odd number
channels = 2

#Miscellaneous simulation parameters
picyzoom = [-1,-1]
timestep = 1000#**12 #average number of photons per "round" of calculations - maybe broken



#dts = [70000,0]
lp = laserpwr
N_A = 6.023*10**23
for bfrate in bfrates:
    for sensitivity in sens:
        f = filename + "bf" + str(int(100*bfrate)) + "sens" + str(sensitivity)
        sim.simulate(filepath, filedir, f, write, analyze, makefig, diffuse, antibunch,
                                pulsed, numlines, maxlines, endtime,temp, concentration, absXsec,
                                k_emission, emwavelength, bfrate, r,eta, n, 
                                accConc, racc, crashtransferprob, tripem, dftime,
                                reprate, wavelength, lp, pulselength, foclen,
                                NA, darkcounts, sensitivity, deadtime, afterpulse, order, 
                                mode, gnpwr, numbins, pulsebins, channels, dopic = 0, normalize = normalize, ac = 10**4, timestep=timestep, units = 'ns')
'''

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

    

    
    


                


