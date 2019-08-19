
import diffuse as d
import calcnextphotondda as calc
import calcnextphoton3 as simpcalc
import time as rt
import math
import numpy

def insert(array, val):
    if array == []:
        return [val]
    if array[-1] < val:
        array.append(val)
        return array
    if array[0] > val:
        array.insert(0, val)
        return array
    index = int(len(array)/2)
    while array[index] < val:
        index = index + 1
    while array[index - 1] > val:
        index = index - 1
    array.insert(index, val)
    return array

def addafterpulse(photons, darks, deadtime, afterpulse):
    for i in range(len(photons)):
        apphoton = photons[i]+ deadtime
        index = 0
        if numpy.random.rand() < afterpulse:#includes probability that went to other detector
            if apphoton < photons[-1]:  #ignore afterpulsing after end of round
                while photons[i+index] < apphoton:
                    index = index + 1
                index = index - 1
                if apphoton - photons[index] > deadtime: #so long as the second detector wasn't already turned off by another photon arrival 
                    photons.insert(index-1, apphoton)
            else:#add it to dark counts instead so it'll get added in next round
                darks = insert(darks, apphoton)
    return photons, darks
            
            

def write(filepath, filedir, fullfilename, antibunch, diffuse, pulsed, numlines, maxlines, endtime,
            temp, concentration, absXsec,k_emission, emwavelength, bfrate, r,
            eta, n, accConc, racc, crashtransferprob, tripem, dftime, 
            reprate,wavelength, laserpwr, pulselength, foclen,
            NA, darkcounts, sensitivity, deadtime, afterpulse, timestep, channels, ac, compare):
    allparams = [filepath, filedir, fullfilename, antibunch, diffuse, pulsed, numlines, maxlines, endtime,
            temp, concentration, absXsec,k_emission, emwavelength, bfrate, r,
            eta, n, accConc, racc, crashtransferprob, tripem, dftime, 
            reprate,wavelength, laserpwr, pulselength, foclen,
            NA, darkcounts, sensitivity, deadtime, afterpulse, timestep, channels]
    
    h = 6.62607004*10**(-34) # m^2 kg / s (Planck's const)
    c = 299792458 # m/s (Speed of light)
    N_A = 6.023*(10**23) # Avogadro
    K_B = 1.38064852*10**(-5) # Boltzmann in kg nm^2 /K s^2
    suffix = ".txt"
    #beamwaist = 2*wavelength*foclen/(math.pi*beamdiam) #nm
    #this gives too small a value (~20 nm) - gets us very fast diffusion times, too fast to be real
    ##beamwaist = 4lambda f / 2 pi D where f is the focal length and D is the diameter of light entering
    ## or lambda/2NA

    beamwaist = wavelength/(2*NA) #~200 nm

    focalVol = (math.pi**(1.5))*(2.33*n/NA)*beamwaist**3 #nm^3

    AvgEms = (concentration * N_A * focalVol / (10**24)) #number not mol
    #AvgEms = 3*10**4
    taurep = (10**12)/(reprate*10**6) #ps

    
    energyperphoton = h*c/(wavelength/(10**9)) #J

    if ac > -1:
        if pulsed == 1:#untested
            exsperpulse = ac/(sensitivity*reprate*10**6)
            experemperpulse = exsperpulse/AvgEms
            phperpulse = experemperpulse/absXsec
            energyperpulse = phperpulse*energyperphoton
            laserpwr = energyperpulse*reprate*1000000000
        else:
            AvgEmEvs = ac/((10**12)*sensitivity)
            k_excitation = AvgEms/AvgEmEvs - k_emission
            phpers = 10**12/(k_excitation*absXsec)
            laserpwr = phpers*(energyperphoton*1000)
            


    energyperpulse = laserpwr/(reprate * 1000000000) #J
    phperpulse = energyperpulse/energyperphoton
    phpers = laserpwr/(energyperphoton*1000)
    #make rounds average timestep emissions per round:
    print(timestep)
    if pulsed == 1:
        k_excitation = 10**6/(phperpulse*absXsec*reprate)#*AvgEms) #ps/excitation
        experemperpulse = absXsec*phperpulse
        if experemperpulse > 1:
            experemperpulse = 1
        exsperpulse = absXsec*phperpulse*AvgEms
        AvgExEvs = exsperpulse*reprate*10**6
        experps = experemperpulse*taurep*AvgEms # comment *avgems out to make per em
        timestep = timestep/experps #ps per round
        if timestep < taurep:
            timestep = taurep
        probex = 1- (1-absXsec)**phperpulse #1-probability of not being excited each pulse
    #CW version:
    
    else:
        k_excitation = 10**12/(phpers*absXsec)#*AvgEms) #ps/excitation
        AvgEmEvs = AvgEms/(k_emission + k_excitation) # average sample emissions per second (in pHz)
        timestep = timestep/AvgEmEvs # in ps
        probex = 1
        exsperpulse = -1
        AvgExEvs = 10**12/k_excitation

    print("kex = " + str(k_excitation))
    print(timestep)
    if timestep < k_excitation + k_emission:
        timestep = 2*(k_excitation + k_emission)
    endround = timestep
    #timestep = timestep * 100
    print(str(timestep) + " ps per round")


    ##focal vol = pi^1.5 kw^3, k is a geometric factor usually taken as optical resolution in the z direction
    ##divided by that in the xy direction, w is the e^-2 beam radius in the xy plane
    ##k >= 2.33 n/NA, r_z = 2nlambda/NA^2, rxy = 0.61 lambda/NA
    ##http://www.fcsxpert.com/classroom/theory/what-is-confocal-volume.html

    diffouttime = (3*((beamwaist)**2)*(math.pi**2)*eta*r/(8*K_B*temp))*10**12
    #= 3 (b pi)^2 eta r / 8 k_B T
    ##sigma = numpy.sqrt(2*beamwaist/(6*math.pi*eta*r))

    #gets a similar order of magnitude if x = focalVol^(2/3)

    ##Phys Chem by Atkins pg 771 - prob of being a distance x from origin after time t is
    ##P = sqrt(2tau/pi t)*e^(-x^2 tau/(2t lambda^2) where tau is the time it takes to go a distance lambda
    ##Diffusion coeff D = lambda^2/2tau = uRT/zF = uK_BT where u is the mobility or the ratio of the 
    ##terminal drift velocity to an applied force
    ##lambda = zuF, also u = 1/ 6 pi eta r Stokes Einstein
    ##
    ##Better: Wikipedia, (also in Atkins) Einstein relation: 
    ## D = k_B T/(6 pi eta r), eta is the viscosity, r is the radius of the spherical particle
    ##
    ##or pg 770, mean distanced travelled by particles of radius r in solvent of viscosity eta is
    ##                       x = sqrt(2kTt/3 pi^2 eta r)

    ##https://www2.gwu.edu/~phy21bio/Reading/randomwalkBerg.pdf gives sigma_x is sqrt(2Dt), P(x) = (sqrt(4 pi Dt)^-1)e^((-x^2)/4Dt))

    # combination of Atkins and Wikipedia gives:
    # mean free path lambda = (root(2)*conc*collision cross sectional area)
    # --> (using equations about D above) Mean free time = 3 eta / 2 conc ^2 pi r^3 Kb T where conc is only acc dependent
    # and r is the sum of the two diameters 

    if not accConc == 0:
        crashtime = 3*eta / (2*accConc**2*numpy.pi*((racc + r)**3)*K_B*temp)  # tau dep on each things conc, mean free path to encounter OPPOSITE species
    else:
        crashtime = float("inf")

    deltat = 0
    line = 0
    sigcts = 0
    nextOut = 0
    nextIn = 0
    diffsIn = 0
    diffsOut = 0
    ndiffsOut = 0
    lastwritten = -deadtime - 1 #both detectors on
    lastlastwritten = lastwritten
    dclen = 20
    channel = numpy.random.randint(channels)

    

    

    numEms = round(AvgEms) #start with the average number of emitters in the focal volume
    if numEms <= 0: #always start with at least one emitter in the focal volume
        numEms = 1
    print("NumEms = " + str(numEms))
    print("AvgEms = " + str(AvgEms))
    


    lastphoton =[]
    bfs = []
    lastcrash = []
    for i in range(numEms):
        lastphoton.append(0)
        lastcrash.append(0)

    lastdiff = 0 # for calculating actual avg ems
    tnumEms = 0

    if diffuse == 1:
        nextOut = d.diffuse(diffouttime, numEms)#the number which have a chance to diffuse out are the number in
        nextIn = d.diffuse(diffouttime, AvgEms)#the number which have a chance to diffuse in are the average number
        nextdiff = min(nextIn, nextOut)
        
    elif diffuse == 0:
        nextdiff = float('inf')
        nextOut = float('inf')
        nextIn = float('inf')

    lastsyncpulse = 0
    if pulsed != 1:
        lastsyncpulse = float('inf')

    avtime = (1/darkcounts)*(10**12) #average time between darkcount photons in ps given darkcounts in s^-1
    darks = calc.darkcounts(avtime, dclen, 0, deadtime)

    dcpointer = 0
    countbright = 0

    pfile = open(filepath +"Figures/"+ filedir + fullfilename +"/" + "params-"+ fullfilename + suffix, 'w')
    pfile.write(fullfilename +" - date and time: " + str(rt.time()) + " \n")
    
    
    pfile.write("Definitive parameters:")
    lvars = locals()
    for i in lvars:
        pfile.write(i + " : " + str(lvars[i]) + "\n")
    pfile.write("\n")
    pfile.close()

    file = open(filepath +"RawData/"+ filedir+ fullfilename +"/" + fullfilename+ suffix, 'w')
    difffile = open(filepath +"RawData/"+ filedir+ fullfilename +"/" + "diffsof" + fullfilename+ suffix, 'w')
    starttime = rt.time()
    wrotediff = [0,0]
    while (line < numlines or lastwritten < endtime):# and line < maxlines:#comment out maxlines term to force to get to endtime
        if line > maxlines:
            break
        photons = []
        
                    
            
        if dcpointer == dclen:
            darks = calc.darkcounts(avtime, dclen, darks[dcpointer-1], deadtime)
            dcpointer = 0

        if numEms < 1:
            print("NumEms = " + str(numEms))
            print("AvgEms = " + str(AvgEms))
            print("Rounded = " + str(round(AvgEms)))
            please = crash
        
        
        if len(lastphoton) != numEms or (len(lastcrash) != numEms and crashtime != float("inf")):
            print("numEms:" + str(numEms))
            print("lastphoton:" +str(lastphoton))
            print("lastcrash:"+str(lastcrash))
            print("write line 250")
        
        if wrotediff != [nextIn, nextOut] and nextIn < nextOut:
            difffile.write(str(nextIn) + ","+ str(numEms + 1)+"\n")
            wrotediff = [nextIn, nextOut]
        elif wrotediff != [nextIn, nextOut] and nextOut < nextIn and numEms >=0:
            difffile.write(str(nextOut) + ","+ str(numEms - 1)+"\n")
            if numEms -1 == 0:
                difffile.write(str(nextIn) + ",1\n")
            wrotediff = [nextIn, nextOut]
        #print("New Round")

        if compare == 0: #dftime == float("inf") and crashtime == float("inf") and compare == 0:
            dataphotons, lastphoton, bfs, numEms, nextIn, nextOut, diffsIn, diffsOut, ndiffsOut, countbright = simpcalc.nextphoton(lastphoton, sensitivity,
                                        k_emission, k_excitation, diffouttime,
                                        numEms, nextOut, nextIn, endround, timestep,
                                        antibunch, pulsed, taurep, bfrate,
                                        absXsec, phperpulse, AvgEms, diffsIn, 
                                        diffsOut, ndiffsOut, countbright, bfs, probex)
            '''if lastwritten > 70000000000:
                print(dataphotons)
                print(lastphoton[0])
                crash = plz'''
        else:
            dataphotons, lastphoton,  numEms, nextIn, nextOut, diffsIn, diffsOut, ndiffsOut, countbright,lastcrash = calc.nextphoton(lastphoton, sensitivity,
                                            k_emission, k_excitation, diffouttime,
                                            numEms, nextOut, nextIn, endround, 
                                            antibunch, pulsed, taurep, bfrate,
                                            absXsec, phperpulse, AvgEms, diffsIn, 
                                            diffsOut, ndiffsOut, countbright, dftime, tripem, crashtime, crashtransferprob, lastcrash, probex)
             
        endround = endround + timestep
        if len(lastphoton) == 1 and endround < lastphoton[0]:
            endround = lastphoton[0]
    ##            print("new dataphotons")
    ##            print(dataphotons)


        if not dataphotons == []: 
            sigcts += len(dataphotons)
            
            datapointer = 0
                        
            while datapointer < len(dataphotons):
                if dcpointer == dclen:
                    darks = calc.darkcounts(avtime, dclen, darks[dcpointer-1], deadtime)
                    dcpointer = 0
                        
                if darks[dcpointer] < dataphotons[datapointer]:
    ##                        print(0)
                    if darks[dcpointer] - lastwritten > deadtime or (darks[dcpointer] - lastlastwritten > deadtime and numpy.random.rand()>0.5): 
                        photons.append(darks[dcpointer])
                        lastlastwritten = lastwritten
                        lastwritten = darks[dcpointer]
    ##                            print("wrote a dc - 3")

                    dcpointer = dcpointer + 1
                    

                else:
    ##                            if datapointer >= len(dataphotons) or photons == []:
    ##                                print(darks)
    ##                                print(dcpointer)
    ##                                print(dataphotons)
    ##                                print(datapointer)
    ##                                print(photons)
    ##                                print(dataphotons[datapointer] + photons[0])
                    if dataphotons[datapointer] - lastwritten > deadtime or (dataphotons[datapointer] - lastlastwritten > deadtime and numpy.random.rand()>0.5):
                        photons.append(dataphotons[datapointer])
                        lastlastwritten = lastwritten
                        lastwritten = dataphotons[datapointer]
    ##                            print("wrote data - 1")
                        
                    datapointer = datapointer + 1
    ##                        print(datapointer)




        #the first time through, we calculate how often to print progress
        if deltat == 0:
            deltat = (rt.time() - starttime)*1000
            print("Delta t for one round: " + str(deltat))
            print("Timestep = " + str(timestep))
            #print(photons)

        if not photons == []:
            if not afterpulse == 0:
                photons, darks = addafterpulse(photons, darks, deadtime, afterpulse)
            #writing into file
            if photons[0] < 0:
                print(line + " is neg - didn't write")
            elif photons[0] == float('inf'):
                print(line + " is inf - didn't write")
            else:
                while photons[0] > lastsyncpulse:
                    file.write(str(channels) + "," + str(lastsyncpulse)+"\n")
                    lastsyncpulse = lastsyncpulse + int(taurep)
                    line = line + 1

                if photons[0] - lastwritten <= deadtime: #leq since afterpulse photons have to go to other det too
                    channel = (channel+1)%channels #went to the other detector
                else:
                    channel = numpy.random.randint(channels)
                file.write(str(channel) + "," + str(int(photons[0]))+"\n")
    ##            print("wrote a line")
                line = line + 1
                if line%numlines==0:
                    print("time = " + str(lastwritten/10**12) + " s")
                    if lastwritten>endtime:
                        break #should never happen
                if line%(numlines/10) == 0:
                    print(fullfilename + ": " + str(line*100/numlines) + "% - " + str(photons[0]/10**9) + " ms")
                elif deltat > 1 and line%(numlines/100)==0:#if it's slow print more often
                    print(str(line*100/numlines) + "% - " + str(photons[0]/10**9) + " ms")


            for j in range(1,len(photons)):
                if photons[j] < 0:
                    print(line + " is neg - didn't write")
                elif photons[j] == float('inf'):
                    print(line + " is inf - didn't write")
                else:
                    while photons[j] > lastsyncpulse:
                        file.write(str(channels) + "," + str(lastsyncpulse)+"\n")
                        lastsyncpulse = lastsyncpulse + int(taurep)
                        line = line + 1

                    if photons[j] - photons[j-1] <= deadtime: #leq since afterpulse photons have to go to other det too
                        channel = (channel+1)%channels #went to the other detector
                    else:
                        channel = numpy.random.randint(channels)
                    file.write(str(channel) + "," + str(int(photons[j]))+"\n")
        ##            print("wrote a line")
                    line = line + 1
                    if line%numlines==0:
                        #p.printProgressBar(sigcts,endsigcts)
                        print("time = " + str(lastwritten/10**12) + " s")
                        if lastwritten>endtime:
                            break #should never happen
                    if line%(numlines/10) == 0:
                        print(fullfilename + ": " + str(line*100/numlines) + "% - " + str(photons[j]/10**9) + " ms")
                    elif deltat > 1 and line%(numlines/100)==0:#if it's slow print more often
                        print(str(line*100/numlines) + "% - " + str(photons[j]/10**9) + " ms")
                
            #print("Lines written: " + str(line))


    
    

    print("Last photon time = " + str(int(max(lastphoton))/(10**12)) + " s")
    print("Signal counts = " + str(sigcts))
    print("Diffusion-in events = " + str(diffsIn))
    print("Diffusion-out events = " + str(diffsOut))
    print("Final emitters in = " + str(numEms))
    file.close()
    difffile.close()
    pfile = open(filepath +"Figures/"+ filedir + fullfilename +"/" + "params-"+ fullfilename + suffix, 'w')
    pfile.write(fullfilename +" - date and time: " + str(rt.time()) + " \n")
    
    
    pfile.write("Definitive parameters:")
    lvars = locals()
    for i in lvars:
        pfile.write(i + " : " + str(lvars[i]) + "\n")
    pfile.write("\n")
    pfile.write("Last photon time = " + str(int(max(lastphoton))/(10**12)) + " s"+" \n")
    pfile.write("Last photon in file = " + str(photons[-1]) + " ps \n")
    pfile.write("Signal counts = " + str(sigcts)+" \n")
    pfile.write("Excitations per pulse " + str(exsperpulse) +" \n")
    pfile.write("Excitations per sec " + str(AvgExEvs) +" \n")
    pfile.write("Probability of excitation per pulse " + str(probex) +" \n")
    pfile.write("Diffusion-in events = " + str(diffsIn)+" \n")
    pfile.write("Diffusion-out events = " + str(diffsOut)+ ", " +str(ndiffsOut) + " \n")
    pfile.write("Final emitters in = " + str(numEms) + " \n")
    pfile.write("Intended Avg Ems = " + str(AvgEms) + " \n")
    pfile.write("Number of records = " + str(numlines)+" \n")
    pfile.write("k_emission = " + str(int(k_emission))+" \n")
    pfile.write("k_excitation = " +str(int(k_excitation))+" \n")
    pfile.write("absXsec = " + str(absXsec) + " \n")
    pfile.write("darkcounts = " + str(darkcounts) + " Hz  \n")
    pfile.write("sensitivity = " + str(sensitivity) + " \n")
    pfile.write("deadtime = " + str(deadtime/1000) + " ns  \n")
    pfile.write("afterpulsing = " + str(afterpulse) +" \n")
    pfile.write("radius = " + str(r)+" M \n")
    pfile.write("concentration = " + str(concentration)+" M \n")
    pfile.write("Solvent viscosity = " + str(eta)+"  \n")
    pfile.write("Solvent refractive index = " + str(n)+"  \n")
    pfile.write("Focal length = " + str(foclen)+"  \n")
    pfile.write("NA = " + str(NA)+"  \n")
    pfile.write("Laser power = " + str(laserpwr) + " mW \n")
    pfile.write("Time per round: " + str(timestep) + " ps  \n")
##    pfile.write("Time limits for photon_gn = " + str(2**gnpwr) + "ps  \n")
##    pfile.write("Pulse bins = " +str(pulsebins) + " \n")
    pfile.write("Average time to diffuse out of " + str(focalVol) + " nm^3 focal volume is " + str(diffouttime/(10**9)) + "  \n")
    pfile.write("Antibunch? " + str(antibunch)+" \n")
    pfile.write("Diffusing as t = -diffouttime*numpy.log(1-numpy.random.rand()) #poissonian with average time diffouttime  \n")
    
    pfile.close()