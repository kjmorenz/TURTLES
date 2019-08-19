import scipy
import numpy

import diffuse as d

def newctsphoton(antibunch, k_emission, k_excitation):
    if antibunch == 1:
        add = scipy.random.exponential(k_emission) #+ scipy.random.exponential(k_excitation)
    #elif antibunch == 2:
    #    add = scipy.random.exponential(k_emission+k_excitation)
    else:
        add = - (k_excitation*numpy.log(1-numpy.random.rand())) #poisson

    return add

def newpulsephoton(antibunch, k_emission, k_excitation):
    if antibunch == 1:
        add = scipy.random.exponential(k_emission)
    else:
        add = - (k_excitation*numpy.log(1-numpy.random.rand())) #poisson
        add = numpy.absolute(numpy.random.normal())
    return add

def quicksort(array):
    less = []
    pivot = []
    more = []
    if len(array) <= 1:
        return array
    else:
        pivotpt = array[int(len(array)/2)]
        for i in array:
            if i < pivotpt:
                less.append(i)
            elif i > pivotpt:
                more.append(i)
            else:
                pivot.append(i)
        less = quicksort(less)
        more = quicksort(more)
    return less + pivot + more
#recursion depth limit is around 1000

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
       
    
def nextphoton(lastphoton, sensitivity,
               k_emission, k_excitation, diffouttime,
               numEms, nextOut, nextIn, endround,
               antibunch, pulsed, taurep,bfrate,
               absXsec, photonsperpulse, AvgEms, diffsIn, 
               diffsOut, ndiffsOut,countbright, probex=1):
    nextphoton = []
    if pulsed == 1:
        newphoton = newpulsephoton
    else:
        newphoton = newctsphoton
    

    while nextOut < endround: #calculate emissions before it diffuses out
        diffOut = numpy.random.randint(numEms)#index identity of emitter which diffused out
        nextem = lastphoton[diffOut]
        while nextem < nextOut:
            bfem = 0
            if pulsed == 1:
                nextem = int(nextem/taurep)*taurep + numpy.random.geometric(p=probex)*taurep #geometric distribution gives us number of pulses it took to get excited
            
            else:
                nextem = nextem + scipy.random.exponential(k_excitation)
                #print(nextem)
            #randomnumber = numpy.random.rand()
#regular emission
            nextem = nextem + newphoton(antibunch, k_emission, k_excitation)
            
            if numpy.random.rand() < sensitivity and nextem < nextOut: # if we didn't miss it due to sensitivity
                nextphoton = insert(nextphoton, nextem)

            if numpy.random.rand() < bfrate: #this emitter had bright fission occur on this excitation
                bfem = nextem + newphoton(antibunch, k_emission, k_excitation) 
                    
                if numpy.random.rand() < sensitivity and bfem < nextOut:
                        nextphoton = insert(nextphoton, bfem)
                        countbright = countbright + 1

            
            nextem = max(bfem, nextem)

        numEms = numEms - 1
        del lastphoton[diffOut]
        diffsOut = diffsOut + 1
        nextOut = nextOut + d.diffuse(diffouttime, numEms)
        if numEms == 0:
            break #we'll want to start a new round
        
        


    while nextIn < endround:
        lastphoton.append(nextIn)
        diffsIn = diffsIn + 1
        numEms = numEms + 1
        nextIn = nextIn + d.diffuse(diffouttime, AvgEms)
        if numEms == 1:
            endround = nextIn
            break #skip to next round if we just diffused out and back in
        
    for i in range(numEms):
        nextem = lastphoton[i]
        bfem = 0
        while nextem < endround: #continue adding photons until end of round, on average 10 emissions
            nextem = max(bfem, nextem)
            bfem = 0
            if pulsed == 1:
                nextem = int(nextem/taurep)*taurep + numpy.random.geometric(p=probex)*taurep #geometric distribution gives us number of pulses it took to get excited
       
            else:
                nextem = nextem + scipy.random.exponential(k_excitation)
            
            nextem = nextem + newphoton(antibunch, k_emission, k_excitation)

            if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                nextphoton = insert(nextphoton, nextem)
            
            if numpy.random.rand() < bfrate: #this emitter had bright fission occur on this excitation
                bfem = nextem + newphoton(antibunch, k_emission, k_excitation) 
                    
                if numpy.random.rand() < sensitivity: # if we didn't miss it due to sensitivity
                    nextphoton = insert(nextphoton, bfem)
                    countbright = countbright + 1

            #regular emission
            
            
        lastphoton[i] = nextem

    if numEms == 0:
        numEms = 1
        while nextOut<nextIn:
            nextOut = nextOut + d.diffuse(diffouttime, 1)
            ndiffsOut = ndiffsOut + 1
        lastphoton = [nextIn]
        nextIn = nextIn + d.diffuse(diffouttime, AvgEms)
        diffsIn = diffsIn + 1
            
    return nextphoton, lastphoton, numEms, nextIn, nextOut, diffsIn, diffsOut, ndiffsOut, countbright

def darkcounts(avtime, length, lastdc, deadtime):
    dc = []
    while len(dc) < length:
        new = lastdc -((avtime)/2)*(numpy.log((1-numpy.random.rand())))#so that the interval is (0,1] instead of [0,1)
        #divide avtime by 2 because assuming 2 detectors
        if not new - lastdc < deadtime:
            dc.append(new)
        elif numpy.random.rand() < 0.5:#if the darkcounts appeared within the deadtime of each other but on different detectors
            dc.append(new)
        lastdc = new
    return dc
