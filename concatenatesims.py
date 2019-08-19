import os
from analyze2 import analyze
from makefig2 import makeafig 
filepath = "C:/Users/Karen/Dropbox (WilsonLab)/WilsonLab Team Folder/Data/Karen/Simulation-cleanedMar92018/"
filedir = "2019-06-06-ACseries/"
filedir = "2019-06-06-concseries/"
filedir = "2019-05-30-redotento9longtaufewerems/"
filedir = "2019-07-09/"
fname = ["Newwrite10to5tau0pmqy-ac", "Newwrite10to5tau1pmqy-ac"]
fnames = ["4000", "3500", "3400", "3300","3200", "3100","3000", "2600","2300","2000","1700","1500","1200","1000"]
filenames = ["1pmqy-sens0.001ac50shorte7","1pmqy-sens0.001ac1000shorte7","1pmqy-sens0.1ac1000shorte7","1pmqy-sens0.1ac50shorte7",
                "0pmqy-sens0.001ac1000shorte7","0pmqy-sens0.001ac1000shorte7","0pmqy-sens0.1ac1000shorte7","0pmqy-sens0.1ac50shorte7"]
filenames = ["RealisticParams1emlongbf0", "RealisticParams1emlongbf1",
            "shortlifetime1emhighsensbf1", "shortlifetime1emhighsensbf0",
            "shortlifetime1embf1", "shortlifetime1embf0",
            "RealisticParams1emSens1bf1", "RealisticParams1emSens1bf0"]
#  \\          "2019-07-09RealisticParamsbf1", "2019-07-09RealisticParamsbf0"]

#for f in range(len(fnames)):
#    filenames.append(fname[0] + fnames[f])
#    filenames.append(fname[1] + fnames[f])
'''filenames = []
concentration = 0.00002#2*10**(-5)
concentrations = [concentration * 10, concentration, concentration / 10, concentration / 100, concentration / 1000]  
for bfrate in [0,1]:
    for conc in concentrations:
        f = "Newwrite10to5tau" + str(bfrate) + "pmqy-conc" + str(conc)
        filenames.append(f)'''


suffix = ".txt"
totime = 1000000 #s

mode = "t2"
gnpwr = 35#5#ms 20#40
numbins = 256#2**8 #
pulsebins = (8)-1#should always be an odd number
channels = 2

makefig = 1
pulsed = 0


def concatenatefile(filepath, filedir, filename, suffix, totime):
    dfilepath = filepath + "RawData/" + filedir
    if not os.path.isdir(dfilepath + filename+ str(totime) + "x"):
        os.mkdir(dfilepath + filename+ str(totime) + "x")
        os.mkdir(filepath + "Figures/"+ filedir+ filename+ str(totime) + "x")

    f = open(dfilepath + filename +"/" + filename + suffix, 'r')
    t = open(dfilepath + filename + str(totime) + "x"+"/" + filename + str(totime) + "x" + suffix, 'w')

    lastphoton = f.readline()
    nextline = f.readline()
    while nextline != '':
        t.write(lastphoton)
        lastphoton = nextline
        nextline = f.readline()
    t.write(lastphoton)
    lastphotont = int(lastphoton[2:-1])
    lastwritten = int(lastphoton[2:-1])
    f.seek(0)
    nextline = f.readline()
    nextline = nextline[0:2]+str(int(nextline[2:-1]) + lastphotont) + '\n'
    r = 1

    while lastwritten < totime*10**12:
        while nextline != '':
            t.write(nextline)
            nextline = f.readline()
            if nextline == '':
                break
            lastwritten = int(nextline[2:-1]) + lastphotont*r
            if lastwritten > totime*10**12:
                break
            nextline = nextline[0:2]+str(lastwritten) + '\n'
        r = r + 1
        f.seek(0)
        nextline = f.readline()
        nextline = nextline[0:2]+str(int(nextline[2:-1]) + lastphotont*r) + '\n'

    f.close()
    t.close()

count = 0
for filename in filenames:
    
    if count > 1:
        concatenatefile(filepath, filedir, filename, suffix, totime)
        analyze(filepath, filedir, filename +str(totime) + "x", 10**7, 2, mode, gnpwr, 
            numbins, pulsebins, channels, makefig, makeafig, pulsed, foutedit = "long")#needs to get set up to run on greg data....
    count = 2
    analyze(filepath, filedir, filename +str(totime) + "x", 10**7, 2, mode, 35, 
            numbins, pulsebins, channels, makefig, makeafig, pulsed, foutedit = "normal", 
            normalize = 1, dopic = 0)