# Richard Darst, 2008
# The University of Texas at Austin

"""Tests for pressure of the LJ fluid.
"""

import random
import sys
import dragunov

N = 400
skip = 100
Nsteps = 1000000
isobarPressure = 1.0
Ps = (.25, .5, 1., 2., 4.)
Ps = (1., )


for isobaricPressure in Ps:
    logname = "tests/isobaric/logfile-P=%s.txt"%isobaricPressure
    logname = "logfile-isobaric.txt"
    logfile = file(logname, "w")
    
    i = 0
    S = dragunov.System(N=N, beta=1/2.0,
                        boxsize=(10,10,10), trialMoveScale=.25,
                        isobaricPressure=isobaricPressure,
                        forceField="lennardjones")
    #S.setMoveProb(shift=N, pressure=1)  # done automatically now
    pairlist = False
    S.fillRandom()
    S.removeOverlaps(10.)
    
    S.makebox()
    S.display()
    
    S.trialMove(n=10000)
    S.display()
    #                    1  2  3  4    5  6  7     8  9     10
    print >> logfile, "# i, T, N, rho, E, P, Pavg, V, Vavg, rho_avg"
    while i <= Nsteps:
        if i%1000 == 0 and pairlist:
            S.pairlistCheck(2.75)
            S.pairlistInit(3.25)
        i += skip
        if i == 250000:
            S.resetStatistics()
        # The below is intergrated into the C trialMove now.
        #if random.random() * S.N <= 1:
        #    S.isobaricTrialMove_py(pressure=1.0, lnVScale=.5)
        S.trialMove(verbose=False, n=skip)
    
        if i % 100 == 0:
            if i%1000 == 0: S.display()
            print >> logfile, i, S.T, S.N, S.density, S.energy(), \
                  S.pressure(add=True), S.pressureAverage(), \
                  S.volume, S.volumeAverage(), \
                  S.N / S.volumeAverage()
                  #, ke
        print "step: %d"%i,"\r",
        sys.stdout.flush()
        
        #print S.pairlist
        logfile.flush()
        #time.sleep(.5)
        #print
