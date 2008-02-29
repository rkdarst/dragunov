# Richard Darst, 2008
# The University of Texas at Austin

"""Tests for pressure of the LJ fluid.
"""

import random
import sys
import dragunov

skip = 100
Nsteps = 10000000
densities = (.6,)


for density in densities:
    #logname = "logfile.txt"
    logname = "tests/widomNVT/logfile-d=%s.txt"%density
    logfile = file(logname, "w")
    
    i = 0
    N = int(1000 * density)
    S = dragunov.System(N=N, beta=1/2.0,
                        boxsize=(10,10,10), trialMoveScale=.25,
                        dt=.001,
                        )
    S.fillRandom()
    S.removeOverlaps(10.)
    
    S.makebox()
    S.display()
    
    #S.trialMove(n=10000)
    S.display()
    #                    1  2  3  4    5  6  7     8  9     10       11
    print >> logfile, "# i, T, N, rho, E, P, Pavg, V, Vavg, rho_avg, mu"
    while i <= Nsteps:
        #if i%10 and pairlist: S.pairlistInit(3.25)
        i += skip
        if i == 100000:
            S.resetStatistics()
        S.trialMove(verbose=False, n=skip)
    
        #if 0 or pairlist: S.pairlistCheck(2.75)
        S.widomInsert(n=10)
        if i % 100 == 0:
            if i%1000 == 0: S.display()
            print >> logfile, i, S.T, S.N, S.density, S.energy(), \
                  S.pressure(add=True), S.pressureAverage(), \
                  S.volume, S.volumeAverage(), \
                  S.N / S.volumeAverage(), \
                  S.widomInsertResults()
        print "step: %d"%i,"\r",
        sys.stdout.flush()
        
        #print S.pairlist
        #S.pairlistCheck(strict=False)
        logfile.flush()
        #time.sleep(.5)
        #print
