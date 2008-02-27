# Richard Darst, 2008
# The University of Texas at Austin

"""Tests for pressure of the LJ fluid.
"""

import sys
import dragunov

skip = 100
Nsteps = 1000000


tests = (.1, .2, .3, .4, .5, .6, .7, .8)

for density in tests:
    print "density:", density
    N = int(density * 1000)
    logname = "tests/lj-pressure/"+"logfile-density=%s.txt"%density
    logfile = file(logname, "w")

    S = dragunov.System(N=N, beta=1/2.0,
                        boxsize=(10,10,10), trialMoveScale=.25,
                        dt=.001)
    pairlist = True
    #if pairlist: S.flags |= dragunov.SVD_USE_PAIRLIST
    S.fillRandom()
    #S.fill()
    S.removeOverlaps(10.)
    #S.zeroVelocity()
    #S.pairlistInit()
    i = 0
    S.makebox()
    S.display()
    S.trialMove(n=10000)
    S.display()
    #time.sleep(5)
    #logfile=sys.stdout
    print >> logfile, "# i, T, rho, E, P, P_avg"
    while i <= Nsteps:
        #if i%10 and pairlist: S.pairlistInit(3.25)
        i += skip
        if i == 150000:
            S.resetStatistics()
        S.trialMove(verbose=False, n=skip)
    
        #ke = dragunov.c_mdStep(S.SD_p, skip, 0)
        #if 0 or pairlist: S.pairlistCheck(2.75)
        #S.widomInsert()
        if i % 100 == 0:
            S.display()
            print >> logfile, i, S.T, S.density, S.energy(), \
                  S.pressure(add=True), S.pressureAverage() #, ke
        print "\rstep: %d"%i,
        sys.stdout.flush()
        
        #print S.pairlist
        #S.pairlistCheck(strict=False)
        logfile.flush()
        #time.sleep(.5)
        #print
