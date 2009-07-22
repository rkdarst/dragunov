# Richard Darst, 2008
# The University of Texas at Austin

"""Tests for pressure of the LJ fluid.
"""

import sys
import dragunov

Nsteps = 20000
N = 500
T = 2.0

tests = (.1, .2, .3, .4, .5, .6, .7, .8)
tests = (.6, )

for density in tests:
    #logname = "tests/lj-pressure/"+"logfile-density=%s.txt"%density
    logname = "logfile-lj.txt"
    logfile = file(logname, "w")

    boxedge = (N / density)**(1/3.)
    S = dragunov.System(N=N, beta=1/T,
                        boxsize=(boxedge, boxedge, boxedge),
                        trialMoveScale=.25,
                        isobaricPressure=1.75,
                        forceField="lennardjones")
    pairlist = True
    S.fillRandom()
    #S.fill()
    S.removeOverlaps(10.)

    S.pairlistEnable(3.5, 2.75, strict=True)
    S.makebox()
    S.display()
    S.trialMove(n=10)
    S.resetTime()
    S.display()

    S.loginfo(logfile)
    while S.mctime < Nsteps:
        if S.mctime in (200, 2000, ):
            S.resetStatistics()
        S.trialMove(verbose=False, n=1)
    
        S.display()
        S.log(logfile)
        print "\rstep: %d"%S.mctime,
        sys.stdout.flush()
        
        logfile.flush()
