# Richard Darst, 2008
# The University of Texas at Austin

"""Tests for pressure of the LJ fluid.
"""

import random
import sys
import dragunov

skip = 100
Nsteps = 1000000
i = 0
isobarPressure = 1.0

N = 400 #int(density * 1000)
logname = "logfile.txt"
logfile = file(logname, "w")

S = dragunov.System(N=N, beta=1/2.0,
                    boxsize=(10,10,10), trialMoveScale=.25,
                    dt=.001,
                    isobarPressure=isobarPressure)
#S.setMoveProb(shift=N, pressure=1)
pairlist = False

S.fillRandom()
S.removeOverlaps(10.)

S.makebox()
S.display()

S.trialMove(n=1000)
S.display()

print >> logfile, "# i, T, N, V, rho, E, P, P_avg"
while i <= Nsteps:
    #if i%10 and pairlist: S.pairlistInit(3.25)
    i += skip
    if i == 150000:
        S.resetStatistics()
    #if random.random() * S.N <= 1:
    #    S.isobaricTrialMove_py(pressure=1.0, lnVScale=.5)
    S.trialMove(verbose=False, n=skip)

    #ke = dragunov.c_mdStep(S.SD_p, skip, 0)
    #if 0 or pairlist: S.pairlistCheck(2.75)
    #S.widomInsert()
    if i % 100 == 0:
        S.display()
        print >> logfile, i, S.T, S.N, S.volume, S.density, S.energy(), \
              S.pressure(add=True), S.pressureAverage() #, ke
    print "step: %d"%i,"\r",
    sys.stdout.flush()
    
    #print S.pairlist
    #S.pairlistCheck(strict=False)
    logfile.flush()
    #time.sleep(.5)
    #print
