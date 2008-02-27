# Richard Darst, 2008
# The University of Texas at Austin

"""Tests for enengy conservation.
"""

import sys
import dragunov

skip = 1
N = 200
#N = 2
S = dragunov.System(N=N, beta=1/2.0,
                    boxsize=(10,10,10), trialMoveScale=.25,
                    dt=.001)
pairlist = False
if pairlist: S.flags |= SVD_USE_PAIRLIST
S.fillRandom()
#S.fill()
S.removeOverlaps(10.)
S.zeroVelocity()
#S.atomtypes[numpy.asarray((numpy.random.uniform(1500, size=500)),
#                          dtype=int)] = 1
#S.pairlistInit()
i = 0
S.makebox()
S.display()
#time.sleep(5)
#logfile = file("logfile.txt", "w")
logfile = file("/dev/stdout", "w")
#logfile=sys.stdout
print >> logfile, "# i, T, rho, E, P, P_avg"
while True:
    if pairlist: S.pairlistInit(3.25)
    i += skip
    if i == 150000:
        S.resetStatistics()
    #S.trialMove(verbose=False, n=skip)

    ke = dragunov.c_mdStep(S.SD_p, skip, 0)
    if pairlist: S.pairlistCheck(2.75)
    #S.widomInsert()
    if i % 100 == 0:
        S.display()
        print >> logfile, i, S.T, S.density, S.energy(), \
              S.pressure(add=True), S.pressureAverage(), ke
    print "\rstep: %d"%i,
    sys.stdout.flush()
    
    #print S.pairlist
    #S.pairlistCheck(strict=False)
    logfile.flush()
    #time.sleep(.5)
    #print
