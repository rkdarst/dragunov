# Richard Darst, 2007
# The University of Texas at Austin
from __future__ import division


import math
import random
random.seed(12345)
import sys
import time
import ctypes
import numpy
numpy.random.seed(123456)
import visual
try:
    from rkddp.interact import interact
except ImportError:
    pass

SVD_ENERGYI_PARTIAL = 1
SVD_VERBOSE_1       = 2
SVD_VERBOSE_2       = 4


class SimData(ctypes.Structure):
    _fields_ = [("q", ctypes.c_void_p),
                ("boxsize", ctypes.c_void_p),
                ("flags", ctypes.c_int),
                ("N", ctypes.c_int),
                ("Nmax", ctypes.c_int),
                ("ndim", ctypes.c_int),
                ("beta", ctypes.c_double),
                ("atomtypes", ctypes.c_void_p),
                ("ei", ctypes.c_void_p),
                ]
SimData_p = ctypes.POINTER(SimData)

dragunov_c = numpy.ctypeslib.load_library('dragunov_c', '.')
dragunov_c.eij.restype = ctypes.c_double
dragunov_c.eij.argtypes = (ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_double, )
dragunov_c.energy_i_4c.restype = ctypes.c_double
dragunov_c.energy_i_4c.argtypes = (SimData_p,       # qdata_p
                                   ctypes.c_int,    # i
                                   ctypes.c_int)    # flags
dragunov_c.trialMove.restype = ctypes.c_int
dragunov_c.trialMove.argtypes = (SimData_p,       # qdata_p
                                 ctypes.c_int)    # n, number of moves

c_eij = dragunov_c.eij
c_energy_i_4c = dragunov_c.energy_i_4c

dragunov_c.init_gen_rand.restype = None
print dragunov_c.init_mt(10)
#print dragunov_c.gen_rand(10)
#print dragunov_c.init_gen_rand(10)
#sys.exit()

class System:    
    def __init__(self):
        SD =      self.simData = SimData()
        self.SD_p = ctypes.pointer(SD)
        SD.N    =   self.N    =   250
        SD.Nmax =   self.Nmax =   250
        SD.beta =   self.beta =   1.

        self.naccept = 0
        self.ntry = 0

        self.q = numpy.zeros(shape=(self.Nmax, 3), dtype=numpy.float)
        SD.q   = self.q.ctypes.data
        self.ei = numpy.zeros(shape=(self.Nmax, ), dtype=numpy.float)
        SD.ei   = self.ei.ctypes.data
        self.pairlist = numpy.zeros(shape=(self.Nmax, 15), dtype=numpy.int)
        self.boxsize = numpy.asarray((10, 10, 10), dtype=numpy.float)
        SD.boxsize   = self.boxsize.ctypes.data
        self.atomtypes = numpy.zeros(shape=(self.Nmax), dtype=numpy.int)
        SD.atomtypes   = self.atomtypes.ctypes.data
        self.atomtypes[0:self.N] = 0   # default to zero

    def fill(self):
        """Fill the box with a crystal structure"""
        for i in range(self.N):
            x = i % 10
            y = (i-x) % 100 // 10
            z = i // 100
            #print (x,y,z)
            #time.sleep(.1)
            self.q[i] = (x,y,z)
        self.E = self.energy()
        for i in range(self.N):
            self.ei[i] = self.energy_i(i)
    def energy_fromall(self):
        """Return total energy by evaluating energy for all atoms"""
        E = 0
        q = self.q
        for i in range(self.N):
            E+= self.energy_i(i, partial=True)
        return E
    def energy_fromi(self):
        """Return total energy by using sum of cached values for each atom"""
        return .5 * sum(self.ei[1:self.N])
    energy = energy_fromi
    def energy_i_2(self, i, partial=False):
        """Return energy of atom i"""
        E = 0
        q = self.q
        boxsize = self.boxsize
        eij = self.eij
        atomtypes = self.atomtypes

        startat = 0
        if partial:
            startat = i+1

        q1 = q[i]
        d = q - q1
        d = numpy.abs(d - (numpy.floor(d/boxsize + .5)) * boxsize)
        d = numpy.sqrt((d*d).sum(axis=1))  # distances from i to j

        for j in range(startat, self.N):
            if i == j:
                continue
            E += eij(atomtypes[i], atomtypes[j], d[j])
        return E
    def energy_i_4(self, i, partial=False):
        """Return energy of atom i, using c function"""
        flags = 0
        if partial:
            flags |= SVD_ENERGYI_PARTIAL
        E = c_energy_i_4c(self.SD_p, i, flags)
        return E
    energy_i = energy_i_4 # select which energ eval function to use

    def eij(self, i, j, dij):
        """Energy of interaction between atomtypes i and j at distance dij"""
        if dij < 1:
            return float("inf")
        elif dij < .7:
            return .7 - d
        else:
            return 0

    def trialMove_py(self, verbose=False):
        """Try a move using metropolis criteria"""







        i = int(math.floor(random.random()*self.N))

        qi_old = self.q[i].copy()


        Eold = self.ei[i]

        # v randn returns normal gaussion distributed points
        vec = numpy.random.randn(3) / 2
        self.q[i] = qi_old + vec


        Enew = self.energy_i(i)

        if Enew <= Eold:     # always accept, E decreases
            accept = True
        else:                # accept with prob


            x = math.exp(self.beta*(Eold-Enew))
            ran = random.random()
            if  ran < x:
                accept = True
            else:
                accept = False

        self.ntry += 1
        if accept:

            self.ei[i] = Enew
            self.naccept += 1

        else:

            self.q[i] = qi_old


            
        if verbose:
            print "%6d %d, %5.3f->%5.3f accF: %.4f\n"%(
                self.n, int(accept), Eold, Enew, self.naccept/self.ntry)


    def trialMove_c(self, n=1, verbose=False):
        """Try n moves using metropolis criteria, using C function"""
        # flags: 4 = verbose
        naccept = dragunov_c.trialMove(self.SD_p, n)
        self.ntry += n
        self.naccept += naccept
        if verbose:
            print "%10s moves"%self.ntry, round(self.naccept/self.ntry, 4)
        
    trialMove = trialMove_c

    def checkContiguous(self):
        """Check if the numpy array are contiguous, if not raise an error."""
        # should be modified to check all arrays
        if not self.q.flags["C_CONTIGUOUS"]:
            print "*** not C contiguous!"
            raise Exception
    def checkIntersected(self):
        # self-test function (not used)
        partial = 1
        for i in xrange(self.N):
            E = c_energy_i_4c(self.SD_p, i, partial|2)
    def checkEnergyConsistency(self):
        """
        """
        E_fromall = self.energy_fromall()
        E_fromi = self.energy_fromi()
        deltaE = E_fromall - E_fromi
        print "E_fromall:",E_fromall, "E_fromi",E_fromi, \
              "E difference:", repr(deltaE)
        assert abs(E_fromall - E_fromi) < .001
    def qWrapped(self):
        """Return coordinates wrapped to fit in boxsize"""
        # wrapped means in interval [0, boxsize)
        q = self.q
        return numpy.mod(q, self.boxsize)
    def display(self):
        """Update visual display of atoms"""
        q = self.qWrapped()
        c = {0: visual.color.white,
             1: visual.color.green,
             }
        if not hasattr(self, "_display"):
            display = [ ]
            for i in range(self.N):
                display.append(visual.sphere(pos=q[i], radius=.5,
                                             color=c[self.atomtypes[i]]))
            self._display = display
        else:
            display = self._display
            for i in range(self.N):
                display[i].pos = q[i]
    def makebox(self):
        """Put the simulation box on the visual display"""
        visual.scene.center = (5, 5, 5)
        radius = .02
        x,y,z = self.boxsize
        c = visual.color.blue
        visual.cylinder(pos=(0,0,0), axis=(x, 0, 0), radius=radius, color=c)
        visual.cylinder(pos=(0,y,0), axis=(x, 0, 0), radius=radius, color=c)
        visual.cylinder(pos=(0,0,z), axis=(x, 0, 0), radius=radius, color=c)
        visual.cylinder(pos=(0,y,z), axis=(x, 0, 0), radius=radius, color=c)

        visual.cylinder(pos=(0,0,0), axis=(0, y, 0), radius=radius, color=c)
        visual.cylinder(pos=(x,0,0), axis=(0, y, 0), radius=radius, color=c)
        visual.cylinder(pos=(0,0,z), axis=(0, y, 0), radius=radius, color=c)
        visual.cylinder(pos=(x,0,z), axis=(0, y, 0), radius=radius, color=c)

        visual.cylinder(pos=(0,0,0), axis=(0, 0, z), radius=radius, color=c)
        visual.cylinder(pos=(x,0,0), axis=(0, 0, z), radius=radius, color=c)
        visual.cylinder(pos=(0,y,0), axis=(0, 0, z), radius=radius, color=c)
        visual.cylinder(pos=(x,y,0), axis=(0, 0, z), radius=radius, color=c)
        
if __name__ == "__main__":
    S = System()
    S.fill()
    if len(sys.argv) > 1 and sys.argv[1] == "timing":
        for i in xrange(10000):
            S.trialMove(n=10)
        #S.trialMove_c(100000)
        S.checkEnergyConsistency()
    elif len(sys.argv) > 1 and sys.argv[1] == "demo":
        S.makebox()
        i = 0
        while True:
            S.display()
            S.trialMove()
            if i%100 == 0:
                print "\r%9d"%i,
                sys.stdout.flush()
            i += 1
        S.checkEnergyConsistency()
    else:
        i = 0
        S.makebox()
        while True:
            #print i,
            if i%10000 == 0:
                S.display()
                #S.checkIntersected()
                #time.sleep(2)
            i+= 1
            S.trialMove(verbose=False)
            if i%1000000 == 0:
                fo = file("output.txt", "w")
                fo.write("title equilibrated 1A hard spheres\n")
                fo.write("natoms 250\n")
                fo.write("boxsize 10 10 10\n")
                fo.write("step %s\n"%i)
                fo.write("data\n")
                for j in range(S.N):
                    fo.write("%r %r %r\n"%tuple(S.q[j]))
                del fo
    #time.sleep(10)
