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

class SimData(ctypes.Structure):
    _fields_ = [("q", ctypes.c_void_p),
                ("boxsize", ctypes.c_void_p),
                ("flags", ctypes.c_int),
                ("N", ctypes.c_int),
                ("Nmax", ctypes.c_int),
                ("ndim", ctypes.c_int),
                ("beta", ctypes.c_double),
                ]
SimData_p = ctypes.POINTER(SimData)

dragunov_c = numpy.ctypeslib.load_library('dragunov_c', '.')
dragunov_c.eij.restype = ctypes.c_double
dragunov_c.eij.argtypes = (ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_double, )
dragunov_c.energy_i_4c.restype = ctypes.c_double
#dragunov_c.energy_i_4c.argtypes = (ctypes.c_void_p, # qdata_p
#                                   ctypes.c_int,    # N
#                                   ctypes.c_int,    # boxsize
#                                   ctypes.c_void_p, # boxsize_p
#                                   ctypes.c_int,    # flags
#                                   )
dragunov_c.energy_i_4c.argtypes = (SimData_p, # qdata_p
                                   ctypes.c_int,    # i
                                   ctypes.c_int,    # flags
                                   )

c_eij = dragunov_c.eij
c_energy_i_4c = dragunov_c.energy_i_4c

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
        SD.q =   self.q.ctypes.data
        self.ei = numpy.zeros(shape=(self.Nmax, ), dtype=numpy.float)
        self.pairlist = numpy.zeros(shape=(self.Nmax, 15), dtype=numpy.int)
        self.boxsize = numpy.asarray((10, 10, 10), dtype=numpy.float)
        SD.boxsize =   self.boxsize.ctypes.data

    def fill(self):
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
        E = 0
        q = self.q
        for i in range(self.N):
            E+= self.energy_i(i, all=False)
        return E
    def energy_fromi(self):
        return sum(self.ei[1:self.N])
    energy = energy_fromi
    def energy_i_2(self, i, all=True):
        # modified to find d forthe whole thing at once
        E = 0
        q = self.q
        boxsize = self.boxsize

        startat = 0
        if not all:
            startat = i+1

        q1 = q[i]
        d = q - q1
        d = numpy.abs(d - (numpy.floor(d/boxsize + .5)) * boxsize)
        d = numpy.sqrt((d*d).sum(axis=1))  # distances from i to j

        for j in range(startat, self.N):
            if i == j:
                continue
            E += self.eij(0, 0, d[j])
        return E
    def energy_i_4(self, i, all=True):
        #print self.q[[1,101]]
        flags = 0
        if all:
            flags = 1
        E = c_energy_i_4c(self.SD_p, i, flags)
        return E
    energy_i = energy_i_4

    def eij(self, i, j, dij):
        if dij < 1:
            return float("inf")
        elif dij < .7:
            return .7 - d
        else:
            return 0

    def trialMove(self, verbose=False):
        #print
        #print
        #print
        #print
        #self.checkContiguous()
        i = int(math.floor(random.random()*self.N))
        #print i
        qi_old = self.q[i].copy()
        #Eold = self.energy_i(i)
        Eold = self.ei[i]

        # v randn returns normal gaussion distributed points
        vec = numpy.random.randn(3) / 2
        self.q[i] = qi_old + vec

        Enew = self.energy_i(i)
        #print "Eold:", Eold, "Enew:", Enew
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
        eprt = "%5.3f -> %5.3f"%(Eold, Enew)
        if accept:
            self.ei[i] = Enew
            self.naccept += 1
            if verbose:
                print "accept move", eprt, \
                      "   ", round(self.naccept/self.ntry, 4)
        else:
            self.q[i] = qi_old
            if verbose:
                print "reject move", eprt, \
                      "   ", round(self.naccept/self.ntry, 4)
        #print "  ", "qiold:", qi_old, "qinew", self.q[i]
        #Etotal = self.energy()
        #if Etotal == float("inf"):
        #    print "Etotal:", Etotal
        #    self.display()
        #    sys.exit(3)
        

    def checkContiguous(self):
        if not self.q.flags["C_CONTIGUOUS"]:
            print "*** not C contiguous!"
            raise Exception
    def checkIntersected(self):
        all = 0
        for i in xrange(self.N):
            E = c_energy_i_4c(self.SD_p, i, all|2)
    def checkEnergyConsistency(self):
        E_fromall = self.energy_fromall()
        E_fromi = self.energy_fromi()
        deltaE = E_fromall - E_fromi
        print "E difference:", repr(deltaE)
        assert abs(E_fromall - E_fromi) < .001
    def qWrapped(self):
        q = self.q
        return numpy.mod(q, self.boxsize)
    def display(self):
        q = self.qWrapped()
        if not hasattr(self, "_display"):
            display = [ ]
            for i in range(self.N):
                display.append(visual.sphere(pos=q[i], radius=.5))
            self._display = display
        else:
            display = self._display
            for i in range(self.N):
                display[i].pos = q[i]
    def makebox(self):
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
        for i in range(100000):
            S.trialMove()
        S.checkEnergyConsistency()
    else:
        S.makebox()
        visual.scene.center = (5, 5, 5)
        i = 0
        while True:
            print i,
            if i%10000 == 0:
                S.display()
                S.checkIntersected()
                #time.sleep(2)
            i+= 1
            S.trialMove(verbose=True)
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
