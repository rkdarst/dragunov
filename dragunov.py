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
import cliargs

time.sleep(0)  # prevent an "imported module not used" error


SVD_ENERGYI_PARTIAL       =   1
SVD_VERBOSE_1             =   2
SVD_VERBOSE_2             =   4
SVD_PAIRLIST_INCREMENTAL  =   8
SVD_USE_PAIRLIST          =  16


class SimData(ctypes.Structure):
    _fields_ = [("q", ctypes.c_void_p),
                ("qold", ctypes.c_void_p),
                ("force", ctypes.c_void_p),
                ("boxsize", ctypes.c_void_p),
                ("flags", ctypes.c_int),
                ("N", ctypes.c_int),
                ("Nmax", ctypes.c_int),
                ("ndim", ctypes.c_int),
                ("dt", ctypes.c_double),
                ("beta", ctypes.c_double),
                ("atomtypes", ctypes.c_void_p),
                ("ei", ctypes.c_void_p),
                ("pairlist", ctypes.c_void_p),
                ("trialMoveScale", ctypes.c_double),
                ("pairlist_minDistance", ctypes.c_double),
                ]
SimData_p = ctypes.POINTER(SimData)

dragunov_c = numpy.ctypeslib.load_library('dragunov_c', '.')

c_eij = dragunov_c.eij
dragunov_c.eij.restype = ctypes.c_double
dragunov_c.eij.argtypes = (ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_double, )

c_energy_i = dragunov_c.energy_i
c_energy_i.restype = ctypes.c_double
c_energy_i.argtypes = (SimData_p,       # qdata_p
                       ctypes.c_int,    # i
                       ctypes.c_int)    # flags

c_energy = dragunov_c.energy
c_energy.restype = ctypes.c_double
c_energy.argtypes = (SimData_p,       # qdata_p
                     ctypes.c_int)    # flags

c_integrate = dragunov_c.integrate
c_integrate.restype = ctypes.c_double
c_integrate.argtypes = (SimData_p,       # qdata_p
                        ctypes.c_int)    # flags

c_mdStep = dragunov_c.mdStep
c_mdStep.restype = ctypes.c_double
c_mdStep.argtypes = (SimData_p,       # qdata_p
                     ctypes.c_int,    # n, number of moves
                     ctypes.c_int)    # flags

c_trialMove = dragunov_c.trialMove
c_trialMove.restype = ctypes.c_int
c_trialMove.argtypes = (SimData_p,       # qdata_p
                        ctypes.c_int,    # n, number of moves
                        ctypes.c_int)    # flags

c_forcedotr_i = dragunov_c.forcedotr_i
c_forcedotr_i.restype = ctypes.c_double
c_forcedotr_i.argtypes = (SimData_p,       # qdata_p
                          ctypes.c_int,    # i
                          ctypes.c_int)    # flags

c_pairlistInit = dragunov_c.pairlistInit
c_pairlistInit.restype = ctypes.c_int
c_pairlistInit.argtypes = (SimData_p,       # SimData_p
                           ctypes.c_double, # cutoff
                           ctypes.c_int)    # flags

c_pairlistCheck = dragunov_c.pairlistCheck
c_pairlistCheck.restype = ctypes.c_int
c_pairlistCheck.argtypes = (SimData_p,       # SimData_p
                            ctypes.c_double, # warn
                            ctypes.c_int)    # flags


dragunov_c.init_gen_rand.restype = None
dragunov_c.init_mt(10)

class System(object):
    def _volume_get(self):
        return numpy.product(self.boxsize)
    volume = property(fget=_volume_get)

    def _density_get(self):
        return self.N / numpy.product(self.boxsize)
    density = property(fget=_density_get)

    def _temperature_get(self):
        return 1 / self.beta
    T = property(fget=_temperature_get)

    def __init__(self, N, beta=1., Nmax=None, boxsize=(10,10,10),
                 trialMoveScale=1,
                 pressureTarget=None,
                 dt=.01):
        if Nmax == None:
            Nmax = N+50
        SD =      self.SD = SimData()
        self.SD_p = ctypes.pointer(SD)
        SD.N    =   self.N    =   N
        SD.Nmax =   self.Nmax =   Nmax
        SD.dt =     self.dt   =   dt
        SD.beta =   self.beta =   beta
        self.pressureTarget = pressureTarget
        self.sizeScalar = 1.
        SD.trialMoveScale = self.trialMoveScale = trialMoveScale
        self.mu_dict = { }
        self._pressureList = [ ]

        self.naccept = 0
        self.ntry = 0

        self.q = numpy.zeros(shape=(self.Nmax, 3), dtype=numpy.float)
        SD.q   = self.q.ctypes.data
        self.qold = numpy.zeros(shape=(self.Nmax, 3), dtype=numpy.float)
        SD.qold   = self.qold.ctypes.data
        self.force = numpy.zeros(shape=(self.Nmax, 3), dtype=numpy.float)
        SD.force   = self.force.ctypes.data
        self.ei = numpy.zeros(shape=(self.Nmax, ), dtype=numpy.float)
        SD.ei   = self.ei.ctypes.data
        self.pairlist = numpy.zeros(shape=(self.Nmax, 152), dtype=numpy.int_)
        SD.pairlist   = self.pairlist.ctypes.data
        self.boxsize = numpy.asarray(boxsize, dtype=numpy.float)
        SD.boxsize   = self.boxsize.ctypes.data
        self.atomtypes = numpy.zeros(shape=(self.Nmax), dtype=numpy.int_)
        SD.atomtypes   = self.atomtypes.ctypes.data
        self.atomtypes[0:self.N] = 0   # default to zero

        self.flags = 0

    def resetStatistics(self):
        self._pressureList = [ ]


    def acceptRatio(self):
        return self.naccept / self.ntry
    def acceptRatio_last(self):
        return self.naccept_last / self.ntry_last
    def fill(self, a=10, b=10, c=10, scale=1.):
        """Fill the box with a crystal structure"""
        for i in range(self.N):
            #x = i % 10
            #y = (i-x) % 100 // 10
            #z = i // 100
            x = i % a
            y = (i // a) % b
            z = (i // (a*b)) % c
            #print (x,y,z)
            #time.sleep(.1)
            self.q[i] = (x*scale, y*scale, z*scale)
        self.E = self.energy()
        for i in range(self.N):
            self.ei[i] = self.energy_i(i)
    def fillRandom(self):
        for i in range(self.N):
            newpos = numpy.random.uniform(size=(3,)) * self.boxsize
            self.q[i] = newpos
    def zeroVelocity(self):
        self.qold[:] = self.q[:]
    def eij(self, i, j, dij):
        """Energy of interaction between atomtypes i and j at distance dij"""
        if j<i:   j,i = i,j  # make i the lower index

        if dij < 0.571428571428571:
            return float("inf")
        elif dij < 1.:
            return 1 - dij
        else:
            return 0.

    def fij(self, i, j, dij, r):
        """Energy of interaction between atomtypes i and j at distance dij"""
        if j<i:   j,i = i,j  # make i the lower index

        if dij < 0.571428571428571:
            return (0,0,0)


        elif dij < 1.:
            return r / dij




        else:
            return (0,0,0)

    def energy_fromall(self):
        """Return total energy by evaluating energy for all atoms"""
        #E = 0
        #for i in range(self.N):
        #    E += c_energy_i(self.SD_p, i, SVD_ENERGYI_PARTIAL)
        #    #E+= self.energy_i(i, partial=True)
        #return E
        return c_energy(self.SD_p, 0)

    def energy_fromi(self):
        """Return total energy by using sum of cached values for each atom"""
        raise Exception("this doesn't work")
        return .5 * sum(self.ei[0:self.N])
    energy = energy_fromall
    def energy_i_py(self, i, partial=False):
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
    def energy_i_c(self, i, partial=False):
        """Return energy of atom i, using c function"""
        flags = self.flags
        if partial:
            flags = flags | SVD_ENERGYI_PARTIAL
        E = c_energy_i(self.SD_p, i, flags)
        return E
    energy_i = energy_i_c # select which energ eval function to use

    ##  def pressure_i_py(self, i, partial=False):
    ##      """Return energy of atom i"""
    ##      fdotr = 0
    ##      q = self.q
    ##      boxsize = self.boxsize
    ##      #atomtypes = self.atomtypes
    ##  
    ##      startat = 0
    ##      if partial:
    ##          startat = i+1
    ##  
    ##      q1 = q[i]
    ##      d = q - q1
    ##      d = numpy.abs(d - (numpy.floor(d/boxsize + .5)) * boxsize)
    ##      d = numpy.sqrt((d*d).sum(axis=1))  # distances from i to j
    ##  
    ##      for j in range(startat, self.N):
    ##          if i == j:
    ##              continue
    ##          #fdotr += fij(atomtypes[i], atomtypes[j], d[j])e
    ##          # XXX
    ##      return fdotr
    def pressure_c(self, add=False):
        flags = self.flags
        flags = flags | SVD_ENERGYI_PARTIAL
        fdotr = 0
        
        for i in range(0, self.N-1):
            fdotr += c_forcedotr_i(self.SD_p, i, flags)
        volume = self.volume
        density = self.density
        dimensions = 3
                 # v-- should be "+" for LJ correct results
                 # v-- this also makes pressure increase with
                 #     increasing density for S2S- 2s units
        pressure = + fdotr / (dimensions * volume)
        if add:
            pressure += density/self.beta
        self._pressureList.append(pressure)
        return pressure
    pressure = pressure_c
    def pressureAverage(self):
        return sum(self._pressureList)/len(self._pressureList)
    def pressureTrialMove(self):
        # get current pressure
        # trial boxsize change
        # - our current size-scalar is L
        # - we need to make a reversible random move,
        #   (random number .9 to 1 and .5 chance of inverting?)
        #
        pass

    def trialMove_py(self, verbose=False):
        """Try a move using metropolis criteria"""





        i = int(math.floor(random.random()*self.N))
        qi_old = self.q[i].copy()


        #Eold = self.ei[i] # ei[i] not updated if something else moves closer
        Eold = self.energy_i(i)

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
                self.N, int(accept), Eold, Enew, self.naccept/self.ntry)

    def trialMove_c(self, n=1, verbose=False):
        """Try n moves using metropolis criteria, using C function"""
        # flags: 4 = verbose
        naccept = c_trialMove(self.SD_p, n, self.flags)
        self.ntry += n
        self.ntry_last = n
        self.naccept += naccept
        self.naccept_last = naccept
        if verbose:
            print "%10s moves"%self.ntry, round(self.naccept/self.ntry, 4)
        
    trialMove = trialMove_c
    def widomInsert(self, type=0, n=None):
        """Save a data point for widom insertion.
        """
        # n is the number of test moves to run
        if n is not None:
            for _ in xrange(n):
                self.widomInsert(type=type)
            return
        newpos = numpy.random.uniform(size=(3,)) * self.boxsize
        newi = self.addAtom(newpos, type=type)
        Eadded = self.energy_i(newi)
        #print Eadded, self.beta
        self.mu_dict.setdefault(type, []).append(math.exp(-Eadded*self.beta))
        self.delAtom(newi)
        
    def widomInsertResults(self, type=0):
        """
        """
        volume = numpy.product(self.boxsize)
        avg = sum(self.mu_dict[type]) / len(self.mu_dict[type])
        #print self.mu_dict[type]
        if avg == 0:
            return 0  # this is technically incorrect, but we need to return
                      # something that won't make parse errors
        # constant related to the deBroglie wavelength
        # the constant should also be adjusted for number of dimensions
        constant = 3./self.beta

        mu = -math.log(volume*avg/(self.N+1))/self.beta + constant
        return mu
        # 

        # constant related to the deBroglie wavelength
        # the constant should also be adjusted for number of dimensions
        #thermWavelength3 = 1.**1.5
        thermWavelength3 = 1.
        mu_ig = -math.log(volume*thermWavelength3/(self.N+1))/self.beta
        mu_excess = -math.log(avg)/self.beta
        return mu_ig + mu_excess
    def widomInsertCorrection1(self):
        # _add_ this to the reported mu tocorrect the mu.
        constant = 3./self.beta
        return math.log(self.volume/(self.N+1))/self.beta - constant
    
    mu = widomInsertResults
    def addAtom(self, pos, type):
        """Add an arbitrary atom to the system.  Return new index"""
        if self.N == self.Nmax:
            raise Exception("can't fit more atoms here")
        self.SD.N    =   self.N    =   self.N+1
        i = self.N-1
        self.q[i] = pos
        self.atomtypes[i] = type
        self.ei[i] = self.energy_i(i)
        return i
    def delAtom(self, i):
        """Delete an arbitry atom from the system."""
        N = self.N
        if i != N-1:
            raise Exception("so far we can only remove the last atom")
        self.SD.N    =   self.N    =   self.N-1
        self.q[i] = 0.
        self.atomtypes[i] = 0.
        self.ei[i] = 0.

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
            c_energy_i(self.SD_p, i, partial|2)
    def checkEnergyConsistency(self):
        """
        """
        E_fromall = self.energy_fromall()
        E_fromi = self.energy_fromi()
        deltaE = E_fromall - E_fromi
        print "E_fromall:",E_fromall, "E_fromi",E_fromi, \
              "E difference:", repr(deltaE)
        assert abs(E_fromall - E_fromi) < .001
    def pairlistInit(self, cutoff):
        flags = self.flags
        return c_pairlistInit(self.SD_p, cutoff, flags)
        
    def pairlistCheck(self, warn, strict=True):
        nviolations = c_pairlistCheck(self.SD_p, warn, self.flags)
        if nviolations:
            print "nviolatios:", nviolations, "at step", self.ntry
        if strict:
            assert nviolations == 0
        #print self.SD.pairlist_minDistance
    def removeOverlaps(self, Emax=float('inf')):
        """Remove overlaps of atoms (make energy non-infinite)

        Iterates through all atoms, and for each one with infinite
        energy, give it random displacements until it has finite
        energy.
        """
        for i in range(self.N):
            Eold = self.energy_i(i)
            #if Eold != float('inf'):
            #if Eold >= Emax:
            #    print Eold, "E too high"
            if Eold < Emax:
                # skip atoms that already have finite energy
                continue
            while Eold >= Emax:
                # until it has finite energy, give it a random displacement.
                # (this is a normal displacement, so may not be the
                #  most efficient)
                self.q[i] += numpy.random.randn(3) * 2
                Eold = self.energy_i(i)
            self.ei[i] = Eold
        print "done removing overlaps"

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
                display.append(visual.sphere(pos=q[i],
                                             radius=getattr(self, "dispDiameter", 4./7.)/2, #r, not d
                                             color=c[self.atomtypes[i]]))
            self._display = display
        else:
            display = self._display
            for i in range(self.N):
                display[i].pos = q[i]
    def makebox(self):
        """Put the simulation box on the visual display"""
        visual.scene.center = self.boxsize / 2.
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
    def writeUghfile(self, filename):
        fo = file(filename, "w")
        fo.write("title equilibrated 1A hard spheres\n")
        fo.write("natoms %s\n"%self.N)
        fo.write("boxsize %s %s %s\n"%(self.boxsize[0],
                                       self.boxsize[1],
                                       self.boxsize[2],))
        fo.write("step %s\n"%(self.ntry+1))
        fo.write("data\n")
        for j in range(self.N):
            fo.write("%r %r %r\n"%tuple(self.q[j]))
        del fo
    def readUghfile(self, filename):
        fo = file(filename)
        while True:
            line = fo.readline()
            if line == "":   # break on end of file
                break
            line = line.split()
            if line == [ ]:  # ignore empty lines
                continue
            if line[0] == "natoms":
                natoms = int(line[1])
            if line[0] == "boxsize":
                boxsize = float(line[1]), float(line[2]), float(line[3])
                if boxsize[0] != self.boxsize[0] or \
                   boxsize[1] != self.boxsize[1] or \
                   boxsize[2] != self.boxsize[2]:
                    raise Exception("wrong boxsize")
            if line[0] == "title":
                print "reading file:", " ".join(line[1:])
            if line[0] == "data":
                for i in range(natoms):
                    line = fo.readline().split()
                    self.q[i] = float(line[0]), float(line[1]), float(line[2])
            

def readUghfile(filename):
    fo = file(filename)
    while True:
        line = fo.readline()
        if line == "":   # break on end of file
            break
        line = line.split()
        if line == [ ]:  # ignore empty lines
            continue
        if line[0] == "natoms":
            natoms = int(line[1])
        if line[0] == "boxsize":
            boxsize = float(line[1]), float(line[2]), float(line[3])
        if line[0] == "title":
            print "reading file:", " ".join(line[1:])
            title = " ".join(line[1:])
        if line[0] == "data":
            S = System(N=natoms, boxsize=boxsize)
            for i in range(natoms):
                line = fo.readline().split()
                S.q[i] = float(line[0]), float(line[1]), float(line[2])
    return S
                
        
        
if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "timing":
        disp = False
        S = System(N=600, trialMoveScale=.25)
        S.flags |= SVD_USE_PAIRLIST
        S.fillRandom()
        S.removeOverlaps()
        disp and S.makebox()
        for i in xrange(100):
            S.pairlistInit(2.5)
            disp and S.display()
            print i,
            S.trialMove(n=1000)
            S.pairlistCheck(1.5)
        print S.pairlist
        #S.trialMove_c(100000)
        #S.checkEnergyConsistency()
    elif len(sys.argv) > 1 and sys.argv[1] == "demo":
        S = System(N=500)
        S.fill()
        S.makebox()
        i = 0
        while True:
            S.trialMove()
            #S.checkEnergyConsistency()
            if i%100 == 0:
                S.display()
                print "\r%9d"%i, S.energy(), S.pressure()
                sys.stdout.flush()
            i += 1
    elif len(sys.argv) > 1 and sys.argv[1] == "exp1":
        opts, args = cliargs.get_cliargs()
        disp = False
        rho = float(args["rho"])
        T = float(args["T"])
        boxedge = 17.5
        S = System(N=int(rho * boxedge**3),
                   boxsize=(boxedge, boxedge, boxedge),
                   beta=1./T,
                   trialMoveScale=float(args["move"]))
        print S.N
        name = "eqled-rho_%s-temp_%s"%(args["rho"], args["T"])
        series = "2scale/"
        #name = "test"
        #S.fill(10, 10, 10, scale=1.72)
        S.fillRandom() 
        S.removeOverlaps()
        disp and S.makebox()
        disp and S.display()

        #logfile = file(series+name+".log", "w")
        logfile = file("plot.txt", "w")
        print >> logfile, "# i accRatio energy accRatio_last"
        def run(n):
            S.trialMove(n=n)
            disp and S.display()
            print "\r       \r", S.ntry,
            sys.stdout.flush()
            print >> logfile, \
                  S.ntry, round(S.acceptRatio(), 4), S.energy(), \
                  round(S.acceptRatio_last(), 4)
            logfile.flush()
            return S.ntry
        while run(500) != 100000:
            pass
        while run(10000) != 6000000:
            pass
        #S.writeUghfile(series+name+".ugh")
        print "accept ratio:", S.acceptRatio()
    elif len(sys.argv) > 1 and sys.argv[1] == "exp1_02":
        opts, args = cliargs.get_cliargs()
        disp = False
        rho = float(args["rho"])
        T = float(args["T"])
        boxedge = 17.5
        S = System(N=int(rho * boxedge**3),
                   boxsize=(boxedge, boxedge, boxedge),
                   beta=1./T,
                   trialMoveScale=float(args["move"]))
        name = "eqled-rho_%s-temp_%s"%(args["rho"], args["T"])
        series = "2scale/"
        #name = "test"
        S.readUghfile(series+name+".ugh")
        disp and S.makebox()
        disp and S.display()

        #insertType = 1
        logfile = file(series+name+".2.insert.log", "w")
        print >> logfile, "# i accRatio energy mu"
        def run2(n):
            S.trialMove(n=n)
            S.widomInsert(type=0, n=500)
            S.widomInsert(type=1, n=500)
            disp and S.display()
            print "\r        \r", S.ntry,
            sys.stdout.flush()
            print >> logfile, \
                  S.ntry, round(S.acceptRatio(), 4), S.energy(), \
                  S.mu(type=0), \
                  S.mu(type=1)
            logfile.flush()
            return S.ntry
        while run2(10000) != 2000000:
            pass
        print "std dev of insert dict:", numpy.std(S.mu_dict[0]), \
              numpy.std(S.mu_dict[1])
        print "RESULT::%s:: %s %s"%(1, args, S.mu(type=1))
        print "RESULT2:::%s:::"%{"args": args,
                                 "mu_0": S.mu(type=0),
                                 "mu_1": S.mu(type=1),
                                 "mu_0_pstddev": numpy.std(S.mu_dict[0]),
                                 "mu_1_pstddev": numpy.std(S.mu_dict[1])
                                 }
        
        print "accept ratio:", S.acceptRatio()
        sys.exit(0)
    elif len(sys.argv) > 1 and sys.argv[1] == "gen_config":
        N = int(sys.argv[2])
        S = System(N=N)
        S.fill()
        outfile = "ugh/equilibrated-%s.ugh"%N
        logfile = outfile + ".log"
        log = file(logfile, "w")
        for i in xrange(3000000):
            S.trialMove(verbose=False)
            if i%10000 == 0 or (i<100000 and i%500 == 0):
                line = "%9s %9.9s %9.9s"%(S.ntry, S.energy() , S.pressure())
                print line
                print >>log, line
                log.flush()
        S.writeUghfile("ugh/equilbrated-%s.ugh"%N)
    else:
        skip = 1
        S = System(N=200, beta=1/2.0,
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
        logfile = file("logfile.txt", "w")
        #logfile=sys.stdout
        print >> logfile, "# i, T, rho, E, P, P_avg"
        while True:
            if pairlist: S.pairlistInit(3.25)
            i += skip
            if i == 150000:
                S.resetStatistics()
            #S.trialMove(verbose=False, n=skip)

            ke = c_mdStep(S.SD_p, skip, 0)
            if pairlist: S.pairlistCheck(2.75)
            #S.widomInsert()
            S.display()
            print >> logfile, i, S.T, S.density, S.energy(), \
                  S.pressure(add=True), S.pressureAverage(), ke
            #print S.pairlist
            #S.pairlistCheck(strict=False)
            logfile.flush()
            #time.sleep(.5)
            print

    #time.sleep(10)
