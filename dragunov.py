# Richard Darst, 2007
# The University of Texas at Austin
from __future__ import division

# Hi, welcome to Dragunov.  I'll be your guide to the code today.
# 
# There are two main components, this python file and the C file.  The
# python file is the only one you'll need to use in your code.  The
# python imports the C file and calls its functions.  If you
# understand everything from the python file, the C file is just small
# additions which you can see later.
#
# Most of the code was originally written in python, and then only
# later transferred to C.  So, there is left over methods in the
# python module which aren't used very often.  There are switches
# which select between the python and the C implementations of
# different things (like making a trial move).  For example, you can see:
#
#    def trialMove_py(self, verbose=False):
#        [...]
#    def trialMove_c(self, n=1, verbose=False):
#        [...]
#    trialMove = trialMove_c
#
# The meaning of this should be clear.  Most C/python functions should
# be similar, but there are still rough edges.
#
# Let's continue the tour.  I'm only going to point out interesting
# things.
#
# COMPILING
#
# Everything you need should be in the build.sh script, just do the
# command "sh build .sh".  If you look in the file, there is a way to
# compile it with a different compiler if you want to try to make it
# go faster.  It makes the module dragunov_c.so which is imported by
# here.
#
# Other dependencies:
#   ctypes
#   python-visual (for fancy displays!)

import math
import random
random.seed(12345)
import os.path
import sys
import time
import ctypes
import numpy
numpy.random.seed(123456)
try:
    import visual
except ImportError:
    # display won't work, but no error raised.  The visual functions
    # will detect that the module wasn't imported, and silently
    # return.  This makes it easier to run jobs on clusters with no
    # modifications.
    # This could be exploited by doing `dragunov.visual = None` to
    # explicitely deactivate it.
    visual = None
try:
    from rkddp.interact import interact
except ImportError:
    pass

# prevent an "imported module not used" error when checked with
# pychecker
time.sleep(0)

# These flags are also defined in the C file, and allow different
# switches to be tuned via the "flags" int paramater which is passed
# to C.

SVD_ENERGYI_PARTIAL       =   1
SVD_VERBOSE_1             =   2
SVD_VERBOSE_2             =   4
SVD_PAIRLIST_INCREMENTAL  =   8
SVD_USE_PAIRLIST          =  16

# This is the primary data structure.  It uses the module `ctypes`.
# Ctypes is black magic, and I'm not going to explain it-- its
# documentation is pretty good.
#
# The thing below is a SimData object.  It represents a C structure
# (`struct SimData` in the C code).  It is passed to all of the C
# functions, and contains pointers to things like positions (q),
# boxsize, and so on.
#
# The structure doesn't really contain the data, but pointers to the
# data.  I make all arrays with numpy and give the SimData struct
# pointers to them.  Since it's all pointers, both the python and C
# always see the same data and are synced.  But, you have to always
# modify arrays in-place.

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
                #("ei", ctypes.c_void_p),
                ("pairlist", ctypes.c_void_p),
                ("trialMoveScale", ctypes.c_double),
                ("trialMoveIsobaricScale", ctypes.c_double),
                ("pairlist_minDistance", ctypes.c_double),
                ("prob_PMove", ctypes.c_double),
                ("isobaricPressure", ctypes.c_double),
                ("ntry_shift",            ctypes.c_int),
                ("naccept_shift",         ctypes.c_int),
                ("ntry_shift_last",       ctypes.c_int),
                ("naccept_shift_last",    ctypes.c_int),
                ("ntry_isobaric",         ctypes.c_int),
                ("naccept_isobaric",      ctypes.c_int),
                ("ntry_isobaric_last",    ctypes.c_int),
                ("naccept_isobaric_last", ctypes.c_int),
                ]
SimData_p = ctypes.POINTER(SimData)

# Now, we load the C module using ctypes.  This is magic.  All the
# code below sets up python interfaces to the C functions.  You can
# make enough sense out of it by just glancing at it, to really
# understand look at the ctypes docs.

CLibraryLookup = { }
ForceFields = {
    "hardsphere": 1,
    "lennardjones": 2,
    "harmonic": 3,
    "2scale-s2s": 10,
    "3scale-s3s": 11,
    "3scale-s2s": 12,
    "2scale-jagla98": 13,
    }
def getCLibrary(forceField):
    """Return the C library corresponding to this force field.

    
    """
    if not ForceFields.has_key(forceField):
        print "unknown forcefield:", forceField
    if not CLibraryLookup.has_key(forceField):
        # add it to the dict
        num = ForceFields[forceField]
        dragunov_c = loadCLibrary(filename="dragunov_%02d_c"%num)
        CLibraryLookup[forceField] = dragunov_c
    return CLibraryLookup[forceField]
def loadCLibrary(filename,
                 path=os.path.dirname(__file__)):
        
    dragunov_c = numpy.ctypeslib.load_library(filename,
                                              os.path.dirname(__file__))
    
    c_eij = dragunov_c.eij
    dragunov_c.eij.restype = ctypes.c_double
    dragunov_c.eij.argtypes = (ctypes.c_int,        # atomtype of i
                               ctypes.c_int,        # atomtype of j
                               ctypes.c_double, )   # distance

    c_fij = dragunov_c.fij
    dragunov_c.fij.restype = ctypes.c_double
    dragunov_c.fij.argtypes = (ctypes.c_int,        # atomtype of i
                               ctypes.c_int,        # atomtype of j
                               ctypes.c_double, )   # distance

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

    c_forcedotr_total = dragunov_c.forcedotr_total
    c_forcedotr_total.restype = ctypes.c_double
    c_forcedotr_total.argtypes = (SimData_p,       # SimData_p
                                  ctypes.c_int)    # flags

    c_nCloserThan = dragunov_c.nCloserThan
    c_nCloserThan.restype = ctypes.c_int
    c_nCloserThan.argtypes = (SimData_p,       # SimData_p
                              ctypes.c_double, # maxdist
                              ctypes.c_int)    # flags

    c_calcForce = dragunov_c.calcForce
    c_calcForce.restype = ctypes.c_double    # always returns zero
    c_calcForce.argtypes = (SimData_p,       # qdata_p
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
    return dragunov_c

class System(object):
    """Simulation Object.

    The main sim-object.  Usually one is initialized per program, but
    there is nothing preventing having multiple of them simultaneously
    (replica exchange, anyone?)
    """
    # These thing define "properties" of the class.  They are like
    # attribute, but are dynamically calculated like methods.
    def _volume_get(self):
        return numpy.product(self.boxsize)
    volume = property(fget=_volume_get)

    def _density_get(self):
        return self.N / numpy.product(self.boxsize)
    density = property(fget=_density_get)

    def _temperature_get(self):
        return 1 / self.beta
    T = property(fget=_temperature_get)
    def _qdot_get(self):
        # warning: this doesn't yield the KE centered at this point,
        # that should be (q(t+1) + q(t-1))/(2duut), not
        # (q(t)-q(t-1))/dt (see the integrate method)
        return (self.q - self.qold) / self.dt
    qdot = property(fget=_qdot_get)

    def __init__(self, N, forceField,
                 beta=1., Nmax=None, boxsize=(10,10,10),
                 trialMoveScale=.25, trialMoveIsobaricScale=.10,
                 isobaricPressure=None,
                 dt=.01):
        """All initilization

        Initilize:
        - numpy arrays
        - C pointers to numpy arrays
        - store other constants of the simulation, like temperature
        - lists for storing time averages.

        Most of these constants are stored in two places: on the
        python object (self) and in the SimData structure (self.SD)
        for access in C.  Be aware of keeping them both updated!
        """
        self.C = getCLibrary(forceField)
        if Nmax == None:
            Nmax = N+50
        SD =      self.SD = SimData()
        self.SD_p = ctypes.pointer(SD)
        SD.N    =   self.N    =   N
        SD.Nmax =   self.Nmax =   Nmax
        SD.dt =     self.dt   =   dt
        SD.beta =   self.beta =   beta
        SD.prob_PMove = self.prob_PMove = 0
        self.mcsteps = 0
        self.mctime = 0

        SD.trialMoveScale = self.trialMoveScale = trialMoveScale
        SD.trialMoveIsobaricScale = self.trialMoveIsobaricScale = \
                                    trialMoveIsobaricScale
        self.mu_dict = { }
        self._pressureList = [ ]
        self._volumeList   = [ ]
        if isobaricPressure:
            SD.isobaricPressure = self.isobaricPressure = isobaricPressure
            # default to 1 V move for every N+1 steps
            self.setMoveProb(shift=self.N, isobaric=1)
            print "info: setting to isobaric ensemble (NPT)"
        else:
            self.setMoveProb(shift=self.N, isobaric=0)


        self.naccept = 0
        self.ntry = 0
        self.SD.ntry_shift            = 0
        self.SD.naccept_shift         = 0
        self.SD.ntry_shift_last       = 0
        self.SD.naccept_shift_last    = 0
        self.SD.ntry_isobaric         = 0
        self.SD.naccept_isobaric      = 0
        self.SD.ntry_isobaric_last    = 0
        self.SD.naccept_isobaric_last = 0


        # numpy.float is ctype double
        self.q = numpy.zeros(shape=(self.Nmax, 3), dtype=numpy.float)
        SD.q   = self.q.ctypes.data
        self.qold = numpy.zeros(shape=(self.Nmax, 3), dtype=numpy.float)
        SD.qold   = self.qold.ctypes.data
        self.force = numpy.zeros(shape=(self.Nmax, 3), dtype=numpy.float)
        SD.force   = self.force.ctypes.data
        #self.ei = numpy.zeros(shape=(self.Nmax, ), dtype=numpy.float)
        #SD.ei   = self.ei.ctypes.data
        self.pairlist = numpy.zeros(shape=(self.Nmax, 152), dtype=numpy.int_)
        SD.pairlist   = self.pairlist.ctypes.data
        self.boxsize = numpy.asarray(boxsize, dtype=numpy.float)
        SD.boxsize   = self.boxsize.ctypes.data
        self.atomtypes = numpy.zeros(shape=(self.Nmax), dtype=numpy.int_)
        SD.atomtypes   = self.atomtypes.ctypes.data
        self.atomtypes[0:self.N] = 0   # default to zero

        self.flags = 0

    def resetStatistics(self):
        """Reset averages

        This method resets the averages for pressure, volume, chemical
        potential, etc.
        """
        self._pressureList = [ ]
        self._volumeList   = [ ]
        self.mu_dict = { }
        self.SD.ntry_shift_last       = self.SD.ntry_shift
        self.SD.naccept_shift_last    = self.SD.naccept_shift
        self.SD.ntry_isobaric_last    = self.SD.ntry_isobaric
        self.SD.naccept_isobaric_last = self.SD.naccept_isobaric
    def resetTime(self):
        """Set MC time and MC steps to zero.
        """
        self.mcsteps = 0
        self.mctime = 0
    def acceptRatioShift(self):
        """Return MC acceptance ratio"""
        if self.SD.ntry_shift == 0:
            return float("nan")
        return self.SD.naccept_shift / self.SD.ntry_shift
    def acceptRatioShiftLast(self):
        """Return MC acceptance ratio"""
        if self.SD.ntry_shift    - self.SD.ntry_shift_last == 0:
            return float("nan")
        return ((self.SD.naccept_shift - self.SD.naccept_shift_last) /
                (self.SD.ntry_shift    - self.SD.ntry_shift_last))
    def acceptRatioIsobaric(self):
        """Return MC acceptance ratio"""
        if self.SD.ntry_isobaric == 0:
            return float("nan")
        return self.SD.naccept_isobaric / self.SD.ntry_isobaric
    def acceptRatioIsobaricLast(self):
        """Return MC acceptance ratio"""
        if self.SD.ntry_isobaric    - self.SD.ntry_isobaric_last == 0:
            return float("nan")
        return ((self.SD.naccept_isobaric - self.SD.naccept_isobaric_last) /
                (self.SD.ntry_isobaric    - self.SD.ntry_isobaric_last))
    def loginfo(self, fileobject=None):
        """Write a header to the logfile.

        Indicates order of fields logged."""
        if fileobject is None:
            fileobject = sys.stdout
        print >> fileobject, \
              "# mcsteps mctime T N density E P Pavg V Vavg AR_shift AR_shift_last AR_isobaric AR_isobaric_last"
    
    def log(self, fileobject=None, pressure=True,):
        """Write a standard line to the logfile

        If `fileobject` is Nono, write to stdout.  The order of fields
        is:

        # 1       2      3 4 5       6 7 8    9 10
        # mcsteps mctime T N density E P Pavg V Vavg
        # 11       12            13          14
        # AR_shift AR_shift_last AR_isobaric AR_isobaric_last
        """
        if fileobject is None:
            fileobject = sys.stdout
        print >> fileobject, self.mcsteps, self.mctime, self.T, self.N, \
              self.density, self.energy(), \
              self.pressure(), self.pressureAverage(), \
              self.volume, self.volumeAverage(), \
              self.acceptRatioShift(), self.acceptRatioShiftLast(), \
              self.acceptRatioIsobaric(), self.acceptRatioIsobaricLast()

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
        #for i in range(self.N):
        #    self.ei[i] = self.energy_i(i)
    def fillRandom(self):
        """Fill the box with all-random particle locations.

        There is no checking for overlaps, so you must call a method
        to remove overlaps/minimize it."""
        for i in range(self.N):
            newpos = numpy.random.uniform(size=(3,)) * self.boxsize
            self.q[i] = newpos
    def zeroVelocity(self):
        """Zero all velocities.

        The implementation is by setting the qold to q, as per the
        verlet integrator used.  This only has any significance if you
        are using the MD integrator, which is a corner case (not what
        dragunov is designed for, really.)
        """
        self.qold[:] = self.q[:]
    def dij(self, i, j):
        d = self.q[i] - self.q[j]
        d = numpy.abs(d - (numpy.floor(d/self.boxsize + .5)) * self.boxsize)
        d = numpy.sqrt((d*d).sum())  # distances from i to j
        return d

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
        """Return total energy by evaluating energy for all atoms

        Uses C function c_energy"""
        #E = 0
        #for i in range(self.N):
        #    E += c_energy_i(self.SD_p, i, SVD_ENERGYI_PARTIAL)
        #    #E+= self.energy_i(i, partial=True)
        #return E
        return self.C.energy(self.SD_p, self.flags)

    #def energy_fromi(self):
    #    """Return total energy by using sum of cached values for each atom"""
    #    raise Exception("this doesn't work")
    #    return .5 * sum(self.ei[0:self.N])
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
        E = self.C.energy_i(self.SD_p, i, flags)
        return E
    energy_i = energy_i_c # select which energ eval function to use

    def calcForce(self):
        """Run the c-force calculation routine.

        This must be done before self.force will be updated.
        """
        self.C.calcForce(self.SD_p, self.flags)
    def pressure_c(self, add=True):
        fdotr = self.C.forcedotr_total(self.SD_p, self.flags)
        volume = self.volume
        dimensions = 3
                 # v-- should be "+" for LJ correct results
                 # v-- this also makes pressure increase with
                 #     increasing density for S2S- 2s units
        pressure = + fdotr / (dimensions * volume)

        if False:
            dr = .005
            R = 1.
            n = self.C.nCloserThan(self.SD_p, R+dr, self.flags)
            print n
            pressure += R * n / (dimensions * volume * dr * self.beta)
        if add:
            pressure += self.density/self.beta
        self._pressureList.append(pressure)
        self._volumeList.append(volume)
        return pressure
    pressure = pressure_c
    def pressureAverage(self):
        return sum(self._pressureList)/len(self._pressureList)
    def volumeAverage(self):
        return sum(self._volumeList)/len(self._volumeList)

    def trialMove_py(self, verbose=False):
        """Try a move using metropolis criteria"""
        raise Exception, "Function fallen out of sync with the C version."




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
            #self.ei[i] = Enew
            self.naccept += 1

        else:
            self.q[i] = qi_old


            
        if verbose:
            print "%6d %d, %5.3f->%5.3f accF: %.4f\n"%(
                self.N, int(accept), Eold, Enew, self.naccept/self.ntry)

    def trialMove_c(self, n=1, verbose=False):
        """Do n move cycles.

        A move cycle is defined as running the number of steps defined
        with the setMoveProb method.  By default, a cycle will run
        Natoms steps.  Thus, each trial move runs approximately the
        same amount of real sampling regardless of the number of atoms
        in the system."""
        # We have to do an explicit loop, or else we won't be able to
        # do the pairlist regeneration at each step.
        for i in xrange(n):
            if self.flags & SVD_USE_PAIRLIST:
                if self._pairlist_warn and self.mctime != 0:
                    self.pairlistCheck()
                self.pairlistInit()
            nsteps = self._movesPerCycle
            naccept = self.C.trialMove(self.SD_p, nsteps, self.flags)
            self.mcsteps += nsteps
            self.mctime += 1
        #if verbose:
        #    print "%10s moves"%self.ntry, round(self.naccept/self.ntry, 4)
        
    trialMove = trialMove_c
    def widomInsert(self, type=0, n=None):
        """Save a data point for widom insertion.

        Insert one atom of `type` at a random location evenly
        distributed throughout the box.  Record the increase of energy
        caused by this atom.  Remove the inserted atom.  Repeat `n`
        times for better sampling.

        Since insertion is a relatively cheap process compared to
        doing 1000 trial moves, it's probably best to do trial moves
        in groups and then multiple insertions at once.
        """
        # n is the number of test moves to run
        #if n is not None:
        #    for _ in xrange(n):
        #        self.widomInsert(type=type)
        #    return
        if n == None:
            n = 1
        for _ in xrange(n):
            newpos = numpy.random.uniform(size=(3,)) * self.boxsize
            newi = self.addAtom(newpos, type=type)
            Eadded = self.energy_i(newi)
            #print Eadded, self.beta
            self.mu_dict.setdefault(type, []).append(math.exp(-Eadded*self.beta))
            self.delAtom(newi)
        
    def widomInsertResults(self, type=0):
        """Return the chemical potential of a species.

        This is a separate function, because it must not only average
        the stored values, but take a log and multiply by constants.

        This only returns the excess chemical potential, not the ideal
        gas contribution (see frenkel and smit).
        """
        volume = numpy.product(self.boxsize)
        avg = sum(self.mu_dict[type]) / len(self.mu_dict[type])
        #print self.mu_dict[type]
        if avg == 0:
            return 0  # this is technically incorrect, but we need to return
                      # something that won't make parse errors

        # must be adjusted to debroglie wavelength, if necessary.
        #mu_ig = -math.log(volume/( (self.N+1)* 1.**3) )/self.beta
        mu_excess = -math.log(avg)/self.beta
        return mu_excess #+ mu_ig
    def widomInsertCorrection1(self):
        """Unused

        Used to correct for the ideal gas contribution to the chemical
        potential.
        """
        # _add_ this to the reported mu tocorrect the mu.
        #constant = 3./self.beta
        #eturn math.log(self.volume/(self.N+1))/self.beta - constant
        pass
    mu = widomInsertResults

    def setMoveProb(self, shift, isobaric=0):
        """Set trial move probabilities and timescale of simulation.

        If we are doing an ensemble like NTP, we need to do both
        particle moves AND moves in volume-space.  This method adjusts
        the relative probabilities of the different types of moves.
        For example, we can have it do (on average) one volume move
        for every N particle moves.  (this is the default, as set in
        the __init__() method).

        The arguments are `shift` and `pressure`, which are for the
        relative probs for each of the different types of moves.  The
        probs are normalized before setting, so you can do::

          setMoveProb(shift=self.N, pressure=1)

        so that a pressure move is done on average once for every N
        shifts.

        The second part of this refers to setting the timescale of
        simulation.  Each time trialMove is called, it'll do a number
        of moves equal to the sum of the numbers passed here.
        """
        sum_ = shift + isobaric
        self._movesPerCycle = sum_
        self.SD.prob_PMove = self.prob_PMove = isobaric / sum_
        
    def trialMove_isobaric_py(self):
        """Volume-adjusting move for the isobaric ensemble.

        Do a random walk in `ln(V)`, maintaining a constant pressure.

        This function has been replaced by the C one, the python one
        isn't used anymore.
        """
        numpy.product(self.boxsize)
        pressure = self.isobaricPressure
        lnVScale = self.trialMoveIsobaricScale

        Vold = self.volume
        Eold = self.energy()

        #Vnew = exp(ln(V + (random.random()-.5) * scale))
        #lengthscale = (Vnew / Vold)**(1./3.) # this may not be right
        linearScale = math.exp(((random.random()-.5) * lnVScale) / 3.)
        volumeScale = linearScale * linearScale * linearScale

        self.q *= linearScale
        self.boxsize *= linearScale

        Enew = self.energy()
        Vnew = Vold * volumeScale


        x = - self.beta * (Enew - Eold + pressure*(Vnew-Vold) - \
                           ((self.N+1)/self.beta)*math.log(volumeScale))
        print x

        if x > 0:
            accept = True
        else:
            x = math.exp(x)
            print x
            ran = random.random()
            if ran < x:  accept = True
            else:        accept = False

        if accept:
            print "+++ %.3f pressure move  +++"%linearScale
            pass
        else:
            print "--- %.3f pressure move  ---"%linearScale
            self.q /= linearScale
            self.boxsize /= linearScale
        

    def addAtom(self, pos, type):
        """Add an arbitrary atom to the system.  Return new index"""
        if self.N == self.Nmax:
            raise Exception("can't fit more atoms here")
        self.SD.N    =   self.N    =   self.N+1
        i = self.N-1
        self.q[i] = pos
        self.atomtypes[i] = type
        #self.ei[i] = self.energy_i(i)
        return i
    def delAtom(self, i):
        """Delete an arbitry atom from the system."""
        N = self.N
        if i != N-1:
            raise Exception("so far we can only remove the last atom")
        self.SD.N    =   self.N    =   self.N-1
        self.q[i] = 0.
        self.atomtypes[i] = 0.
        #self.ei[i] = 0.

    def checkContiguous(self):
        """Check if the numpy array are contiguous, if not raise an error."""
        # should be modified to check all arrays
        if not self.q.flags["C_CONTIGUOUS"]:
            print "*** not C contiguous!"
            raise Exception
    def checkIntersected(self):
        # self-test function (not used)
        partial = 1
        energy_i = self.C.energy_i
        for i in xrange(self.N):
            energy_i(self.SD_p, i, partial|2)
    def checkEnergyConsistency(self):
        """Check of energy self-consistency.

        Do the functinos which calculate total energy of all atoms, vs
        those that do it for only one atom, give same answer for the
        energy?

        Used as a check, not in the main simulation.
        """
        E_fromall = self.energy_fromall()
        E_fromi = self.energy_fromi()
        deltaE = E_fromall - E_fromi
        print "E_fromall:",E_fromall, "E_fromi",E_fromi, \
              "E difference:", repr(deltaE)
        assert abs(E_fromall - E_fromi) < .001
    def forceFromNDeriv(self):
        """Numerically differentiate the energy, return force matrix

        This provides a check for the force calculation.  Not used in
        simulation.
        """
        # Use this algorithm: - (E(x+dx) - E(x-dx)) / 2*dx
        dx = .01
        q = self.q
        force = self.force.copy()
        force[:] = 0
        for n in range(self.N):
            for i in range(3):
                qorig = q[n,i]
                q[n,i] -= dx
                Elow = self.energy_i(n)
                q[n,i] = qorig + dx
                Ehigh = self.energy_i(n)
                q[n,i] = qorig
                force[n,i] = - (Ehigh-Elow) / (2 * dx)
        return force
    # There are various functions I've made for calculating pairlists,
    # but they haven't been tested fully.  (I have tested them, and
    # they seem to work, but I need to double check them.)
    def pairlistInit(self, cutoff=None):
        """Fill the pairlist with distance information.

        See `pairlistEnable`.
        """
        if cutoff is None:
            cutoff = self._pairlist_cutoff
        flags = self.flags
        return self.C.pairlistInit(self.SD_p, cutoff, flags)
        
    def pairlistCheck(self, warn=None, strict=False):
        """Check pairlists for too-close violations

        See `pairlistEnable`.
        """
        if warn is None:
            warn = self._pairlist_warn
            strict = self._pairlist_strict
        nviolations = self.C.pairlistCheck(self.SD_p, warn, self.flags)
        if nviolations:
            print "nviolatios:", nviolations, "at step", self.mctime
        if strict:
            assert nviolations == 0
        #print self.SD.pairlist_minDistance
    def pairlistEnable(self, cutoff, warn=None, strict=False):
        """Enable and store data about pairlists.

        This function is your one stop for enabling everything related
        to pairlists.  This function sets the SVD_USE_PAIRLIST flag
        for the system object, which makes pairlist get regerated
        automatically every cycle when trialMove is called (once per
        cycle, not once per trialMove function call).

        Use the `cutoff` parameter to set the distance at which atoms
        are recorded.

        If the parameter `warn` is set, at each cycle the `
        pairlistCheck` method will be called, to see if there are any
        atoms at a distance less than `warn`.  If so, a message will
        be printed.  Furthermore, if `strict` is set to be true, then
        an AssertionError will be raised, abruptly halting all
        computations.
        """
        self.flags |= SVD_USE_PAIRLIST
        self._pairlist_cutoff = cutoff
        self._pairlist_warn = warn
        self._pairlist_strict = strict
        
        
    def removeOverlaps(self, Emax=float('inf')):
        """Remove overlaps of atoms (make energy non-infinite)

        Iterates through all atoms, and for each one with infinite
        energy (or energy above `Emax`), give it random displacements
        until it has finite energy.
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
            #self.ei[i] = Eold
        #print "done removing overlaps"

    def qWrapped(self):
        """Return coordinates wrapped to fit in boxsize

        Used for display."""
        # wrapped means in interval [0, boxsize)
        q = self.q
        return numpy.mod(q, self.boxsize)
    def display(self):
        """Update visual display of atoms.

        Requires python-visual to be installed.  Call this function
        again to update positions"""
        if visual is None:
            return
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
        if visual is None:
            print "visual module was not imported, visual siletly deactivated."
            return
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
    def hash(self):
        """Return a checksum of system state

        Note that Python defines a hash very specifically--
        specifically, a hash must stay the same for anything with the
        same data, and a hash must not change (so that it can be used
        for keys for things like dictionaries.)  This hash
        implementation violates both of these:
        - There is data that can change in this object and not be
          reflected in the hash
        - The object is mutable, and the hash can change.
        As such, this is defined as `hash`, not the magic method
        `__hash__`.

        The data this hash is computed by is: q, qold, boxsize,
        atomtypes, N
        """
        x = ( hash(tuple( tuple(i) for i in self.q)),
              hash(tuple( tuple(i) for i in self.qold)),
              hash(tuple( i for i in self.boxsize)),
              hash(tuple( i for i in self.atomtypes)),
              self.N
              )
        return hash(x)
            

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
    import cliargs
    if len(sys.argv) > 1 and sys.argv[1] == "timing":
        # test timings of the simulation.
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
        # show a nice demo to impress your boss.
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
        # something to do an experiment
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
        # another experiment code
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
        # generate equilibrated state points.
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
    elif len(sys.argv) > 1 and sys.argv[1] == "test_ff":
        # I guess this tests the force fields.
        #N = 200
        N = 2
        S = System(N=N, beta=1/2.0,
                   boxsize=(10,10,10), trialMoveScale=.25,
                   dt=.001)
        #S.fillRandom()
        S.fill()
        S.zeroVelocity()
        #S.removeOverlaps(50.)
        S.makebox() ; S.display()
        skip = 1
        i = 0
        while True:
            i += skip
            if i == 150000:
                S.resetStatistics()
            #S.trialMove(verbose=False, n=skip)
            #S.calcForce()
            print S.force[:N]
            print S.forceFromNDeriv()[:N]
            print
            ke = S.C.mdStep(S.SD_p, skip, 0)
            S.display()
            #print >> logfile, i, S.T, S.density, S.energy(), \
            #      S.pressure(add=True), S.pressureAverage(), ke
    else:
        # random playing around with it.
        skip = 1
        N = 200
        #N = 2
        S = System(N=N, beta=1/2.0,
                   boxsize=(10,10,10), trialMoveScale=.25,
                   dt=.0001)
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

            ke = S.C.mdStep(S.SD_p, skip, 0)
            if pairlist: S.pairlistCheck(2.75)
            #S.widomInsert()
            S.display()
            print >> logfile, i, S.T, S.density, S.energy(), \
                  S.pressure(add=True), S.pressureAverage(), ke
            #print S.pairlist
            #S.pairlistCheck(strict=False)
            logfile.flush()
            #time.sleep(.5)
            #print

    #time.sleep(10)
