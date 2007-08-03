runs = [
    {"rho":1.00, "T":.03, "move":.15 },
    {"rho":1.00, "T":.04, "move":.15 },
    {"rho":1.00, "T":.05, "move":.20 },
    {"rho":1.00, "T":.06, "move":.20 },
    {"rho":1.00, "T":.07, "move":.25 },

    {"rho":1.21, "T":.03, "move":.15 },
    {"rho":1.21, "T":.04, "move":.15 },
    {"rho":1.21, "T":.05, "move":.20 },
    {"rho":1.21, "T":.06, "move":.20 },
    {"rho":1.21, "T":.07, "move":.20 },

    {"rho":1.66, "T":.03, "move":.12 },
    {"rho":1.66, "T":.04, "move":.12 },
    {"rho":1.66, "T":.05, "move":.15 },
    {"rho":1.66, "T":.06, "move":.15 },
    {"rho":1.66, "T":.07, "move":.15 },

    {"rho":2.19, "T":.03, "move":.10 },
    {"rho":2.19, "T":.04, "move":.10 },
    {"rho":2.19, "T":.05, "move":.10 },
    {"rho":2.19, "T":.06, "move":.10 },
    {"rho":2.19, "T":.07, "move":.10 },
    ]

#runs = [
#    # for 3-scale-scale in S-2S units
#    {"rho":1.00, "T":.03, "move":.15 },
#    {"rho":1.00, "T":.06, "move":.20 },
#    {"rho":1.00, "T":.09, "move":.30 },
#    {"rho":1.00, "T":.12, "move":.35 },
#    {"rho":1.00, "T":.15, "move":.35 },
#    
#    {"rho":1.50, "T":.03, "move":.10 }, 
#    {"rho":1.50, "T":.06, "move":.18 },
#    {"rho":1.50, "T":.09, "move":.25 },
#    {"rho":1.50, "T":.12, "move":.30 },
#    {"rho":1.50, "T":.15, "move":.30 },
#    
#    {"rho":2.00, "T":.03, "move":.10 },
#    {"rho":2.00, "T":.06, "move":.15 },
#    {"rho":2.00, "T":.09, "move":.20 },
#    {"rho":2.00, "T":.12, "move":.20 },
#    {"rho":2.00, "T":.15, "move":.25 },
#
#    {"rho":2.00, "T":0.3, "move":.25 },
#    {"rho":2.00, "T":0.7, "move":.25 },
#    {"rho":2.00, "T":1.0, "move":.25 },
#    {"rho":2.00, "T":1.5, "move":.25 },
#    {"rho":2.00, "T":3.0, "move":.25 },
#    {"rho":2.00, "T":6.0, "move":.25 },
#    {"rho":2.00, "T":12., "move":.25 },
#]

 
def makeArray():
    a = { }
    for d in runs:
        a.setdefault(d["rho"], {}).setdefault(d["T"], {})
    return a

import math
import os
import sys
import dragunov
import numpy


class Experiment:
    disp = False
    def __init__(self, args):
        """args dict must have keys 'rho', 'T', and 'move'
        """
        self.args = args
        self.rho = float(args["rho"])
        self.T = float(args["T"])
        self.move = float(args["move"])
        self.name = "exp-rho_%s-temp_%s"%(self.rho, self.T)
        self.series = "2scaleB/"
        self.boxedge = 10

    def makeSys(self):
        S = dragunov.System(N=int(self.rho * self.boxedge**3),
                            boxsize=(self.boxedge, self.boxedge, self.boxedge),
                            beta=1./self.T,
                            trialMoveScale=self.move)
        return S

    def equilibrate(self, steps):
        S = self.makeSys()
        print S.N
        #name = "test"
        #S.fill(10, 10, 10, scale=1.72)
        S.fillRandom()
        S.removeOverlaps()
        self.disp and S.makebox()
        self.disp and S.display()

        logfile = file(self.series+self.name+".eq.0.log", "w")
        #logfile = file("plot.txt", "w")
        print >> logfile, "# i accRatio energy accRatio_last"
        def run(n):
            S.trialMove(n=n)
            self.disp and S.display()
            print "\r       \r", S.ntry,
            sys.stdout.flush()
            print >> logfile, \
                  S.ntry, round(S.acceptRatio(), 4), S.energy(), \
                  round(S.acceptRatio_last(), 4)
            logfile.flush()
            return S.ntry
        while run(500) != 100000:
            pass
        while run(10000) != steps:
            pass
        print
        S.writeUghfile(self.series+self.name+".eq.0.ugh")
        print "accept ratio:", S.acceptRatio()

    def sample(self, steps):
        S = self.makeSys()
        S.readUghfile(self.series+self.name+".eq.0.ugh")
        self.disp and S.makebox()
        self.disp and S.display()

        #insertType = 1
        logfile = file(self.series+self.name+".insert.0.log", "w")
        print >> logfile, "# i accRatio energy accRatioLast mu0 mu1 pressure"
        def run2(n):
            S.trialMove(n=n)
            S.widomInsert(type=0, n=500)
            S.widomInsert(type=1, n=500)
            self.disp and S.display()
            print "\r        \r", S.ntry,
            sys.stdout.flush()
            print >> logfile, \
                  S.ntry, round(S.acceptRatio(), 4), \
                  S.energy(), \
                  round(S.acceptRatio_last(), 4), \
                  S.mu(type=0), \
                  S.mu(type=1), \
                  S.pressure()
            logfile.flush()
            return S.ntry
        while run2(10000) != steps:
            pass
        print
        #print "std dev of insert dict:", numpy.std(S.mu_dict[0]), \
        #      numpy.std(S.mu_dict[1])
        #print "RESULT::%s:: %s %s"%(1, args, S.mu(type=1))
        resultLine= "RESULT2:::%s:::"%{"T": self.T,
                                 "rho": self.rho,
                                 "move": self.move,
                                 "args": self.args,
                                 "mu_0": S.mu(type=0),
                                 "mu_1": S.mu(type=1),
                                 "mu_0_pstddev": numpy.std(S.mu_dict[0]),
                                 "mu_1_pstddev": numpy.std(S.mu_dict[1]),
                                 "pressureAverage": S.pressureAverage(),
                                 "ntry": S.ntry,
                                 "naccept": S.naccept,
                                 "Paccept": S.naccept/S.ntry,
                                 }
        print resultLine
        print >> logfile, resultLine
        print "accept ratio:", S.acceptRatio()
        print >> logfile, resultLine

if __name__ == "__main__":
    if sys.argv[1] == "1":
        print "building:"
        os.system("sh ./build.sh")
        for R in runs:
            print "running:", R
            E = Experiment(R)
            print E
            E.equilibrate(3000000)
        #    r = os.system("time python main.py exp1 "
        #              "--T=%(T)s --rho=%(rho)s --move=%(move)s "%
        #              R)
        #    print r
        #    if r: sys.exit(r)
        #    print
        #    print
        

    if sys.argv[1] == "2":
        print "building:"
        os.system("sh ./build.sh")
        for R in runs:
            print R
            print "running:", R
            E = Experiment(R)
            print E
            E.sample(1000000)
            #r = os.system("time python main.py exp1_02 "
            #          "--T=%(T)s --rho=%(rho)s --move=%(move)s "%
            #          R)
            #if r>>8:
            #    print "exiting, r =", r>>8
            #    sys.exit(r>>8)
            #print
            #print

    if sys.argv[1] == "3":
        import re
        from rpy import r
        from rkddp.interact import interact
        from rkddp import myr

        #logfile = file("2scale-log2.txt").read()
        logfile = file(sys.argv[2]).read()
        array = makeArray()
        Legend = myr.Legend()
        #species = "K_1"
        species = "mu_1"
        plottitle = "2-scale insertion, species="+species,
        if False:
            pass
            ##  # Get results from the log file
            ##  result_rec = re.compile(r'RESULT::([0-9]+):: (\{.*?\}) ([^ ]*?)\n')
            ##  start = 0
            ##  #for R in runs:
            ##  #logfile = file("exp1_02b.log").read()
            ##  f = file("plot02b.txt", "w")
            ##  while True:
            ##      m = result_rec.search(logfile, start)
            ##      if not m:
            ##          break
            ##      #print m.group(), m.groups()
            ##      args = eval(m.group(2))
            ##      mu = eval(m.group(3))
            ##      print args, mu
            ##      rho= float(args["rho"])
            ##      T = float(args["T"])
            ##      print args["rho"]
            ##      array[rho][T][species] = mu
            ##      #print >> f, r["T"], r["rho"], mu
            ##      start = m.end()
        else:
            # Get results from the log file using the new output format
            result_rec = re.compile(r'RESULT2:::(\{.*\}):::\n')
            start = 0
            #for R in runs:
            #logfile = file("exp1_02b.log").read()
            while True:
                m = result_rec.search(logfile, start)
                #print m.group()
                if not m:
                    break
                #print m.group(), m.groups()
                d = eval(m.group(1))
                rho= float(d["rho"])
                T = float(d["T"])
                move = d["move"]
                E = Experiment({"T":T, "rho":rho, "move":d["move"]})
                S = E.makeSys()

                mu_0 = d["mu_0"]
                array[rho][T]["mu_0"] = mu_0

                mu_1 = d["mu_1"]
                array[rho][T]["mu_1"] = mu_1

                SP = math.exp(-S.beta*mu_0)
                K_0 = (S.N/S.volume)/(SP*S.beta)
                print S.beta
                array[rho][T]["K_0"] = K_0

                SP = math.exp(-S.beta*mu_1)
                K_1 = (S.N/S.volume)/(SP*S.beta)
                array[rho][T]["K_1"] = K_1

                print rho, T, mu_1, T/mu_1, K_1
                start = m.end()

        def rplot(data, par={}, title={}, xlab="", ylab=""):
            r.plot_new()
            xmax = -1e6
            xmin = 1e6
            ymax = -1e6
            ymin = 1e6
            for x in data:
                print x[0]
                print x[1]
                xmax = max(xmax, max(x[0]))
                xmin = min(xmin, min(x[0]))
                ymax = max(ymax, max(x[1]))
                ymin = min(ymin, min(x[1]))
                print xmin, xmax, ymin, ymax
            r.plot_window(xlim=(xmin, xmax), ylim=(ymin, ymax))
            r.title(
                main=title.get("main", ""),
                xlab=title.get("xlab", ""),
                ylab=title.get("ylab", ""))
            r.box(col="black") ; r.axis(1) ; r.axis(2)
            for i, x in enumerate(data):
                r.lines(x[0], x[1], col=i+1, **par)
                                        
        data = [ ]
        for i, rho in enumerate(sorted(array.iterkeys())):
            l = [ ]
            Legend.add(legend="rho="+str(rho), lty=1, col=i+1)
            for T in sorted(array[rho].iterkeys()):
                l.append((T, array[rho][T][species]))
            data.append(zip(*l))
        rplot(data, par={"type": "l"},
              title={"main": plottitle,
                     "xlab":"T",
                     "ylab":species+" in natural units"})
        Legend.make(x="topleft")
        interact()


    if sys.argv[1] == "test":
        print "building:"
        os.system("sh ./build.sh")
        T = 1.446 / 10.75
        rho = .102 * 5.088448
        #boxedge = 17.2
        boxedge = 10.
        S = dragunov.System(N=int(rho*(boxedge**3)),
                            boxsize=(boxedge, boxedge, boxedge),
                            beta=1./T,
                            trialMoveScale=.5)
        S.dispRadius = 1/1.72
        print S.N
        #S.fill(10, 10, 10, scale=1.7)
        S.fillRandom()
        S.removeOverlaps()
        S.makebox()
        S.display()
        while True:
            S.trialMove(500)
            print S.energy()
            S.display()
