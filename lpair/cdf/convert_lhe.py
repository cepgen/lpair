#!/bin/python

#infile = open('events.mumu.el-el.pt15.m0to500.ascii')
infile = open('events.tautau.el-el.pt15.m0to500.ascii')

isProton1 = False
isProton2 = False
isDiMuon = False
isLepton1 = False
isLepton2 = False

beamEnergy = 3500.

class Event:
    def __init__(self, arr):
        self._eventId = int(arr[0])
        self._numParticles = int(arr[1])
        self._numParticlesAdded = 0
        self.Particles = []
    def AddParticle(self, part):
        self.Particles.append(part)
        self._numParticlesAdded += 1
        if part.Type == 1:
            self.InP1 = part
        elif part.Type == 2:
            self.InP2 = part
        elif part.Type == 3:
            self.OutP1 = part
        elif part.Type == 4:
            self.DiMuon = part
        elif part.Type == 5:
            self.OutP2 = part
        elif part.Type == 6:
            part.SetPDGId(15)
            self.Lepton1 = part
        elif part.Type == 7:
            part.SetPDGId(-15)
            self.Lepton2 = part
            #self.CloseEvent()
            return "closed"
    def CloseEvent(self):
        if self._eventId%1000==0:
            print "Processing event #"+str(self._eventId)
            print "  --> "+str(self._numParticles)+"/"+str(self._numParticlesAdded)+" particles"
        

"""    def AddParticle(self, arr):
        part = Particle(arr[2])
        if arr[0]=='3': # outgoing proton 1
            part.SetP()
            self.proton1 = part"""

class Particle:
    def __init__(self, arr):
        self.P = [0., 0., 0.]
        self.P4 = [0., 0., 0., 0.]
        self.M = 0.
        self.E = 0.
        self.SetPDGId(arr[2])
        self.SetType(arr[0])
    def SetType(self, typ):
        self.Type = int(typ)
    def SetP(self, arr):
        self.P = [float(arr[0]), float(arr[1]), float(arr[2])]
        self.Px = float(arr[0])
        self.Py = float(arr[1])
        self.Pz = float(arr[2])
    def SetP4(self, arr):
        self.P4 = [float(arr[0]), float(arr[1]), float(arr[2]), float(arr[3])]
        self.SetP(arr[0:3])
        self.SetE(arr[3])
    def SetE(self, energy):
        self.E = float(energy)
    def SetM(self, mass):
        self.M = float(mass)
    def SetPDGId(self, pdgId):
        self.PDGId = int(pdgId)
    def PrintLHE(self):
        print "\t"+str(self.PDGId)+"\t\t1\t1\t2\t0\t0\t"+str(self.Px)+"\t"+str(self.Py)+"\t"+str(self.Pz)+"\t"+str(self.E)+"\t"+str(self.M)+" 0.\t1."


null  = ['0.', '0.', '0.', '0.']
evts = []

print "<LesHouchesEvents version=\"1.0\">"
print "<header>"
print "This file was created from the output of the LPAIR generator"
print "</header>"
print "<init>"
print "2212  2212  0.35000000000E+04  0.35000000000E+04 0 0 10042 10042 2  1"
print "0.10508723460E+01  0.96530000000E-02  0.26731120000E-03   0"
print "</init>"
for line in infile:
    ln = line.split()
    if len(ln) == 2:
        ev = Event(ln) # new event invoked
        print "<event>"
        print "\t6\t0\t0.2983460E-04\t0.9118800E+02\t0.7546772E-02\t0.1300000E+00"
        #print "\t6\t661\t0.2983460E-04\t.9118800E+02\t0.7821702E-02\t0.1300000E+00"
        continue
    if ln == null:
        continue
    if len(ln) == 7: # particle definition
        part = Particle(ln)
        continue
    elif len(ln) == 5: # kinematic information
        part.SetP4(ln[0:4])
        part.SetM(ln[4])
        if ev.AddParticle(part) == 'closed':
            # photons
            print "\t22\t-1\t0\t0\t0\t0\t0.00000000000E+00\t0.00000000000E+00\t0.00000000000E+02\t0.10000000000E+02\t0.00000000000E+00\t0.\t1."
            print "\t22\t-1\t0\t0\t0\t0\t0.00000000000E+00\t0.00000000000E+00\t0.00000000000E+00\t0.10000000000E+02\t0.00000000000E+00\t0.\t-1."
            ev.OutP1.PrintLHE()
            ev.OutP2.PrintLHE()
            ev.Lepton1.PrintLHE()
            ev.Lepton2.PrintLHE()
            print ""+str()
            print "</event>"
print "</LesHouchesEvents>"
