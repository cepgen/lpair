#!/bin/python

#infile = open('events.mumu.inel-inel.pt15.m0to500.ascii')
infile = open('events.tautau.inel-inel.pt15.m0to500.ascii')

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
    def AddParticle(self, part):
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
            self.Mum = part
        elif part.Type == 7:
            self.Mup = part
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
        

null  = ['0.', '0.', '0.', '0.']
evts = []

for line in infile:
    ln = line.split()
    if len(ln) == 2:
        ev = Event(ln) # new event invoked
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
            outstr = str(beamEnergy-ev.OutP1.E)+" "+str(beamEnergy-ev.OutP2.E)+" "
            outstr+= str(ev.OutP1.Px)+" "+str(ev.OutP1.Py)+" "+str(ev.OutP1.Pz)+" "
            outstr+= str(ev.OutP2.Px)+" "+str(ev.OutP2.Py)+" "+str(ev.OutP2.Pz)+" "
            outstr+= str(ev.Mum.Px)+" "+str(ev.Mum.Py)+" "+str(ev.Mum.Pz)+" "+str(ev.Mum.E)+" "
            outstr+= str(ev.Mup.Px)+" "+str(ev.Mup.Py)+" "+str(ev.Mup.Pz)+" "+str(ev.Mup.E)+" "
            outstr+= str(ev.DiMuon.Pz)+" "+str(ev.DiMuon.E)
            print outstr
    #print ln


"""    if ln[1]=='1': # particle type tag
        if ln[0]=='3':
            isProton1 = True
            Pp1 = [float(ln[1]), float(ln[2]), float (ln[3])]
            Eg1 = beamEnergy-float(ln[4])
            print Eg1
            continue
        elif ln[0]=='4':
            isDiMuon = True
            continue
        elif ln[0]=='5':
            isProton2 = True
            continue
        elif ln[0]=='6':
            isLepton1 = True
            continue
        elif ln[0]=='7':
            isLepton2 = True
            continue
    else:
        if isProton1:
            #print ln
            #print Eg1
            isProton1 = False         
"""
