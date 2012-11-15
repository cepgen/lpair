#!/usr/bin/python
from os import path
from sys import argv, exit
from glob import glob
import subprocess
from array import array
from ROOT import gStyle, TCanvas, TPie, TH1D
### UGLY!! ###
import warnings
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )
##############

runTimesBinsSize = 60

res = list()


def Colourize(string, status='g', bold=True):
    attr = []
    if status=='g':
        attr.append('32')
    elif status=='r':
        attr.append('31')
    if bold:
        attr.append('1')
    return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), string)

def GetRunTimes(dir):
    count = []
    for file in glob(path.join(dir, 'logs', '*.out')):
        for line in reversed(open(file).readlines()):
            if "CPU time" in line:
                count.append(float(line.split(" : ")[1].split('sec')[0].strip()))
                break
    return count

def PlotRunTimes(dir, res):
    if res==[]:
        return False
    name = dir.split('/')[-1]
    if name=='':
        name = dir.split('/')[-2]
    maxVal = max(res)*1.1
    #gStyle.SetOptStat(111111)
    gStyle.SetOptStat(111)
    canv = TCanvas()
    h = TH1D('h'+name, "Running time distribution for "+name, int(maxVal/runTimesBinsSize), 0., maxVal)
    for d in res:
        h.Fill(d)
    h.Draw()
    h.GetXaxis().SetTitle("Running time [s]")
    h.GetYaxis().SetTitle("Number of jobs")
    canv.SetLogy()
    canv.SaveAs('time_'+name+'.png')

def GetJobsStatus(dir):
    for file in glob(path.join(dir, 'logs', '*.out')):
        for line in open(file):
            if "Successfully completed" in line:
                res.append([int(file.split("/")[-1].split(".")[-2]), 0])
                break
            if "Exited with exit code" in line:
                code = int(line.split("code ")[1].split(".")[0])
                res.append([int(file.split("/")[-1].split(".")[-2]), code])
                break
    return res
            
def PlotExitCodes(dir, res):
    if res==[]:
        return False
    name = dir.split('/')[-1]
    if name=='':
        name = dir.split('/')[-2]
    ls = array("d")
    out = []
    codes = []
    codesArr = []
    for it in res:
        if it[1] not in codes:
            codes.append(it[1])
            codesArr.append(array('c',str(it[1])+'\0')) # ugly! thx to pyROOT!
            ls.append(1)
            out.append([it[1], 1])
        else:
            ls[codes.index(it[1])] += 1
            for i in range(0,len(out)):
                if out[i][0]==it[1]:
                    out[i][1] += 1

    canv = TCanvas("canv", "", 400, 400)
    pie = TPie("pie"+dir, "Jobs exit codes for "+name, len(ls), ls)
    pie.SetLabels(array('l', map(lambda x: x.buffer_info()[0], codesArr)))
    if 0 in codes:
        pie.SetEntryFillColor(codes.index(0), 8)
    #pie.SetLabelFormat("%txt (%frac)")
    pie.SetLabelsOffset(-.1)
    pie.Draw("r")
    canv.SaveAs('exitcodes_'+name+".png")
    return out


p = argv[1]
r = GetJobsStatus(p)
numCompleted = len(r)
failedJobs = sorted([res[0] for res in r if res[1]!=0])
numFailed = len(failedJobs)
numSuccess = numCompleted-numFailed
print "Completed jobs : "+str(numCompleted)
print "  "+Colourize(str(numFailed), 'r')+" ("+str(numFailed*100/numCompleted)+"%) failed"
print "  "+Colourize(str(numSuccess), 'g')+" ("+str(numSuccess*100/numCompleted)+"%) successful !"
#print str(numFailed)+" failed jobs ("+str(numFailed*100/numCompleted)+"%) out of "+str(numCompleted)+" completed jobs :"
print ",".join([str(fj) for fj in failedJobs])
e = PlotExitCodes(p, r)
t = GetRunTimes(p)
PlotRunTimes(p, t)
