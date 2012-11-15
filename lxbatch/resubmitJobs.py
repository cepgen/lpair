#!/usr/bin/python
from sys import argv, exit
from subprocess import call

if len(argv)==3:
    workDir = argv[1]
    jobsList = argv[2].split(',')
else:
    print "Usage : "+str(argv[0])+" <work directory> <jobs list>"
    exit(1)

for j in jobsList:
    try:
        log = open(workDir+"/logs/"+str(j).zfill(4)+".out")
    except:
        print "[ERROR] log file not found for job "+str(j)+" ! Is the job finished ?"
        exit(2)
    try:
        firstResubmit = call(["ls "+workDir+"/logs/"+str(j).zfill(4)+"_*.log > /dev/null 2>&1"], shell=True)!=0
        if firstResubmit:
            resubNum = 0
    except:
        print "[ERROR] Failed to determine whether the job "+str(j)+" has already been resubmited"
    for line in log.readlines():
        if "Successfully completed" in line:
            print "[ERROR] job "+str(j)+" seems to be successful"
            exit(3)
        if "Exited with exit code" in line:
            exitCode = line.split("code")[1].split('.')[0].strip()
            print "[INFO] job "+str(j)+" failed with exit code "+str(exitCode)

    command = "mv "+workDir+"/logs/"+str(j).zfill(4)+".out "+workDir+"/logs/"+str(j).zfill(4)+"_"+str(resubNum)+".log"
    call(command, shell=True)
    command = "bsub -q 2nd -R \"tmp>50&&mem>200&&swp>400&&pool>1000\" -o "+workDir+"/logs/"+str(j).zfill(4)+".out -e "+workDir+"/logs/"+str(j).zfill(4)+".err `echo sh "+workDir+"/inputs/*_"+str(j).zfill(4)+".sh`"
    call(command, shell=True)
    print "[INFO] Job "+str(j)+" successfully submitted !"
