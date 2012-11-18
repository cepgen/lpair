#!/usr/bin/python
from sys import argv, exit
from subprocess import call
from os import listdir

forceResubmit = False
quietMode = False

def main(argv):
    if len(argv)>2:
        workDir = argv[1]
        jobsList = argv[2].split(',')
        if len(argv)>3:
            if argv[3]=='1':
                forceResubmit = True
            if len(argv)>4:
                if argv[4]=='1':
                    quietMode = True
        print jobsList
        print "\n"
    else:
        print "Usage : "+str(argv[0])+" <work directory> <jobs list> [<force resubmit> <quiet mode>]"
        exit(1)

    for j in jobsList:
        try:
            log = open(workDir+"/logs/"+str(j).zfill(4)+".out")
        except:
            if not forceResubmit:
                if not quietMode:
                    print "[ERROR] Log file not found for job "+str(j)+" ! Is the job finished ?"
                exit(2)
            elif not quietMode:
                print "[INFO] Log file not found for job "+str(j)+". Forced resubmit"
        try:
            firstResubmit = call(["ls "+workDir+"/logs/"+str(j).zfill(4)+"_*.log > /dev/null 2>&1"], shell=True)!=0
            if firstResubmit:
                print "First resubmit for job "+str(j)+" !"
                resubNum = 0
            else:
                listOfResubmits = []
                for f in listdir(workDir+"/logs/"):
                    if str(j).zfill(4)+"_" in f:
                        num = f.split('_')[-1].split('.')[-2]
                        listOfResubmits.append(int(num))
                resubNum = max(listOfResubmits)+1
        except:
            if not quietMode:
                print "[ERROR] Failed to determine whether the job "+str(j)+" has already been resubmited"
            exit(3)

        for line in log.readlines():
            if "Successfully completed" in line:
                if not quietMode:
                    print "[ERROR] Job "+str(j)+" seems to be successful"
                exit(3)
            if "Exited with exit code" in line:
                exitCode = line.split("code")[1].split('.')[0].strip()
                if not quietMode:
                    print "[INFO] Job "+str(j)+" failed with exit code "+str(exitCode)

        command = "rm "+workDir+"/logs/"+str(j).zfill(4)+".err"
        call(command, shell=True)
        command = "mv "+workDir+"/logs/"+str(j).zfill(4)+".out "+workDir+"/logs/"+str(j).zfill(4)+"_"+str(resubNum)+".log"
        call(command, shell=True)
        queue = "8nh"
        if "inelinel" in workDir:
            queue = "1nd"
        command = "bsub -q "+queue+" -R \"tmp>50&&mem>200&&swp>400&&pool>1000\" -o "+workDir+"/logs/"+str(j).zfill(4)+".out -e "+workDir+"/logs/"+str(j).zfill(4)+".err `echo sh "+workDir+"/inputs/*_"+str(j).zfill(4)+".sh`"
        call(command, shell=True)
        if not quietMode:
            print "[INFO] Job "+str(j)+" successfully submitted !"

if __name__ == "__main__":
        main(argv)
