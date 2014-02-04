#!/usr/bin/env python

import sys
import optparse
import commands
import os
import glob
import time
import datetime
import tarfile

#######################
# Get options
#######################

parser = optparse.OptionParser("usage: %prog [options]")

parser.add_option ('--o', dest='baseOutdir', type='string',
                   default = 'MCJob',
                   help="Base output directory for the job. Today's date will be appended (eg TestJob_Jul11).")
parser.add_option ('--jdl', dest='jdlBase', type='string',
                   default = 'scripts/go_mc.jdl',
                   help="Condor jdl base. Only change this if you know what you're doing!")
parser.add_option ('--script', dest='scriptBase', type='string',
		   default = 'scripts/go_mc.sh',
		   help="Condor executable script base. Only change this if you know what you're doing!")
parser.add_option ('--ana', dest='anaBase', type='string',
		   default = 'scripts/ana_mc.C',
		   help="Root macro, your ana.C for running over data.")
parser.add_option ('--test', action="store_true",
                   dest="test", default=False,
                   help="Create the job, but don't submit it. For testing purposes.")

parser.add_option ('--name', dest='datasetName', type='string',
                   default = 'bkgMC',
                   help="Name of the MC sample. Output ROOT files and internal trees will be named by this to keep track of things afterwards.")

parser.add_option ('--inputFolder', dest='inputFolder', type='string',
		   default='',
		   help='Path to input ROOT files.')

options, args = parser.parse_args()

baseOutdir = options.baseOutdir
jdlBase = options.jdlBase
scriptBase = options.scriptBase
anaBase = options.anaBase
njobs = options.njobs
test = options.test
datasetName = options.datasetName
inputFolder = options.inputFolder
runStaged = options.runStaged

baseOutdir = baseOutdir+"_"+datetime.datetime.now().strftime("%b%d")

cwd = os.getcwd()

if not os.path.isdir(baseOutdir) :
    os.system("mkdir -p "+baseOutdir)
    print "Making directory %s." % baseOutdir
else:
    print "Output directory %s already exists.  Exiting." % baseOutdir
    sys.exit()

allFiles = []

if len(inputFolder) == 0:
    print "Zero files given, already done! Exiting."
    sys.exit()
else:
    if inputFolder[-1] != '/':
        inputFolder = inputFolder+'/'
    allFiles += glob.glob(inputFolder+"*.root")

nfiles = len(allFiles)
# If -1, make 1 job per file:
if njobs == -1: njobs = nfiles

# Split files into jobs as evenly as possible:

baseFilesPerJob = nfiles/njobs
extraFiles = nfiles%njobs

numberOfFiles = {}
for ijob in range(njobs):
    numberOfFiles[ijob] = baseFilesPerJob

for ijob in range(extraFiles):
    numberOfFiles[ijob] += 1


files_for_jobs = {}
IFILE = 0
for ijob in range(njobs):
    nfiles_ijob = numberOfFiles[ijob]
    files_for_jobs[ijob] = []
    for ifile in range(nfiles_ijob):
        files_for_jobs[ijob].append( allFiles [IFILE] )
        IFILE += 1

if IFILE != nfiles:
    print "Mismatch.", IFILE, nfiles-1,"  Exiting."
    sys.exit()

for ijob in range(njobs):

    jobid = str(ijob)
#    if ijob == njobs-1:
#	print files_for_jobs[ijob]

    filelistname = "filelist_"+jobid
    filelist = open(baseOutdir+"/"+filelistname, 'w')
	
    for ifile in files_for_jobs[ijob]:
        filelist.write(ifile+"\n")
    filelist.close()

if not os.path.isdir(baseOutdir+"/JobOut"): os.system("mkdir "+baseOutdir+"/JobOut")

name_jdl = jdlBase.replace('scripts/', '')
name_script = scriptBase.replace('scripts/', '')
name_ana = anaBase.replace('scripts/', '')

commandList = []
commandList.append("cp "+jdlBase+" "+baseOutdir)
commandList.append("cp "+scriptBase+" "+baseOutdir)

if njobs != 1 or runStaged:
    commandList.append("cp "+jdlBase+" "+baseOutdir+"/pileup_"+name_jdl)
    commandList.append('replace STAGING  pileup -- '+baseOutdir+"/pileup_"+name_jdl)
    commandList.append('replace NJOBS  '+str(njobs)+' -- '+baseOutdir+"/pileup_"+name_jdl)
    commandList.append('replace SCRIPT  '+name_script+' -- '+baseOutdir+"/pileup_"+name_jdl)
    commandList.append('replace ANALYZER  '+name_ana+' -- '+baseOutdir+"/pileup_"+name_jdl)

    commandList.append("cp "+jdlBase+" "+baseOutdir+"/btag_"+name_jdl)
    commandList.append('replace STAGING  btag -- '+baseOutdir+"/btag_"+name_jdl)
    commandList.append('replace NJOBS  '+str(njobs)+' -- '+baseOutdir+"/btag_"+name_jdl)
    commandList.append('replace SCRIPT  '+name_script+' -- '+baseOutdir+"/btag_"+name_jdl)
    commandList.append('replace ANALYZER  '+name_ana+' -- '+baseOutdir+"/btag_"+name_jdl)

    commandList.append("cp "+jdlBase+" "+baseOutdir+"/acceptance_"+name_jdl)
    commandList.append('replace ANALYZER  ANALYZER,pileupReweighting_'+datasetName+'.root,btagEfficiency_'+datasetName+'.root -- '+baseOutdir+"/acceptance_"+name_jdl)
    commandList.append('replace STAGING  acceptance -- '+baseOutdir+"/acceptance_"+name_jdl)
    commandList.append('replace NJOBS  '+str(njobs)+' -- '+baseOutdir+"/acceptance_"+name_jdl)
    commandList.append('replace SCRIPT  '+name_script+' -- '+baseOutdir+"/acceptance_"+name_jdl)
    commandList.append('replace ANALYZER  '+name_ana+' -- '+baseOutdir+"/acceptance_"+name_jdl)

    commandList.append("cp scripts/combineStagedInput.sh "+baseOutdir+"/combineStagedInput.sh")
    commandList.append('replace DATASETNAME  '+datasetName+' -- '+baseOutdir+"/combineStagedInput.sh")

    commandList.append("cp scripts/combineStagedOutput.sh "+baseOutdir+"/combineStagedOutput.sh")
    commandList.append('replace DATASETNAME  '+datasetName+' -- '+baseOutdir+"/combineStagedOutput.sh")

    commandList.append("cp scripts/combineHistograms.C "+baseOutdir+"/combineHistograms.C")

else:
    commandList.append("cp "+jdlBase+" "+baseOutdir+"/"+name_jdl)
    commandList.append('replace STAGING  all -- '+baseOutdir+"/"+name_jdl)
    commandList.append('replace NJOBS  '+str(njobs)+' -- '+baseOutdir+"/"+name_jdl)
    commandList.append('replace SCRIPT  '+name_script+' -- '+baseOutdir+"/"+name_jdl)
    commandList.append('replace ANALYZER  '+name_ana+' -- '+baseOutdir+"/"+name_jdl)

commandList.append("cp "+scriptBase+" "+baseOutdir)
commandList.append('replace DATASETNAME  '+datasetName+' -- '+baseOutdir+"/"+name_script)
commandList.append('replace ANALYZER  '+name_ana+' -- '+baseOutdir+"/"+name_script)

commandList.append("cp "+anaBase+" "+baseOutdir)
commandList.append('replace DATASETNAME  '+datasetName+' -- '+baseOutdir+"/"+name_ana)

for command in commandList:
    os.system(command)

if not test:

    os.system('./scripts/remakeTarball.sh')
    
    os.chdir(baseOutdir)

    tar = tarfile.open('fileLists.tgz', 'w:gz')
    jIter = 0
    for job in range(njobs):
        tar.add('filelist_'+str(jIter))
        jIter += 1
    tar.close()

    if njobs != 1 or runStaged:
        os.system("condor_submit pileup_"+name_jdl)
        os.system("condor_submit btag_"+name_jdl)
        os.system("cat ../doc/staged_mc_warning.txt")
    else:
        os.system("condor_submit "+name_jdl)

    os.chdir(cwd)





