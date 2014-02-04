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
                   default = 'DataJob',
                   help="Base output directory for the job. Today's date will be appended (eg TestJob_Jul11).")
parser.add_option ('--jdl', dest='jdlBase', type='string',
                   default = 'scripts/go_data.jdl',
                   help="Condor jdl base. Only change this if you know what you're doing!")
parser.add_option ('--script', dest='scriptBase', type='string',
		   default = 'scripts/go_data.sh',
		   help="Condor executable script base. Only change this if you know what you're doing!")
parser.add_option ('--ana', dest='anaBase', type='string',
		   default = 'scripts/ana_data.C',
		   help="Root macro, your ana.C for running over data.")
parser.add_option ('--njobs', dest='njobs', type='int',
                   default = '1',
                   help="Number of jobs, default = 1 job for all files. A value of -1 will submit one job per file. Otherwise, data nTuples will be split as evenly as possibly among this many jobs.")
parser.add_option ('--test', action="store_true",
                   dest="test", default=False,
                   help="Create the job, but don't submit it. For testing purposes.")
parser.add_option ('--json', dest='json', type='string',
                   default = 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt',
                   help="JSON file used to select events. Default is the 22Jan2013 re-reco.")
parser.add_option ('--inputFolders', dest='inputFolders', type='string',
		   default='',
		   help='Colon-separated list of folders containing nTuple files to run over. Using ALL will give the full 22Jan2013 DoublePhoton re-reco (tag cms538v1).')

options, args = parser.parse_args()

baseOutdir = options.baseOutdir
jdlBase = options.jdlBase
scriptBase = options.scriptBase
anaBase = options.anaBase
njobs = options.njobs
test = options.test
json = options.json
inputFolders = options.inputFolders

baseOutdir = baseOutdir+"_"+datetime.datetime.now().strftime("%b%d")

if inputFolders == 'ALL':
    inputFolders = '/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012A-22Jan2013-v1/Photon/:/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/:/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012C-22Jan2013-v2/DoublePhoton/:/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012D-22Jan2013-v1/DoublePhoton/'

cwd = os.getcwd()

if not os.path.isdir(baseOutdir) :
    os.system("mkdir -p "+baseOutdir)
    print "Making directory %s." % baseOutdir
else:
    print "Output directory %s already exists.  Exiting." % baseOutdir
    sys.exit()

allFiles = []

if len(inputFolders) != 0:
    for folder in inputFolders.split(':'):
        if folder[-1] != '/':
	    folder = folder+'/'
	allFiles += glob.glob(folder+"*.root")
else:
    print "Zero files given, already done! Exiting."
    sys.exit()

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

    filelistname = "filelist_"+jobid
    filelist = open(baseOutdir+"/"+filelistname, 'w')

    for ifile in files_for_jobs[ijob]:
        filelist.write(ifile+"\n")
    filelist.close()

if not os.path.isdir(baseOutdir+"/JobOut"): os.system("mkdir "+baseOutdir+"/JobOut")

name_jdl = jdlBase.replace('scripts/', '')
name_script = scriptBase.replace('scripts/', '')
name_ana = anaBase.replace('scripts/', '')

path_jdl = baseOutdir+"/"+name_jdl
path_script = baseOutdir+"/"+name_script
path_ana = baseOutdir+"/"+name_ana

commandList = []
commandList.append("cp "+jdlBase+" "+baseOutdir)
#commandList.append("cp "+json+" "+baseOutdir)
commandList.append("cp "+scriptBase+" "+baseOutdir)
commandList.append("cp "+anaBase+" "+baseOutdir)

commandList.append('replace NJOBS  '+str(njobs)+' -- '+path_jdl)
commandList.append('replace SCRIPT  '+name_script+' -- '+path_jdl)
commandList.append('replace ANALYZER  '+name_ana+' -- '+path_jdl)
commandList.append('replace ANALYZER  '+name_ana+' -- '+path_script)

commandList.append('replace JSON  '+json+' -- '+path_jdl)
commandList.append('replace JSON  '+json+' -- '+path_script)
commandList.append('replace JSON  '+json+' -- '+path_ana)

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

    os.system("condor_submit "+name_jdl)
    os.chdir(cwd)





