#!/usr/bin/env python

import argparse
import os
import glob
from _ctypes import sizeof

def getListOfSets():
    list  = [name for name in os.listdir('.') if (os.path.isdir(os.path.join('.', name)) and 'Set_' in name)]
    return sorted(list, key = lambda a: int(a.split("_")[1]))

parser = argparse.ArgumentParser(
                                 description='Check if root file exists in each Set and resubmit if necessary.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--className', dest="className", default="*", help="Name of the used analysis class to identify the root input files.")
parser.add_argument('--uncertType', dest="uncertType", default="default", help="Uncertainty type to identify the root input files.")
parser.add_argument('--noResubmit', action='store_true', help="Update jobs_submittedOrComplete, but do not run Submit again.")
args = parser.parse_args()

filename = "GRID_ANALYSIS_" + args.className + "_" + args.uncertType + ".root"
nFiles = -1 # keep track of how many root files are found

badSets = []

# Check if all Sets contain the root files
lSets = getListOfSets()
for set in lSets:
    nFilesInSet = len(glob.glob(set + "/" + filename))
    
    if nFiles == -1:
        # determine expected number of root files from first set
        if nFilesInSet == 0:
            print "WARNING:", set, "contains", nFilesInSet, "root files. Using next set to determine expected number of root files."
        else:
            nFiles = nFilesInSet
    
    if nFilesInSet == 0:
        badSets.append(set)
    elif nFilesInSet != nFiles:
        print "WARNING:", set, "contains", nFilesInSet, "root files, while", nFiles, "are expected. Please take care of this yourself."

# Delete bad Sets from log file
logFile = open("jobs_submittedOrComplete","r") # read only
lines = logFile.readlines()
logFile.close()

logFile = open("jobs_submittedOrComplete","w") # overwrite
for line in lines:
    s = line.split("/")[-1].strip()
    if s not in badSets:
        logFile.write(line)
logFile.close()

print len(badSets), "sets do not contain any root files and have been removed from jobs_submittedOrComplete:"
print badSets

if not args.noResubmit:
    print "Resubmit failed sets..."
    os.system("source Submit --Submit")
else:
    print "To resubmit the jobs, run 'source Submit --Submit'"
