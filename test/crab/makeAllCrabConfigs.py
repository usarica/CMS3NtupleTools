
# interacting with the os
import subprocess
from subprocess import Popen, PIPE, STDOUT
from optparse import OptionParser

# get the input file name
optParser = OptionParser()
optParser.add_option("-f", "--file", dest="file",
                  help="Input file containing dataset names", metavar=" datasets.txt")
(options, args) = optParser.parse_args()
if options.file == None:
        optParser.error("specify the -f argument!")

datasets = open(options.file, "r")
tag = "V02-00-12"
config = "../testSingle_cfg.py"

for dataset in datasets:
	print "doing... ", dataset.rstrip("\n")
	cmd = "../../../NtupleMacros/NtupleTools/makeCrabFiles.py -CMS2cfg "
	cmd += config + " -d " + dataset.rstrip("\n") + " -t " + tag
	p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
	print p.stdout.read()

datasets.close()

