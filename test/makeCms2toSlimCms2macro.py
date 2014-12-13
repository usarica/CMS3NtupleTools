#! /usr/bin/env python

import sys,re,os

CMSSWpath = os.environ['CMSSW_BASE']
infilename = CMSSWpath+'/src/CMS3/NtupleMaker/python/SlimCms2_cff.py'
outfilename = 'Cms2toSlimCms2macro.C'

if len(sys.argv) < 1:
    print 'Usage: makeCms2toSlimCms2config.py [OPTIONS]\n'
    print 'where the required options are:'
    print '\t--infile\tinput SlimCms2_cff.py file\n'
    print 'Optional arguments:'
    print '\t--macro\tname of root macro to output; default is Cms2toSlimCms2macro.C'
    sys.exit()
elif len(sys.argv) > 5:
    print 'Usage: makeCms2toSlimCms2config.py [OPTIONS]\n'
    print 'where the required options are:'
    print '\t--infile\input SlimCms2_cff.py file\n'
    print 'Optional arguments:'
    print '\t--macro\tname of root macro to output; default is Cms2toSlimCms2macro.C'    
    sys.exit()

for i in range(0,len(sys.argv)):
    if sys.argv[i] == '--infile':
        infilename = sys.argv[i+1]
    if sys.argv[i] == '--macro':
        outfilename = sys.argv[i+1]

infile = open(infilename, 'r')
outfile = open(outfilename, 'w')

outfile.write('#include "TFile.h"\n')
outfile.write('#include "TTree.h"\n')
outfile.write('#include <string>\n')
outfile.write('#include <iostream>\n\n')

outfile.write('int Cms2toSlimCms2macro(std::string ifname)\n')
outfile.write('{\n')

outfile.write('\tTFile *oldfile = TFile::Open(ifname.c_str());\n');
outfile.write('\tif (not oldfile)\n')
outfile.write('\t{\n')
outfile.write('\t\tstd::cout << "Could not open input file.  Exiting." << std::endl;\n')
outfile.write('\t\treturn 1;\n')
outfile.write('\t}\n\n')

outfile.write('\tTTree *oldtree = (TTree*)oldfile->Get("Events");\n')
outfile.write('\tif (not oldtree)\n')
outfile.write('\t{\n')
outfile.write('\t\tstd::cout << "Input file does not contain a good tree.  Exiting." << std::endl;\n')
outfile.write('\t\treturn 2;\n')
outfile.write('\t}\n\n')

outfile.write('\toldtree->SetBranchStatus("*", 1);\n')

for line in infile:
    if line[0] == '#':
        continue
    substr = re.search("\('drop (.*)'\)", line)
    if substr == None:
        continue
    outfile.write('\toldtree->SetBranchStatus("')
    outfile.write(substr.group(1))
    outfile.write('", 0);\n')    

outfile.write('\n')
outfile.write('\tstd::string ofname = ifname;\n')
outfile.write('\tofname.erase(ofname.length()-5,5);\n')
outfile.write('\tofname.append("_slim.root");\n')
outfile.write('\tTFile *newfile = TFile::Open(ofname.c_str(),"recreate");\n')
outfile.write('\tTTree *newtree = oldtree->CloneTree(0);\n')
outfile.write('\tnewtree->CopyEntries(oldtree);\n')
outfile.write('\tnewfile->Write();\n')
outfile.write('\tdelete oldfile;\n')
outfile.write('\tdelete newfile;\n')
outfile.write('\treturn 0;\n')
outfile.write('}\n')
