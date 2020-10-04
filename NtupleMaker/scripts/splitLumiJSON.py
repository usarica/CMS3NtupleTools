#!/bin/env python

import json
import sys


if len(sys.argv) < 3:
   print("  usage: {0} <JSON file> <first run> <last run>".format(sys.argv[0]))
   exit(0)


input_json_file = sys.argv[1]
first_run = int(sys.argv[2])
last_run = int(sys.argv[3])
output_json_file = "{}_firstRun_{}_lastRun_{}.txt".format(input_json_file.replace(".txt","").replace(".json",""), first_run, last_run)

json_dict = None
with open(input_json_file) as fin:
   json_dict = json.load(fin)

new_json_dict=dict()
for dd in json_dict:
   run_number = int(dd)
   if (first_run<0 or run_number>=first_run) and (last_run<0 or run_number<=last_run):
      new_json_dict[dd] = json_dict[dd]

with open(output_json_file, 'w') as fout:
    json.dump(new_json_dict, fout)

