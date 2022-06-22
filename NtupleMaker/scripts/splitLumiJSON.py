#!/bin/env python

import json
import sys
import argparse



def run(args):
   input_json_file = args.json
   first_run = args.first_run
   last_run = args.last_run
   incruns = args.include_run
   excruns = args.exclude_run
   strincruns=""
   strexcruns=""
   strfirstrun=""
   strlastrun=""
   if incruns is not None and len(incruns)>0:
      strincruns="_inc_{}".format('_'.join([str(x) for x in incruns]))
   if excruns is not None and len(excruns)>0:
      strexcruns="_exc_{}".format('_'.join([str(x) for x in excruns]))
   if first_run>0:
      strfirstrun="_firstRun_{}".format(first_run)
   if last_run>0:
      strlastrun="_lastRun_{}".format(last_run)

   if "".join([strincruns,strexcruns,strfirstrun,strlastrun])=="":
      raise RuntimeError("The JSON file cannot just be copied.")

   output_json_file = "{}{}{}{}{}.txt".format(input_json_file.replace(".txt","").replace(".json",""), strfirstrun, strlastrun, strincruns, strexcruns)

   json_dict = None
   with open(input_json_file) as fin:
      json_dict = json.load(fin)

   new_json_dict=dict()
   for dd in json_dict:
      run_number = int(dd)
      if excruns is not None and run_number in excruns:
         continue
      if ((first_run<0 or run_number>=first_run) and (last_run<0 or run_number<=last_run)) or (incruns is not None and run_number in incruns):
         new_json_dict[dd] = json_dict[dd]

   with open(output_json_file, 'w') as fout:
      json.dump(new_json_dict, fout)



if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument("--json", type=str, help="JSON file to split", required=True)
   parser.add_argument("--first_run", type=int, help="First run to include", required=False, default=-1)
   parser.add_argument("--last_run", type=int, help="Last run to include", required=False, default=-1)
   parser.add_argument("--include_run", type=int, help="Include additional runs outside of the specified range", required=False, action='append')
   parser.add_argument("--exclude_run", type=int, help="Exclude runs within the specified range", required=False, action='append')
   args = parser.parse_args()
   run(args)
