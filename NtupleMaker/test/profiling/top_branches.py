#!/usr/bin/env python

import argparse
import sys
import os


def get_total_size(br):
    f = r.TMemFile("buffer","CREATE")
    if br.GetTree().GetCurrentFile():
        f.SetCompressionSettings(br.GetTree().GetCurrentFile().GetCompressionSettings())
    f.WriteObject(br,"thisbranch")
    key = f.GetKey("thisbranch");
    basket_size_zip, basket_size_tot = get_basket_size(br)
    return [int(basket_size_zip + key.GetNbytes()), int(basket_size_tot + key.GetNbytes())]

def get_basket_size(br):
    brs = [br]
    subbrs = br.GetListOfBranches()
    if subbrs: brs += subbrs

    ret_zip, ret_tot = 0, 0
    for br in brs:
        ret_zip += br.GetZipBytes()
        ret_tot += br.GetTotBytes()

    return ret_zip, ret_tot


# def main_edm(fname_in, treename, maxnum, precision, groupmakers):
#     # print fname_in



def main(fname_in, treename, maxnum, precision, groupmakers):


    use_edm_tool = True
    use_aliases = False

    from root_utils import get_treename_from_file

    d_bname_to_alias = {}
    if use_aliases:
        import ROOT as r
        f = r.TFile(fname_in)
        if not treename: treename = get_treename_from_file(f)
        tree = f.Get(treename)
        aliases = tree.GetListOfAliases()
        branches = tree.GetListOfBranches()

        for ialias, alias in enumerate(aliases):
            aliasname = alias.GetName()
            branch = tree.GetBranch(tree.GetAlias(aliasname))
            branchname = branch.GetName()
            d_bname_to_alias[branchname.rsplit("obj",1)[0]] = aliasname

    if use_edm_tool:
        d_info = {}
        events = -1
        uid = os.getuid()

        import ROOT as r
        f = r.TFile(fname_in)
        if not treename: treename = get_treename_from_file(f)

        tmp_name = "/tmp/branchsizes_{}_{}.txt".format(uid,str(hash(fname_in))[:10])
        os.system("edmEventSize -n {} -a -d {} -o {} >& /dev/null".format(treename,fname_in,tmp_name))
        with open(tmp_name, "r") as fhin:
            line1 = fhin.readline()
            events = int(line1.split()[-1])
            fhin.readline()
            for line in fhin:
                line = line.strip()
                if not line: continue
                name, uncomp_bpere, comp_bpere = line.split()
                if use_aliases and not groupmakers:
                    name = d_bname_to_alias.get(name,name)
                uz = float(uncomp_bpere)*events
                z = float(comp_bpere)*events
                if groupmakers and "_" in name:
                    name = name.split("_")[1]
                if name not in d_info: d_info[name] = [z,uz]
                else: d_info[name] = [z+d_info[name][0],uz+d_info[name][1]]
    else:
        ## Bugged. Doesn't show pfcands p4 which is largest contribution
        import ROOT as r
        f = r.TFile(fname_in)
        tree = f.Get(treename)
        events = tree.GetEntries()
        d_info = {}
        for br in tree.GetListOfBranches():
            z,uz = get_total_size(br) # zipped and unzipped sizes
            name = br.GetName()
            # print name,z,uz
            if groupmakers and "_" in name:
                name = name.split("_")[1]
            if name not in d_info: d_info[name] = [z,uz]
            else: d_info[name] = [z+d_info[name][0],uz+d_info[name][1]]

    tot_z = sum(x[0] for x in d_info.values())
    # tot_uz = sum(x[1] for x in d_info.values())
    branches = []
    for n, (z,uz) in d_info.items():
        branches.append({"bname": n, "frac": 1.0*z/tot_z, "uncompBytes": uz, "compBytes": z})

    sizeperevt = int(os.stat(fname_in).st_size / events)

    print
    print " Tree name: %s" % treename
    print "      nevents: %i" % events
    print "      bytes/evt: %i" % sizeperevt
    print
    top_branches = sorted(branches, key=lambda x: x.get("frac",-1), reverse=True)[:maxnum]
    maxcols = max([len(b["bname"]) for b in top_branches[:maxnum]])+4
    print ("%-{0}s %s".format(maxcols)) % ("branchname", "size (%) [compression factor]")
    print "-" * (maxcols+15)
    for b in top_branches:
        print ("%-{0}s %2.{1}f [%03.1f]".format(maxcols, precision)) % (b["bname"], 100.0*b["frac"], 1.0*b["uncompBytes"]/b["compBytes"])

if __name__ == "__main__":
    pass

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="name of file to make classfile on")
    parser.add_argument("-t", "--tree", help="treename (default: Events)", default=None)
    parser.add_argument("-n", "--num", help="number of top branches to show", default=30)
    parser.add_argument("-p", "--precision", help="number of decimal places", default=1)
    parser.add_argument("-m", "--groupmakers", help="group branches from same makers together", default=False, action="store_true")
    args = parser.parse_args()
    fname_in = args.filename
    treename = args.tree
    maxnum = int(args.num)
    main(fname_in, treename, maxnum, int(args.precision), args.groupmakers)

    # edmEventSize -n Events -a -d /hadoop/cms/store/group/snt/run2_moriond17/WZ_TuneCUETP8M1_13TeV-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/merged_ntuple_1.root
    # main_edm("/hadoop/cms/store/group/snt/run2_moriond17/WZ_TuneCUETP8M1_13TeV-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/merged_ntuple_1.root", "Events", 30, 1, False)
    # main("/hadoop/cms/store/group/snt/run2_moriond17/WZ_TuneCUETP8M1_13TeV-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/V08-00-16/merged_ntuple_1.root", "Events", 30, 1, True)
