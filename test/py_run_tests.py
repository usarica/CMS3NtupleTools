import os
import sys
import time
import commands

from multiprocessing import Queue, Process, Pool, Semaphore
from tqdm import tqdm

def do_cmd(cmd):
    commands.getoutput(cmd)
    return True

def follow(fh):
    """
    Iterator generator that returns lines as data is added to the file.
    copied from https://github.com/six8/pytailer/blob/master/src/tailer/__init__.py
    """
    trailing = True
    while 1:
        where = fh.tell()
        line = fh.readline()
        if line:
            if trailing and line in ['\r\n','\n','\r']:
                # This is just the line terminator added to the end of the file
                # before a new line, ignore.
                trailing = False
                continue
            if line[-1] in ['\r\n','\n','\r']:
                line = line[:-1]
                if line[-1:] == '\r\n':
                    # found crlf
                    line = line[:-1]
            trailing = False
            yield line
        else:
            trailing = True
            fh.seek(where)
            return

def nlines_back(n):
    """
    return escape sequences to move character up `n` lines
    and to the beginning of the line
    """
    return "\033[{0}A\r".format(n+1)

class CampaignTest(object):

    def __init__(self, tag=None, exe="cmsRun", pset="main_pset.py", outdir="outputs", **kwargs):
        self.tag = tag
        self.exe = exe
        self.pset = pset
        self.outdir = outdir
        self.kwargs = kwargs
        self.fhlog = None
        self.status = "idle"
        self.total_events = -1
        self.current_event = -1
        self.event_rate = -1
        self.t0 = time.time()

    def get_logname(self):
        return "{outdir}/log_{tag}.txt".format(outdir=self.outdir,tag=self.tag)

    def get_outputname(self):
        return "{outdir}/ntuple_{tag}.root".format(outdir=self.outdir,tag=self.tag)

    def get_cmd(self):
        # return "{} >& {}".format(self.cmd, self.get_logname())
        return "{exe} {pset} {args} output={output} >& {logname}".format(
                exe=self.exe,
                pset=self.pset,
                args=" ".join("{}={}".format(k,v) for k,v in self.kwargs.items()),
                output=self.get_outputname(),
                logname=self.get_logname(),
                )

    def init_log(self):
        if not os.path.exists(self.outdir): os.system("mkdir -p {}".format(self.outdir))
        os.system("rm {logname}; touch {logname};".format(logname=self.get_logname()))
        self.fhlog = open(self.get_logname(),"r")

    def tail_new_lines(self):
        return follow(self.fhlog)

    def get_elapsed(self):
        return time.time() - self.t0

    def update_status(self):
        for line in self.tail_new_lines():
            if "nevents = " in line:
                self.total_events = int(line.split("=")[-1].strip())
                self.status = "started"
            elif "Successfully opened file" in line:
                self.status = "opened input file"
            elif "Begin processing the 1st record" in line:
                self.status = "first event setup"
            elif "Begin processing the" in line:
                if self.status != "looping over events":
                    # reset time to exclude initial overhead. this
                    # wrongly rate because we've thrown out 100 events
                    # (or whatever the report frequency is), but the rate is an
                    # EMA and will average out eventually
                    self.t0 = time.time() 
                    self.status = "looping over events"
                self.current_event = int(line.split(" ",4)[3][:-2])
            elif "Event Throughput" in line:
                self.current_event = self.total_events
                self.event_rate = round(float(line.split(":")[1].strip().split()[0]),2)
                self.status = "done"
                self.t1 = time.time()

    def get_prefix(self):
        return (" "*120)+"\r"+"{:<25}".format(self.tag)

    def get_status_str(self):
        if self.total_events < 0 or self.current_event < 0:
            return self.get_prefix()+self.status
        elif self.status == "done":
            return self.get_prefix()+"done in {:.0f}s @ {:.0f}Hz".format(self.t1-self.t0, (self.total_events-100)/(self.t1-self.t0))
        else:
            # format_meter(n, total, elapsed, ncols=None, prefix='', ascii=False, unit='it', unit_scale=False, rate=None, bar_format=None, postfix=None, unit_divisor=1000)
            return self.get_prefix()+tqdm.format_meter(self.current_event,self.total_events,self.get_elapsed(),unit="evt",postfix=self.status).encode('utf-8')

    def __repr__(self):
        return self.get_cmd()

def run_tests(cts):
    os.nice(10)
    pool = Pool(11)
    [ct.init_log() for ct in cts]
    res = pool.map_async(do_cmd, [ct.get_cmd() for ct in cts])
    pool.close()
    done = False
    print "\n"*(len(cts))
    while not done:
        # vals = res._value
        [ct.update_status() for ct in cts]
        print nlines_back(len(cts))
        print "\n".join([ct.get_status_str() for ct in cts])
        done = res.ready()
        time.sleep(2)
    pool.join()


if __name__ == "__main__":

    nevents = 2000
    cts = [
            CampaignTest(tag="data_2016_80x_v2",
                globaltag="80X_dataRun2_2016SeptRepro_v7",
                nevents=nevents,
                year=2016,
                is80x=True,
                data=True,
                inputs="file:/hadoop/cms/store/user/namin/localcache/data/Run2016F/DoubleMuon/MINIAOD/03Feb2017-v1/100000/201B07A6-57EB-E611-B690-0CC47A57D066.root"),
            CampaignTest(tag="data_2016_94x_v3",
                globaltag="94X_dataRun2_v10",
                nevents=nevents,
                year=2016,
                data=True,
                inputs="file:/hadoop/cms/store/user/namin/localcache/data/Run2016C/DoubleMuon/MINIAOD/17Jul2018-v1/50000/D229CC30-1E8B-E811-844A-A0369FD0B228.root"),
            CampaignTest(tag="data_2017_94x_rereco1",
                globaltag="94X_dataRun2_v11",
                nevents=nevents,
                year=2017,
                data=True,
                metrecipe=True,
                inputs="file:/hadoop/cms/store/user/namin/localcache/data/Run2017C/DoubleMuon/MINIAOD/31Mar2018-v1/80000/04FCFB0D-FF39-E811-94C7-AC162DA6D2F8.root"),
            CampaignTest(tag="data_2017_94x_rereco2",
                globaltag="94X_dataRun2_v11",
                nevents=nevents,
                year=2017,
                data=True,
                inputs="file:/hadoop/cms/store/user/namin/localcache/data/Run2017F/DoubleEG/MINIAOD/09May2018-v1/10000/444E03EB-B75F-E811-AFBA-F01FAFD8F16A.root"),
            CampaignTest(tag="data_2018_102x_rereco1",
                globaltag="102X_dataRun2_Sep2018Rereco_v1",
                nevents=nevents,
                year=2018,
                data=True,
                inputs="file:/hadoop/cms/store/user/namin/localcache/data/Run2018A/DoubleMuon/MINIAOD/17Sep2018-v2/00000/7B954B49-BE06-B64C-89DC-F568513B41A3.root"),
            CampaignTest(tag="data_2018_102x_prompt",
                globaltag="102X_dataRun2_Prompt_v11",
                nevents=nevents,
                year=2018,
                data=True,
                inputs="file:/hadoop/cms/store/user/namin/localcache/data/Run2018D/EGamma/MINIAOD/PromptReco-v2/000/322/204/00000/F09A218A-71B3-E811-9A04-02163E013E33.root"),
            CampaignTest(tag="mc_2016_80x_v2",
                globaltag="80X_mcRun2_asymptotic_2016_TrancheIV_v8",
                nevents=nevents,
                year=2016,
                is80x=True,
                inputs="file:/hadoop/cms/store/user/namin/localcache/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1C8BFCBE-26B6-E611-80E5-A0000420FE80.root"),
            CampaignTest(tag="mc_2016_80x_v2_fastsim",
                globaltag="80X_mcRun2_asymptotic_2016_miniAODv2_v0",
                nevents=nevents,
                year=2016,
                is80x=True,
                fastsim=True,
                inputs="file:/hadoop/cms/store/user/namin/localcache/mc/RunIISpring16MiniAODv2/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/5A5E2A72-4F2D-E611-B2F5-02163E017638.root"),
            CampaignTest(tag="mc_2016_94x_v3",
                globaltag="94X_mcRun2_asymptotic_v3",
                nevents=nevents,
                year=2016,
                inputs="file:/hadoop/cms/store/user/namin/localcache/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/100000/6053E770-34E9-E811-98D2-246E96D10C28.root"),
            CampaignTest(tag="mc_2017_94x_v2",
                globaltag="94X_mc2017_realistic_v17",
                nevents=nevents,
                year=2017,
                metrecipe=True,
                inputs="file:/hadoop/cms/store/user/namin/localcache/mc/RunIIFall17MiniAODv2/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/50000/00436CBC-6B70-E811-850C-00259075D70C.root"),
            CampaignTest(tag="mc_2018_102x_v1",
                globaltag="102X_upgrade2018_realistic_v12",
                nevents=nevents,
                year=2018,
                inputs="file:/hadoop/cms/store/user/namin/localcache/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/042C8EE9-9431-5443-88C8-77F1D910B3A5.root"),
            ]

    run_tests(cts)

    # TODO:
    # actually check the branches to see if values are reasonable/valid
