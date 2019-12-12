## Grid submission

1. Checkout ProjectMetis and make sure it is in the PATH (via the `setup.sh` script):
```
git clone https://github.com/aminnj/ProjectMetis
cd ProjectMetis
source setup.sh
```
2. Edit `samples_*.csv` appropriately. And after making sure everything compiles properly, make a tarball for the worker node:
```bash
mtarfile tarball_v1.tar.xz --xz --xz_level 3 -x "ZZMatrixElement/MELA/data/Pdfdata" "*ZZMatrixElement/MELA/data/*.root"
```
3. Edit `submit_jobs.py` to consider the right samples. In particular, update the `tarfile`  variable
and `tag` to uniquely identify the submission campaign.
4. Run `python submit_jobs.py` in a screen to (re)submit jobs.
5. Visit the monitoring page to view progress and output location.
