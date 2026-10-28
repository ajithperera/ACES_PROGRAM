#!/bin/bash
# OO-pCCD slow tier. From the directory that should get the report:
#   sbatch ~/Develop/ACES_PROGRAM/aces2/tests/pccd_slow/submit.sbatch [--long] [--case NAME ...]
# The default set takes ~15 min on 10 cores; --long adds C2H4 cc-pVTZ
# SYMMETRY=OFF (~1 h). Work directories go to node-local $TMPDIR; failed
# cases leave run.log in --workdir if one is given.
#SBATCH --job-name=pccd_slow
#SBATCH --output=pccd_slow_%j.out
#SBATCH --partition=hpg-milan
#SBATCH --account=ufhpc
#SBATCH --qos=ufhpc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=80gb
#SBATCH --time=4:00:00
module load intel/2025.1.0 openmpi/5.0.7
export PATH=$HOME/Develop/bin:$PATH OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
md5sum $(which xaces2 xjoda xvscf xvtran xintprc xpccd xoorot xrunpccd)
python3 $HOME/Develop/ACES_PROGRAM/aces2/tests/pccd_slow/run_slow.py \
        --jobs $SLURM_CPUS_PER_TASK "$@"
