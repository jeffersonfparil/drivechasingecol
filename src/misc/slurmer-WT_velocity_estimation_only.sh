#!/bin/bash

### Input parameter description
# (1) ${HOME_DIR}          directory where qunatinemo2 folder, plink, emmax, gemma, gcta, and genomic_prediction folder are located
# (2) ${repN}              replicate number
# (3) SPARTAN_ACCOUNT      Spartan account name
# (4) SPARTAN_PARTITION    Sparton node partition
# (5) SPARTAN_CORES        Partition-specific number of cores
# (6) SPARTAN_RAM          Partition-specific RAM

### Hardcoded input parameters
DIR=/data/cephfs/punim1173/driveChase/
SRC_DIR=${DIR}drivechasingecol/src/misc
SPARTAN_ACCOUNT=punim1173
SPARTAN_PARTITION=snowy
SPARTAN_TIMELIMIT="0-1:0:00"
SPARTAN_CORES=16
SPARTAN_RAM=32

### Iterate across all 131,220 combinations of the simulation variables plus the replications 
###     - iterating because job arrays is clunky, i.e. submitting all my jobs to physical instead of snowy where I specified
###     - also in-script parallelisation with foreach and doParallel is not working
for rep in $(seq 1 100)
do
### build the slurm script
echo -e '#!/bin/bash' > driveChase_rep-${rep}.slurm
echo -e "
# Partition for the job:
#SBATCH --partition=${SPARTAN_PARTITION}
# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1
# The name of the job:
#SBATCH --job-name=driveChase_${rep}REP
# The project ID which this job should run under:
#SBATCH --account=${SPARTAN_ACCOUNT}
# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SPARTAN_CORES}
# The amount of memory in megabytes per process in the job:
#SBATCH --mem=${SPARTAN_RAM}GB
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=${SPARTAN_TIMELIMIT}
# Send yourself an email when the job:
# aborts abnormally (fails)
# #SBATCH --mail-type=FAIL
# begins
# #SBATCH --mail-type=BEGIN
# ends successfully
# #SBATCH --mail-type=END
# Use this email address:
#SBATCH --mail-user=parilj@unimelb.edu.au

### load modules
module load r/4.0.0
module load parallel/20170422

### Execute
rep=${rep}
time \\
parallel --jobs ${SPARTAN_CORES} \\
Rscript ${SRC_DIR}/driveChaseEcolRun-WT_velocity_estimation_only.R \\
    ${DIR} \\" >> driveChase_rep-${rep}.slurm
echo '{} ::: $(seq $(echo ${rep}*260 - 260 + 1 | bc) $(echo ${rep}*260 | bc))' >> driveChase_rep-${rep}.slurm
### submit into queue
sbatch driveChase_rep-${rep}.slurm
done

##################
### MONITORING ###
##################
# ### add to ~/.bashrc
# DIR=/data/cephfs/punim1173/driveChase/
# SRC_DIR=${DIR}drivechasingecol/src/
# alias Q='squeue -u parilj'
# alias ERROR='grep "error" ${SRC_DIR}slurm-*.out'
# alias DISK='check_project_usage'
# alias PROGRESS='FINISHED1=$(ls ${DIR}driveChaseEcol-OUTPUT/ | grep "\.csv$" | wc -l); FINISHED2=$(ls ${DIR}driveChaseEcol-OUTPUT/ONEROW_csv_from_spartan/ | grep "\.csv$" | wc -l); PERC=$(echo "scale=2; ( ${FINISHED1} + ${FINISHED2} ) * 100 / 26000" | bc); echo "${PERC}% finished; (${FINISHED}/26000)"'

### monitoring on nectar VMs
# DIR=/homevol/speel
# SRC_DIR=${DIR}/drivechasingecol/src
# # rep=1
# # time nohup parallel -j 15 Rscript ${SRC_DIR}/driveChaseEcolRun.R ${DIR}/ {} ::: $(seq $(echo ${rep}*260 - 260 + 1 | bc) $(echo ${rep}*260 | bc)) &
# nRun=$(grep "^tar" ${DIR}/nohup.out | wc -l); nFin=$(ls ${DIR}/driveChaseEcol-OUTPUT/ | grep "\.csv$" | wc -l); echo RUN:${nRun}-FINISHED:${nFin}
