#!/bin/bash

# run_rmd.sh
# Run rmarkdown script in batch mode
# Alexey Larionov, 15Mar2016
# Use: sbatch run.rmd rmd_script

# ---------------------------------------- #
#           sbatch instructions            #
# ---------------------------------------- #

#SBATCH -J run_rmd
#SBATCH --time=12:00:00
#SBATCH -A TISCHKOWITZ-SL2
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH -p sandybridge
##SBATCH --qos=INTR

# Modules section (required, do not remove)
. /etc/profile.d/modules.sh
module purge
module load default-impi
module load gcc/5.2.0
module load boost/1.50.0
module load texlive/2015
module load pandoc/1.15.2.1

# Set initial working folder
cd "${SLURM_SUBMIT_DIR}"

# Report settings
echo "Job name: ${SLURM_JOB_NAME}"
echo "Allocated node: $(hostname)"
echo "Initial working folder:"
echo "${SLURM_SUBMIT_DIR}"
echo ""
echo " ------------------ Output ------------------ "
echo ""

# Files and folders
root_folder="/scratch/medgen/scripts/rscripts_05.16"
rmd_file=${1}
rmd_script="${root_folder}/scripts/${rmd_file}"
html_report="${root_folder}/scripts/${rmd_file%.Rmd}.html"

# Check that Rmd script exists
if [ ! -f "${rmd_script}" ] 
then
  echo "Use: sbatch run.rmd rmd_script"
  echo "Rmd script does not exist"
  echo "Script terminated"
  exit
fi

# Make list of R expressions (avoid blank spaces)
r_expressions="
 -e library('rmarkdown',lib='/scratch/medgen/tools/r/R-3.2.2/lib64/R/library/')
 -e render('"${rmd_script}"',output_file='"${html_report}"')"

# Run the list of R expressions
/scratch/medgen/tools/r/R-3.2.2/bin/R "${r_expressions}"

# Completion message
echo "Script completed"