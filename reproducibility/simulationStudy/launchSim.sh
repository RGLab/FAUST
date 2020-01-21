#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --partition=[PUT YOUR PARTITION NAME HERE]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=[PUT YOUR TIME LIMITS HERE]
#SBATCH --mem=60000
#SBATCH -o ./logs/sim$(($1))$(($2))$(($3))$(($4))
#SBATCH -J sim$(($1))$(($2))$(($3))$(($4))
#SBATCH --threads-per-core=1
#SBATCH --dependency=singleton

echo "Start of program at `date`" 
Rscript --no-save --no-restore ./01_runSimulation.R -a $1 -b $2 -c $3 -d $4
echo "End of program at `date`"
EOT
