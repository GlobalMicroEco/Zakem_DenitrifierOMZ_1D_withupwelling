#!/bin/bash
## now loop through the above array

#SBATCH --time=10:10:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10G   # memory per CPU core
#SBATCH -J "SimMixing1"   # job name
#SBATCH --mail-user=liangxu@caltech.edu   # email address


## /SBATCH -p general # partition (queue)
## /SBATCH -o slurm.%N.%j.out # STDOUT
## /SBATCH -e slurm.%N.%j.err # STDERR

module load julia/1.8.1
julia run_addW.jl



# for mbid in {1..19}
# do

# echo "#!/bin/bash" > exp$mbid
# echo "#SBATCH --time=1-00:00:00" >> exp$mbid
# echo "#SBATCH --ntasks=1" >> exp$mbid
# echo "#SBATCH --nodes=1" >> exp$mbid
# echo "#SBATCH --mem-per-cpu=10G" >> exp$mbid
# echo "#SBATCH -J 'expnb'$mbid " >> exp$mbid
# echo "#SBATCH --mail-user=liangxu@caltech.edu" >> exp$mbid

# echo "module load julia/1.8.1" >> exp$mbid

# echo "julia Cluster_sim_change_s.jl $mbid " >> exp$mbid
# echo "rm exp$mbid" >> exp$mbid

# sbatch exp$mbid


# done
