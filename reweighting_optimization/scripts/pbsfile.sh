#! /bin/bash
#PBS -W group_list=shirtsgroup
#PBS -q shirtscluster
#PBS -l walltime=20:00:00
#PBS -l select=1:mpiprocs=8:ncpus=8
#PBS -m abe
#PBS -M bsz9ur@virginia.edu
#PBS -r n

# submit a single PBS job for a free energy calc
# location of GROMACS binaries on cluster: /h3/n1/shirtsgroup/gromacs_46/install/bin/
# location of GROMACS top files on cluster: /h3/n1/shirtsgroup/gromacs_46/install/share/gromacs/top

module load mpich3-gnu
cat $$PBS_NODEFILE

cd $$PBS_O_WORKDIR

NP=`wc -l < $$PBS_NODEFILE`
NN=`sort -u $$PBS_NODEFILE | wc -l`
echo Number of nodes is $$NN
echo Number of processors is $$NP

##important for some reason
export OMP_NUM_THREADS=1

##$LOCATION/grompp_d -f $NAME.$NUM.mdp -c $NAME.gro -p $NAME.top -o $NAME.$NUM.tpr -maxwarn 2

sleep 1

$LOCATION/mdrun_d -ntmpi 8 -pin off -dlb yes -deffnm $NAME.$NUM -dhdl $NAME.$NUM.dhdl.xvg -o $NAME.$NUM.trr

# print end time
echo
echo "Job Ended at `date`"
echo "###################################################################"
