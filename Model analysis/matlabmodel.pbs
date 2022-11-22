#!/bin/bash
#PBS -r n
#PBS -q default
#PBS -N testyuqi
#PBS -l nodes=1:ppn=16
#PBS -l walltime=800:00:00
# Load the MATLAB module
# module load matlab_R2008b
# module load matlab
#
# cd to the directory where the job was submitted
cd $PBS_O_WORKDIR
echo Current directory is $PBS_O_WORKDIR
echo This job run on the following processors:
echo `cat $PBS_NODEFILE` 
#
# Run MATLAB
/u1/local/MATLAB/R2015a/bin/matlab -nodisplay -nosplash -nojvm < /home/kiki/kikifolder/model/surrogateFFTmodel23.m

