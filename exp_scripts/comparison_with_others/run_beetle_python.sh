#!/bin/bash
 
#BSUB -P csc143
#BSUB -W 01:59
#BSUB -nnodes 1

#BSUB -J run_beetle_python
#BSUB -o run_beetle_python.%J.out
#BSUB -e run_beetle_python.%J.err 

CURRDIR=$(pwd)
LOGDIRNAME=run_beetle_python_log

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cp $CURRDIR/beetle*.py $LOGDIRNAME

cd $LOGDIRNAME

DATANAME=stagbeetle832x832x494.dat
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME

cp $DATASETPATH .

# one process cases
# using default options for all
jsrun python3 ./beetle_uniform.py &> beetle_uniform.log

jsrun python3 ./beetle_gaussian.py &> beetle_gaussian.log

# copy things back
cp -r *.log $CURRDIR/$LOGDIRNAME

# clean the run dir
#cd ..

#rm -r $LOGDIRNAME
