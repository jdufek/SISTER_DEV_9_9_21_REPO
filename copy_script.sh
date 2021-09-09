
#!/bin/bash

RUNNAME='LDM_F7_R1_HR'
RUNDIR='/projects/dufeklab/shared/LDM/'$RUNNAME

mkdir $RUNDIR

cp ./COMPRESS* $RUNDIR/
cp ./DENSITY* $RUNDIR/
cp ./INPUT $RUNDIR/
cp ./MELT* $RUNDIR/
cp ./MU* $RUNDIR/
cp ./PRESSURE* $RUNDIR/
cp ./start_time $RUNDIR/
cp ./TEMP $RUNDIR/
cp ./TIMESTEPS $RUNDIR/
cp ./VOLUME* $RUNDIR/
cp ./M_FLUX_AVG $RUNDIR/

cd $RUNDIR
cd ..
chmod -R 775 ./$RUNNAME
