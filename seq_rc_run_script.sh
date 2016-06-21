#!/bin/bash

bin_parallel=seq_rc;
#bin_parallel=build_generic/seq_rc;
path=input;

export LD_LIBRARY_PATH=/home/mbasaran/lib;
#export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64;

########################################################################
################### CURRENT MATRIX LIST ################################
########################################################################
: '
matrices[0]="atmosmodd";isSymmetric[0]=0;hasValues[0]=1;
matrices[3]="circuit5M_dc"; isSymmetric[3]=0; hasValues[3]=1;
matrices[5]="Freescale1"; isSymmetric[5]=0; hasValues[5]=1;
matrices[11]="rajat31"; isSymmetric[11]=0; hasValues[11]=1;
matrices[12]="webbase-1M"; isSymmetric[12]=0; hasValues[12]=1;
matrices[13]="wheel_601"; isSymmetric[13]=0; hasValues[13]=1;
matrices[8]="memchip"; isSymmetric[8]=0; hasValues[8]=1;
matrices[14]="pre2"; isSymmetric[14]=0; hasValues[14]=1;
matrices[15]="tmt_unsym"; isSymmetric[15]=0; hasValues[15]=1;
matrices[17]="FullChip"; isSymmetric[17]=0; hasValues[17]=1;
matrices[18]="atmosmodl"; isSymmetric[18]=0; hasValues[18]=1;
matrices[21]="cage14"; isSymmetric[18]=0; hasValues[18]=1;
matrices[22]="cage15"; isSymmetric[18]=0; hasValues[18]=1;

matrices[2]="bcsstk35"; isSymmetric[2]=1; hasValues[2]=1;
matrices[4]="crankseg_2"; isSymmetric[4]=1; hasValues[4]=1;
matrices[7]="ldoor"; isSymmetric[7]=1; hasValues[7]=1;
matrices[19]="delaunay_n23"; isSymmetric[19]=1; hasValues[19]=0;

matrices[1]="bcspwr05"; isSymmetric[1]=1; hasValues[1]=0;
matrices[9]="pcrystk03"; isSymmetric[9]=1; hasValues[9]=0;
matrices[10]="pkustk09"; isSymmetric[10]=1; hasValues[10]=0;

matrices[6]="in-2004"; isSymmetric[6]=0; hasValues[6]=0; 
matrices[16]="wikipedia-20061104"; isSymmetric[16]=0; hasValues[16]=0;
matrices[20]="hugebubbles-00010"; isSymmetric[20]=1; hasValues[20]=0;
'
#######################################################################


mtx_name="atmosmodl";
mtx_isSymmetric=0;
mtx_hasValues=1;

echo "=============================== ${mtx_name} ===============================";

row_col=rownet;
./${bin_parallel} matrix=${path}/${mtx_name}/${mtx_name}.mtx isSymmetric=${mtx_isSymmetric} hasValues=${mtx_hasValues};

