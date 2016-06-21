#!/bin/bash

export LD_LIBRARY_PATH=/home/mbasaran/lib;
#export LD_LIBRARY_PATH=/mnt/compgen/inhouse/icc/maggie/lib/intel64;


matrices[0]="atmosmodd"; isSymmetric[0]="0"; hasValues[0]="1";						# symmetric
matrices[1]="atmosmodl"; isSymmetric[1]="0"; hasValues[1]="1";						# symmetric
matrices[2]="circuit5M_dc"; isSymmetric[2]="0"; hasValues[2]="1";					# unsymmetric %91
matrices[3]="Freescale1"; isSymmetric[3]="0"; hasValues[3]="1";						# unsymmetric %91
matrices[4]="rajat31"; isSymmetric[4]="0"; hasValues[4]="1";						# symmetric
matrices[5]="webbase-1M"; isSymmetric[5]="0"; hasValues[5]="1";						# unsymmetric %10
matrices[6]="wheel_601"; isSymmetric[6]="0"; hasValues[6]="0";						# unsymmetric %0
matrices[7]="memchip"; isSymmetric[7]="0"; hasValues[7]="1";						# unsymmetric %91
matrices[8]="pre2"; isSymmetric[8]="0"; hasValues[8]="1";							# unsymmetric %36
matrices[9]="tmt_unsym"; isSymmetric[9]="0"; hasValues[9]="1";						# symmetric
matrices[10]="FullChip"; isSymmetric[10]="0"; hasValues[10]="1";					# symmetric
matrices[11]="cage14"; isSymmetric[11]="0"; hasValues[11]="1";						# symmetric
matrices[12]="cage15"; isSymmetric[12]="0"; hasValues[12]="1";						# symmetric
matrices[13]="bcsstk35"; isSymmetric[13]="1"; hasValues[13]="1";
matrices[14]="crankseg_2"; isSymmetric[14]="1"; hasValues[14]="1";					# symmetric
matrices[15]="ldoor"; isSymmetric[15]="1"; hasValues[15]="1";						# symmetric
matrices[16]="delaunay_n23"; isSymmetric[16]="1"; hasValues[16]="0";				# symmetric
matrices[17]="pkustk09"; isSymmetric[17]="1"; hasValues[17]="0";
matrices[18]="in-2004"; isSymmetric[18]="0"; hasValues[18]="0";						# unsymmetric %36
matrices[19]="wikipedia-20061104"; isSymmetric[19]="0"; hasValues[19]="0";			# unsymmetric %12
matrices[20]="hugebubbles-00010"; isSymmetric[20]="1"; hasValues[20]="0";			# symmetric # too big
matrices[21]="rajat31"; isSymmetric[21]="0"; hasValues[21]="1";						# symmetric
matrices[22]="web-Google"; isSymmetric[22]="0"; hasValues[22]="0";					# unsymmetric %31
matrices[23]="hugetrace-00000"; isSymmetric[23]="1"; hasValues[23]="0";				# symmetric # too big
matrices[24]="channel-500x100x100-b050"; isSymmetric[24]="1"; hasValues[24]="0";	# symmetric
matrices[25]="StocF-1465"; isSymmetric[25]="1"; hasValues[25]="1";					# symmetric
matrices[26]="hugetric-00020"; isSymmetric[26]="1"; hasValues[26]="0";				# symmetric
matrices[27]="as-Skitter"; isSymmetric[27]="1"; hasValues[27]="0";					# symmetric
matrices[28]="rel9"; isSymmetric[28]="0"; hasValues[28]="1";						# colPar - no transpose
matrices[29]="germany_osm"; isSymmetric[29]="1"; hasValues[29]="0";					# symmetric
matrices[30]="hugetric-00010"; isSymmetric[30]="1"; hasValues[30]="0";				# symmetric
matrices[31]="kkt_power"; isSymmetric[31]="1"; hasValues[31]="1";					# symmetric
matrices[32]="atmosmodj"; isSymmetric[32]="0"; hasValues[32]="1";					# symmetric
matrices[33]="thermal2"; isSymmetric[33]="1"; hasValues[33]="1";					# symmetric
matrices[34]="G3_circuit"; isSymmetric[34]="1"; hasValues[34]="1";					# symmetric
matrices[35]="Rucci1"; isSymmetric[35]="0"; hasValues[35]="1"; 						# unsymmetric %0
matrices[36]="ljournal-2008"; isSymmetric[36]="0"; hasValues[36]="0";				# unsymmetric %73
matrices[37]="soc-LiveJournal1"; isSymmetric[37]="0"; hasValues[37]="0";			# unsymmetric %75
matrices[38]="circuit5M"; isSymmetric[38]="0"; hasValues[38]="1";					# symmetric
matrices[39]="wb-edu"; isSymmetric[39]="0"; hasValues[39]="0";						# unsymmetric %33
matrices[40]="spal_004"; isSymmetric[40]="0"; hasValues[40]="1";					# unsymmetric %0
matrices[41]="wikipedia-20070206"; isSymmetric[41]="0"; hasValues[41]="0";			# unsymmetric %12
matrices[42]="patents"; isSymmetric[42]="0"; hasValues[42]="0";						# unsymmetric %0
matrices[43]="RM07R"; isSymmetric[43]="0"; hasValues[43]="1";						# unsymmetric %93
matrices[44]="mouse_gene"; isSymmetric[44]="1"; hasValues[44]="1";					# symmetric
matrices[45]="nd24k"; isSymmetric[45]="1"; hasValues[45]="1";						# symmetric
matrices[46]="F1"; isSymmetric[46]="1"; hasValues[46]="1";							# symmetric
matrices[47]="human_gene1"; isSymmetric[47]="1"; hasValues[47]="1";					# symmetric
matrices[48]="12month1"; isSymmetric[48]="0"; hasValues[48]="0";					# unsymmetric %0
matrices[49]="kron_g500-logn18"; isSymmetric[49]="1"; hasValues[49]="1";			# symmetric
matrices[50]="eu-2005"; isSymmetric[50]="0"; hasValues[50]="0";						# unsymmetric %28
matrices[51]="Ga41As41H72"; isSymmetric[51]="1"; hasValues[51]="1";					# symmetric
matrices[52]="human_gene2"; isSymmetric[52]="1"; hasValues[52]="1";					# symmetric
matrices[53]="cit-Patents"; isSymmetric[53]="0"; hasValues[53]="0";					# unsymmetric %0
matrices[54]="Hamrle3"; isSymmetric[54]="0"; hasValues[54]="1";						# unsymmetric %0
matrices[55]="Stanford_Berkeley"; isSymmetric[55]="0"; hasValues[55]="0";			# unsymmetric %25
# matrices[56]=""; isSymmetric[56]=""; hasValues[56]="";
# matrices[57]=""; isSymmetric[57]=""; hasValues[57]="";
# matrices[58]=""; isSymmetric[58]=""; hasValues[58]="";
# matrices[59]=""; isSymmetric[59]=""; hasValues[59]="";





input_path=input/ozan_mtx;
bin_name=build_generic/mtx_transposer;

for i in {0..55}
do
	matrix=${matrices[${i}]};
	isSymmetric=${isSymmetric[${i}]};
	hasValues=${hasValues[${i}]};
		
	echo ""
	echo "=============================== ${matrix} ===============================";
	
	./${bin_name} matrix=${input_path}/${matrix}/${matrix}.mtx orderings=${input_path}/${matrix}/rownet/${matrix}.mtx.order.R.1.1 blocks=${input_path}/${matrix}/rownet/${matrix}.mtx.outdim.R.1.1 isSymmetric=${isSymmetric} hasValues=${hasValues};
done