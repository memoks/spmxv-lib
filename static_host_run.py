
import os
import sys
import subprocess

# check argument count
if(not len(sys.argv) > 5):
    print "Incorrect number of arguments!!"
    print "Usage: " + sys.argv[0] + " <binary> <mtx_dir> <cahce_size> <block_count> <threads_per_block> <output_file>"
    sys.exit("Example python " + sys.argv[0] + " /build_generic/hybrid_bipartition /home/endoplasmic/workspace_c/spmxv-lib/input 1024 60 4 results.csv")

# check binary
binary = sys.argv[1]
if(not os.path.isfile(binary)):
    sys.exit(binary + " binary not found!!")

# check matrix directory
mtxDir = sys.argv[2]
if(not os.path.isdir(mtxDir)):
    sys.exit("Invalid matrix directory, see if the path is correct path!!")

cacheSize = sys.argv[3]
numBlocks = sys.argv[4]
threadsPerBlock = sys.argv[5]
outputFile = sys.argv[6]

# set environment variables
subprocess.call(["export", "LD_LIBRARY_PATH=/opt/intel/composer_xe_2015/mkl/lib/intel64:/opt/intel/composer_xe_2015/lib/intel64"], shell=True)
subprocess.call(["export", "MIC_ENV_PREFIX=MIC"], shell=True)
subprocess.call(["export", "MIC_LD_LIBRARY_PATH=/opt/intel/composer_xe_2015/mkl/lib/mic:/opt/intel/composer_xe_2015/lib/mic"], shell=True)
subprocess.call(["export", "MIC_KMP_PLACE_THREADS=" + numBlocks + "c," + threadsPerBlock + "t"], shell=True)
subprocess.call(["export", "MIC_OMP_NUM_THREADS=" + str(int(numBlocks) * int(threadsPerBlock))], shell=True)

print mtxDir

mtxListStr = "Na5 c-60 c-57 circuit_4 c-61 graham1 a2nnsnsl ckt11752_dc_1 ckt11752_tr_0 a0nsdsil vsp_south31_slptsk bcsstk38 vsp_model1_crew1_cr42_south31 t0331-4l fp email-Enron ship_001 c-65 c-54 bundle1 c-56 ca-AstroPh c-55 blockqp1 rajat15 c-63 rajat18 c-59 TSOPF_FS_b162_c1 c-66 c-66b Si10H16 vsp_c-60_data_cti_cs4 c-67 c-67b soc-Epinions1 c-64 c-64b pdb1HYS nemeth20 c-58 c-62 c-62ghs Zd_Jac2 rajat23 c-68 TSOPF_FS_b39_c7 rajat20 rajat25 rajat28 nsct rajat16 rajat17 dc1 dc2 dc3 trans4 trans5 c-69 Si5H12 mip1 cbuckle c-70 boyd1 pkustk02 cyl6 c-72 TSC_OPF_300 av41092 sme3Da ASIC_100k vsp_mod2_pgp2_slptsk case39 soc-sign-epinions c-71 soc-Slashdot0811 HTC_336_9129 HTC_336_4438 transient soc-Slashdot0902 gyro_k gyro nemsemm1 vsp_msc10848_300sep_100in_1Kout torso1 kron_g500-logn16 3D_51448_3D ibm_matrix_2 SiO dbir1 Maragal_8 gupta1 fem_filter dbir2 c-73 c-73b Maragal_7 12month1 ns3Da Raj1 raefsky4 trdheim net100 boyd2 TSOPF_FS_b162_c3 web-NotreDame lp1 TSOPF_FS_b39_c19 kron_g500-logn17 net125 invextr1_new sme3Db Ga3As3H12 TSOPF_FS_b162_c4 rajat21 rajat24 bcsstk30 GaAsH6 cnr-2000 net150 matrix_9 ins2 gupta2 ASIC_320k net4-1 RM07R vanbody citationCiteseer c-big TSOPF_FS_b39_c30 nasasrb sme3Dc matrix-new_3 srb1 Si34H36 kron_g500-logn18 Ge99H100 3dtube Ge87H76 pkustk04 ASIC_680k webbase-1M Ga10As10H30 barrier2-1 barrier2-2 barrier2-3 barrier2-4 barrier2-10 barrier2-11 barrier2-12 barrier2-9 vsp_bcsstk30_500sep_10in_1Kout rajat29 Ga19As19H42 F2 m_t1 mono_500Hz SiO2 rajat30 para-4 wiki-Talk para-10 para-5 para-6 para-7 para-8 para-9 Si41Ge41H72 x104 pkustk13 pkustk12 bmw7st_1 kron_g500-logn19 flickr CO pkustk14 Ga41As41H72 bmwcra_1 Si87H76 ohne2 eu-2005 bmw3_2 in-2004 troll wikipedia-20051105 coPapersCiteseer as-Skitter F1 FullChip coPapersDBLP dielFilterV3clx wikipedia-20060925 Fault_639 wikipedia-20061104 inline_1 wikipedia-20070206 tp-6 degme Stanford_Berkeley stat96v3 neos3 Stanford Chebyshev4 karted EternityII_E language 2D_54019_highK ts-palko rail4284 3D_28984_Tetra rail2586 bas1lp stat96v4 Zd_Jac6 Zd_Jac3 Maragal_6 std1_Jac2 std1_Jac3 connectus cont1_l web-BerkStan web-Google IMDB amazon0601 amazon0505 LargeRegFile amazon0312 sls web-Stanford neos IG5-18 ESOC NotreDame_www TSOPF_RS_b39_c30 IG5-17 neos1 neos2 mri2 TSOPF_RS_b162_c4 soc-sign-Slashdot090216 soc-sign-Slashdot090221 soc-sign-Slashdot081106 IG5-16 TSOPF_RS_b39_c19 cit-HepPh TSOPF_RS_b162_c3 email-EuAll HEP-th-new cit-HepTh HEP-th"
matrixList = mtxListStr.split(" ")
for mtxName in matrixList:
    
    if(not os.path.isdir(mtxDir + "/" + mtxName)):
       continue

    mtxPath = mtxDir + "/" + mtxName + "/" + mtxName + ".mtx"
    
    if(not os.path.isfile(mtxPath)):
       continue

    print mtxPath

    args = []
    args.append(binary)
    args.append("matrix=" + mtxPath)
    args.append("num_blocks=" + numBlocks)
    args.append("threads_per_block=" + threadsPerBlock)
    args.append("storage_format=HYBRID_JDS_CSR")
    args.append("partition_type=1D_row_wise")
    args.append("partition_method=recursive_bipartition")
    args.append("ordering_type=column_net")
    args.append("targeted_cache_size_kb=128")
    args.append("simd_length=8")
    args.append("run_count_warm_up=10")
    args.append("run_count_measure=100")

    # Do not disturb other coprocessors
    subprocess.call(["export", "OFFLOAD_INIT=on_offload"], shell=True)
    run_binary = subprocess.Popen(args, stdout=subprocess.PIPE)
    run_status = run_binary.wait()
    run_output = run_binary.stdout.read()

    # Append results to general output file
    if(run_status == 0):
       with open(outputFile, "a") as file: file.write(run_output)
    else:
       with open(outputFile, "a") as file: file.write(mtxPath + " error occurred!!\n")
