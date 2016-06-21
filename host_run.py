
import sys
import os
import subprocess

# sanity check first
# ------------------------------------------------------------------------------------------------

# check argument count
if(not len(sys.argv) > 6):
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

matrixList = []
for mtxName in os.listdir(mtxDir):

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
