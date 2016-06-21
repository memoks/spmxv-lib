
bin_parallel=binary/fast_colNet_hybrid_JDS_CSR;
path=input;

export LD_LIBRARY_PATH=/home/mbasaran/lib;
#export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64;

core_count=60;
thread_count=4;
partition_type="colnet";
cache_size_kb=64;

# inter compiler thread assignment variables
KMP_PLACE_THREADS=${core_count}c,${thread_count}t;
OMP_NUM_THREADS=$((${core_count}*${thread_count}));

export KMP_PLACE_THREADS=${KMP_PLACE_THREADS};
export OMP_NUM_THREADS=${OMP_NUM_THREADS};

mtx_name=$1;
./${bin_parallel} matrix=${path}/${mtx_name}/${mtx_name}.mtx numBlocks=${core_count} threadsPerBlock=${thread_count} partitionType=${partition_type} cacheSizeKB=${cache_size_kb};

