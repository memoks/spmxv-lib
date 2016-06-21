# BUILD_ARCH = generic|mic_offload|mic_native 
# Must be specified from command line to build for particular architechture
# PROGRAM_MODE = 0|1|2|3 (no_output|measure|debug|trace)
# When PROGRAM_MODE is not specified default option is 1:measure
################################################################################################

# Default program mode is "measure-mode" (minimum output mode)
ifeq ($(PROGRAM_MODE), )
PROGRAM_MODE = 1
endif

# SELF_DIR = $(dir $(lastword $(MAKEFILE_LIST)))
PROJECT_DIR = $(shell pwd)

# all other builds will/can use generic routines, so define the path regardless
# in addition to generic routines, some builds can have extra routines defined  
# under make_macros directory in that specific architecture's macro file.  
GENERIC_PATH = arch/generic

# src directory tree
SRC = src
DS = data_structure
SC = scheduler
TM = timer
IN = include
TS = test
INT = interface
MA = main
IO = io
ALG = algorithm
TD = task_decomposition
UT = util
CU = control_unit

# Choose a run configuration with argument RUN_MODE
# Usage: RUN_MODE=<NAME_OF_THE_FILE_INCLUDING_MAIN_METHOD>
#################################################################################################

MAIN_FILE = $(RUN_MODE)

# For architecture specific settings (includes additional codes and compiler options)
#################################################################################################

include make_macros/$(BUILD_ARCH)_macros

#################################################################################################

FILES = $(MA)/$(MAIN_FILE) kernel \
	spmxv_sequential spmxv_omp_loop spmxv_omp_task spmxv_gws spmxv_dws spmxv_static \
	$(IO)/cli $(IO)/input_parser $(IO)/converter $(IO)/mmio \
	$(SC)/job_batch $(SC)/job_queue $(SC)/block_info $(SC)/ring_scheduler $(SC)/ring_sched $(SC)/tree_scheduler $(SC)/tree_sched \
	$(TM)/custom_timer \
	$(DS)/sub_mtx $(DS)/comm_tree $(DS)/spm_storage $(DS)/vector $(DS)/tree $(DS)/list_generic $(DS)/quintet $(DS)/stack_generic \
	$(ALG)/algorithm \
	$(TD)/partitioning \
	$(UT)/utility \
	$(CU)/dws $(CU)/gws $(CU)/static $(CU)/omp_loop $(CU)/fast_run

OBJS = $(foreach file, $(FILES), $(BUILD_DIR)/$(file).o) $(STATIC_LIBS)

HEADERS = $(SRC)/include/config.h \
	$(SRC)/include/$(DS)/list.h \
	$(SRC)/include/control_unit/cu_options.h

all: clean_up $(BUILD_DIR)/$(RUN_MODE)
	@echo ""
	@echo "PROGRAM_MODE = (0:no_output | 1:measure | 2:debug | 3:trace)"
	@echo "Project directory path: $(PROJECT_DIR)"
	@echo "Program mode of current build is $(PROGRAM_MODE)"
	@echo "Run mode: $(RUN_MODE)"
	@echo "Build successful!!"
	@echo ""

$(BUILD_DIR)/$(RUN_MODE) : $(OBJS) $(ARCH_OBJS) $(HEADERS)
	$(LD) $(OBJS) $(ARCH_OBJS) $(LDFLAGS) -o $(BUILD_DIR)/$(RUN_MODE)

clean_up :
	rm -rf $(BUILD_DIR)
	mkdir $(BUILD_DIR)
	mkdir $(BUILD_DIR)/$(IO)
	mkdir $(BUILD_DIR)/$(DS)
	mkdir $(BUILD_DIR)/$(SC)
	mkdir $(BUILD_DIR)/$(TM)
	mkdir $(BUILD_DIR)/$(TS)
	mkdir $(BUILD_DIR)/$(INT) 
	mkdir $(BUILD_DIR)/$(MA)
	mkdir $(BUILD_DIR)/$(ALG)
	mkdir $(BUILD_DIR)/$(TD)
	mkdir $(BUILD_DIR)/$(UT)
	mkdir $(BUILD_DIR)/$(CU)
	@echo "Build directory created!"
	@echo ""
			

clean : 
	rm -rf $(ALL_OBJS) $(BUILD_DIR)/$(MAIN_FILE).o $(BUILD_DIR)/$(RUN_MODE)
	rm -rf $(BUILD_DIR)/$(DS) $(BUILD_DIR)/$(SC) $(BUILD_DIR)/$(TM) $(BUILD_DIR)/$(TS) $(BUILD_DIR)/$(INT) $(BUILD_DIR)/$(MA)

# Main File
#################################################################################################

$(BUILD_DIR)/$(MA)/$(MAIN_FILE).o : $(SRC)/$(MA)/$(MAIN_FILE).c
	$(CC) $(CFLAGS) $(SRC)/$(MA)/$(MAIN_FILE).c -o $(BUILD_DIR)/$(MA)/$(MAIN_FILE).o

# Algorithm
#################################################################################################

$(BUILD_DIR)/$(ALG)/algorithm.o : $(SRC)/$(ALG)/algorithm.c $(SRC)/$(IN)/$(ALG)/algorithm.h
	$(CC) $(CFLAGS) $(SRC)/$(ALG)/algorithm.c -o $(BUILD_DIR)/$(ALG)/algorithm.o

# Task Decomposition
#################################################################################################

$(BUILD_DIR)/$(TD)/partitioning.o : $(SRC)/$(TD)/partitioning.c $(SRC)/$(IN)/$(TD)/partitioning.h
	$(CC) $(CFLAGS) $(SRC)/$(TD)/partitioning.c -o $(BUILD_DIR)/$(TD)/partitioning.o

# Conjugate Gradined Interface (Spmxv routines defined here for each implementation)
#################################################################################################

$(BUILD_DIR)/$(INT)/common.o : $(SRC)/$(INT)/common.c $(SRC)/$(IN)/$(INT)/common.h
	$(CC) $(CFLAGS) $(SRC)/$(INT)/common.c -o $(BUILD_DIR)/$(INT)/common.o

$(BUILD_DIR)/$(INT)/dws.o : $(SRC)/$(INT)/dws.c $(SRC)/$(IN)/$(INT)/dws.h
	$(CC) $(CFLAGS) $(SRC)/$(INT)/dws.c -o $(BUILD_DIR)/$(INT)/dws.o

$(BUILD_DIR)/$(INT)/sequential.o : $(SRC)/$(INT)/sequential.c $(SRC)/$(IN)/$(INT)/sequential.h
	$(CC) $(CFLAGS) $(SRC)/$(INT)/sequential.c -o $(BUILD_DIR)/$(INT)/sequential.o

# IO
#################################################################################################

$(BUILD_DIR)/$(IO)/cli.o : $(SRC)/$(IO)/cli.c $(SRC)/$(IN)/$(IO)/cli.h
	$(CC) $(CFLAGS) $(SRC)/$(IO)/cli.c -o $(BUILD_DIR)/$(IO)/cli.o

$(BUILD_DIR)/$(IO)/input_parser.o : $(SRC)/$(IO)/input_parser.c $(SRC)/$(IN)/$(IO)/input_parser.h
	$(CC) $(CFLAGS) $(SRC)/$(IO)/input_parser.c -o $(BUILD_DIR)/$(IO)/input_parser.o 
	
$(BUILD_DIR)/$(IO)/converter.o : $(SRC)/$(IO)/converter.c $(SRC)/$(IN)/$(IO)/converter.h
	$(CC) $(CFLAGS) $(SRC)/$(IO)/converter.c -o $(BUILD_DIR)/$(IO)/converter.o

$(BUILD_DIR)/$(IO)/mmio.o : $(SRC)/$(IO)/mmio.c $(SRC)/$(IN)/$(IO)/mmio.h
	$(CC) $(CFLAGS) $(SRC)/$(IO)/mmio.c -o $(BUILD_DIR)/$(IO)/mmio.o

# Util
#################################################################################################

$(BUILD_DIR)/$(UT)/utility.o : $(SRC)/$(UT)/utility.c $(SRC)/$(IN)/$(UT)/utility.h
	$(CC) $(CFLAGS) $(SRC)/$(UT)/utility.c -o $(BUILD_DIR)/$(UT)/utility.o

# Scheduler
#################################################################################################

$(BUILD_DIR)/$(SC)/job_batch.o : $(SRC)/$(SC)/job_batch.c $(SRC)/$(IN)/$(SC)/job_batch.h
	$(CC) $(CFLAGS) $(SRC)/$(SC)/job_batch.c -o $(BUILD_DIR)/$(SC)/job_batch.o

$(BUILD_DIR)/$(SC)/job_queue.o : $(SRC)/$(SC)/job_queue.c $(SRC)/$(IN)/$(SC)/job_queue.h
	$(CC) $(CFLAGS) $(SRC)/$(SC)/job_queue.c -o $(BUILD_DIR)/$(SC)/job_queue.o

$(BUILD_DIR)/$(SC)/block_info.o : $(SRC)/$(SC)/block_info.c $(SRC)/$(IN)/$(SC)/block_info.h
	$(CC) $(CFLAGS) $(SRC)/$(SC)/block_info.c -o $(BUILD_DIR)/$(SC)/block_info.o

$(BUILD_DIR)/$(SC)/ring_sched.o : $(SRC)/$(SC)/ring_sched.c $(SRC)/$(IN)/$(SC)/ring_sched.h
	$(CC) $(CFLAGS) $(SRC)/$(SC)/ring_sched.c -o $(BUILD_DIR)/$(SC)/ring_sched.o

$(BUILD_DIR)/$(SC)/ring_scheduler.o : $(SRC)/$(SC)/ring_scheduler.c $(SRC)/$(IN)/$(SC)/ring_scheduler.h
	$(CC) $(CFLAGS) $(SRC)/$(SC)/ring_scheduler.c -o $(BUILD_DIR)/$(SC)/ring_scheduler.o

$(BUILD_DIR)/$(SC)/tree_sched.o : $(SRC)/$(SC)/tree_sched.c $(SRC)/$(IN)/$(SC)/tree_sched.h
	$(CC) $(CFLAGS) $(SRC)/$(SC)/tree_sched.c -o $(BUILD_DIR)/$(SC)/tree_sched.o

$(BUILD_DIR)/$(SC)/tree_scheduler.o : $(SRC)/$(SC)/tree_scheduler.c $(SRC)/$(IN)/$(SC)/tree_scheduler.h
	$(CC) $(CFLAGS) $(SRC)/$(SC)/tree_scheduler.c -o $(BUILD_DIR)/$(SC)/tree_scheduler.o

# Timer
#################################################################################################

$(BUILD_DIR)/$(TM)/custom_timer.o : $(SRC)/$(TM)/custom_timer.c $(SRC)/$(IN)/$(TM)/custom_timer.h
	$(CC) $(CFLAGS) $(SRC)/$(TM)/custom_timer.c -o $(BUILD_DIR)/$(TM)/custom_timer.o

# Data Structures
#################################################################################################

$(BUILD_DIR)/$(DS)/sub_mtx.o : $(SRC)/$(DS)/sub_mtx.c $(SRC)/$(IN)/$(DS)/sub_mtx.h
	$(CC) $(CFLAGS) $(SRC)/$(DS)/sub_mtx.c -o $(BUILD_DIR)/$(DS)/sub_mtx.o
	
$(BUILD_DIR)/$(DS)/comm_tree.o : $(SRC)/$(DS)/comm_tree.c $(SRC)/$(IN)/$(DS)/comm_tree.h
	$(CC) $(CFLAGS) $(SRC)/$(DS)/comm_tree.c -o $(BUILD_DIR)/$(DS)/comm_tree.o

$(BUILD_DIR)/$(DS)/tree.o : $(SRC)/$(DS)/tree.c $(SRC)/$(IN)/$(DS)/tree.h
	$(CC) $(CFLAGS) $(SRC)/$(DS)/tree.c -o $(BUILD_DIR)/$(DS)/tree.o

$(BUILD_DIR)/$(DS)/spm_storage.o : $(SRC)/$(DS)/spm_storage.c $(SRC)/$(IN)/$(DS)/spm_storage.h
	$(CC) $(CFLAGS) $(SRC)/$(DS)/spm_storage.c -o $(BUILD_DIR)/$(DS)/spm_storage.o

$(BUILD_DIR)/$(DS)/vector.o : $(SRC)/$(DS)/vector.c $(SRC)/$(IN)/$(DS)/vector.h
	$(CC) $(CFLAGS) $(SRC)/$(DS)/vector.c -o $(BUILD_DIR)/$(DS)/vector.o
	
$(BUILD_DIR)/$(DS)/list_generic.o : $(SRC)/$(DS)/list_generic.c $(SRC)/$(IN)/$(DS)/list_generic.h
	$(CC) $(CFLAGS) $(SRC)/$(DS)/list_generic.c -o $(BUILD_DIR)/$(DS)/list_generic.o

$(BUILD_DIR)/$(DS)/quintet.o : $(SRC)/$(DS)/stack_generic.c $(SRC)/$(IN)/$(DS)/quintet.h
	$(CC) $(CFLAGS) $(SRC)/$(DS)/quintet.c -o $(BUILD_DIR)/$(DS)/quintet.o
	
$(BUILD_DIR)/$(DS)/stack_generic.o : $(SRC)/$(DS)/stack_generic.c $(SRC)/$(IN)/$(DS)/stack_generic.h
	$(CC) $(CFLAGS) $(SRC)/$(DS)/stack_generic.c -o $(BUILD_DIR)/$(DS)/stack_generic.o

# SpMxV Control Units / Interfaces
#################################################################################################

$(BUILD_DIR)/$(CU)/fast_run.o : $(SRC)/$(CU)/fast_run.c $(SRC)/$(IN)/$(CU)/fast_run.h
	$(CC) $(CFLAGS) $(SRC)/$(CU)/fast_run.c -o $(BUILD_DIR)/$(CU)/fast_run.o

$(BUILD_DIR)/$(CU)/dws.o : $(SRC)/$(CU)/dws.c $(SRC)/$(IN)/$(CU)/dws.h
	$(CC) $(CFLAGS) $(SRC)/$(CU)/dws.c -o $(BUILD_DIR)/$(CU)/dws.o
	
$(BUILD_DIR)/$(CU)/gws.o : $(SRC)/$(CU)/gws.c $(SRC)/$(IN)/$(CU)/gws.h
	$(CC) $(CFLAGS) $(SRC)/$(CU)/gws.c -o $(BUILD_DIR)/$(CU)/gws.o
	
$(BUILD_DIR)/$(CU)/static.o : $(SRC)/$(CU)/static.c $(SRC)/$(IN)/$(CU)/static.h
	$(CC) $(CFLAGS) $(SRC)/$(CU)/static.c -o $(BUILD_DIR)/$(CU)/static.o

$(BUILD_DIR)/$(CU)/omp_loop.o : $(SRC)/$(CU)/omp_loop.c $(SRC)/$(IN)/$(CU)/omp_loop.h
	$(CC) $(CFLAGS) $(SRC)/$(CU)/omp_loop.c -o $(BUILD_DIR)/$(CU)/omp_loop.o

# Generic spmxv routines
#################################################################################################

$(BUILD_DIR)/spmxv_sequential.o : $(SRC)/$(GENERIC_PATH)/spmxv_sequential.c $(SRC)/$(GENERIC_PATH)/inc/spmxv_sequential.h
	$(CC) $(CFLAGS) $(SRC)/$(GENERIC_PATH)/spmxv_sequential.c -o $(BUILD_DIR)/spmxv_sequential.o

$(BUILD_DIR)/spmxv_omp_loop.o : $(SRC)/$(GENERIC_PATH)/spmxv_omp_loop.c $(SRC)/$(GENERIC_PATH)/inc/spmxv_omp_loop.h
	$(CC) $(CFLAGS) $(SRC)/$(GENERIC_PATH)/spmxv_omp_loop.c -o $(BUILD_DIR)/spmxv_omp_loop.o

$(BUILD_DIR)/spmxv_omp_task.o : $(SRC)/$(GENERIC_PATH)/spmxv_omp_task.c $(SRC)/$(GENERIC_PATH)/inc/spmxv_omp_task.h
	$(CC) $(CFLAGS) $(SRC)/$(GENERIC_PATH)/spmxv_omp_task.c -o $(BUILD_DIR)/spmxv_omp_task.o

$(BUILD_DIR)/spmxv_gws.o : $(SRC)/$(GENERIC_PATH)/spmxv_gws.c $(SRC)/$(GENERIC_PATH)/inc/spmxv_gws.h
	$(CC) $(CFLAGS) $(SRC)/$(GENERIC_PATH)/spmxv_gws.c -o $(BUILD_DIR)/spmxv_gws.o

$(BUILD_DIR)/spmxv_dws.o : $(SRC)/$(GENERIC_PATH)/spmxv_dws.c $(SRC)/$(GENERIC_PATH)/inc/spmxv_dws.h
	$(CC) $(CFLAGS) $(SRC)/$(GENERIC_PATH)/spmxv_dws.c -o $(BUILD_DIR)/spmxv_dws.o

$(BUILD_DIR)/spmxv_static.o : $(SRC)/$(GENERIC_PATH)/spmxv_static.c $(SRC)/$(GENERIC_PATH)/inc/spmxv_static.h
	$(CC) $(CFLAGS) $(SRC)/$(GENERIC_PATH)/spmxv_static.c -o $(BUILD_DIR)/spmxv_static.o

$(BUILD_DIR)/kernel.o : $(SRC)/$(GENERIC_PATH)/kernel.c $(SRC)/$(GENERIC_PATH)/inc/kernel.h
	$(CC) $(CFLAGS) $(SRC)/$(GENERIC_PATH)/kernel.c -o $(BUILD_DIR)/kernel.o
