from Matrix import Matrix
import subprocess
import sys
import os

# sanity check first
# ------------------------------------------------------------------------------------------------

# check argument count
if(not len(sys.argv) > 3):
    print "Incorrect number of arguments!!"
    print "Usage: " + sys.argv[0] + " matrixDirPath (partitionSize)+"
    sys.exit("Example python " + sys.argv[0] + " /mtx/Freescale1.mtx 1024 512 256 128")

# check matrix directory
rootDir = sys.argv[1]
if(not os.path.isdir(rootDir)):
    sys.exit("Invalid matrix directory, see if the path is correct path!!")

# check mtx_info_reader binary
if(not os.path.isfile("mtx_info_reader")):
    sys.exit("mtx_info_reader binary not found. Put it in the same folder as this script!!")

# check kway (kpatoh) binary
if(not os.path.isfile("kway")):
    sys.exit("kway binary not found. Put it in the same folder as this script!!")

# ------------------------------------------------------------------------------------------------


# get partition sizes from command line parameters
partitionSizes = sys.argv[2:]

# directory that this script resides in
scriptDir = os.path.dirname(os.path.abspath(__file__))

# path of mtx_info_reader binary
mtxInfoReaderPath = scriptDir + "/mtx_info_reader"

#path of kway partition code
kwayPath = scriptDir + "/kway"

# generate matrix list from the output of mtx_info_reader
# ------------------------------------------------------------------------------------------------
matrixList = []
for mtxName in os.listdir(rootDir):

    if(not os.path.isdir(mtxName)):
        continue

    mtxPath = rootDir + mtxName + "/" + mtxName + ".mtx"        

    if(not os.path.isfile(mtxPath)):
        continue

    args = []
    args.append(mtxInfoReaderPath)
    args.append(mtxPath)
    args.extend(partitionSizes)

    mtxInfoReader = subprocess.Popen(args, stdout=subprocess.PIPE)
    status = mtxInfoReader.wait()
    externalPythonScript = mtxInfoReader.stdout.read()
    
    subprocess.call(["mkdir", rootDir + mtxName + "/colnet"])

    if(status == 0):
        exec(externalPythonScript)
    else:
        print "Error: ", mtxName

# call kpatoh to partition each matrix
# ------------------------------------------------------------------------------------------------
for i, matrix in enumerate(matrixList):
    print "============================== " + matrix.name + " =============================="
    matrix.kWayPartition("kway")
    print "------------------------------ " + matrix.name + " ------------------------------"
    print ""
