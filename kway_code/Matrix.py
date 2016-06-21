from subprocess import call
import commands
import os

class Matrix:
        def __init__(self, mmfPath, isSymmetric, isBinary, netType, cachePart):
                self.mmfPath = mmfPath
		self.dir = os.path.dirname(mmfPath)
		self.name = os.path.basename(mmfPath[0:len(mmfPath) - 4])
		self.isSymmetric = isSymmetric
		self.isBinary = isBinary
		self.netType = netType
		self.cachePart = cachePart

	def getVertexOrderPath(self, cacheSize):
		return self.dir + "/" + self.netType + "/" + self.name + "_vertexOrder_" + cacheSize

	def getVertexDimensionPath(self, cacheSize):
		return self.dir + "/" + self.netType + "/" + self.name + "_vertexDimension_" + cacheSize

	def getVertexPartVectorPath(self, cacheSize):
		return self.dir + "/" + self.netType + "/" + self.name + "_vertexPartVector_" + cacheSize

	def getNetOrderPath(self, cacheSize):
		return self.dir + "/" + self.netType + "/" + self.name + "_netOrder_" + cacheSize

	def getNetDimensionPath(self, cacheSize):
		return self.dir + "/" + self.netType + "/" + self.name + "_netDimension_" + cacheSize

	def getNetPartVectorPath(self, cacheSize):
		return self.dir + "/" + self.netType + "/" + self.name + "_netPartVector_" + cacheSize

	def printSubMtxInfo(self):
		print self.name
		for k, v in self.cachePart.iteritems():
			print k, ': ', v, 'sub-matices'

		commandList = self.getKPatohCommandList("kway")
		for i, v in enumerate(commandList):
			print i, ": ", v

	def getKPatohCommandList(self, kWayBinaryPath):
		commandList = []

		for k, v in self.cachePart.iteritems():
			commandList.append("./" + kWayBinaryPath + " " +  self.mmfPath + " COLNET " + str(self.isSymmetric) + " " + str(self.isBinary) + " " + str(v) + " " + self.getVertexOrderPath(k) + " " + self.getVertexDimensionPath(k) + " " + self.getNetOrderPath(k) + " " + self.getNetDimensionPath(k) + " " + self.getVertexPartVectorPath(k) + " " + self.getNetPartVectorPath(k))
		return commandList

	def kWayPartition(self, kWayBinaryPath):
		commandList = self.getKPatohCommandList(kWayBinaryPath)
		for i, v in enumerate(commandList):
			print v
			print commands.getoutput(v)
