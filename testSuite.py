import multivector
import matplotlib
import time
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d import Axes3D
multi = multivector.multi

class varyingBFieldTestSession:
    
    initialRotor = multi(bincodes = [0], quickScales = [1])
    initialSpinByRotor = multi(bincodes = [4], quickScales = [1])
    initialSpinByCross = [0, 0, 1]
    I = multi(bincodes = [7], quickScales = [1])
    e1 = multi(bincodes = [1], quickScales = [1])
    e2 = multi(bincodes = [2], quickScales = [1])
    e3 = multi(bincodes = [4], quickScales = [1])
    omega = 1
    B0 = 1
    B1 = 1
    
    
    def __init__(self, timestepList, iterations):
        self.iterations = iterations
        self.timestepList = timestepList

    def runAllTestsAndPlotAllGraphs(self, plt):
        self.runAllErrorTests()
        self.plotAllSpins(plt)
        self.plotErrorGraphs(plt)
        plt.show()
        
    def runAllErrorTests(self):
        self.rotorErrors = []
        self.crossProductErrors = []
        self.exactRotorSpins = []
        self.approximateRotorSpins = []
        self.approximateCrossProductSpins = []
        
        for timestep in self.timestepList:
            testResults = self.runRotorAndCrossProductTest(timestep, self.iterations)
            self.rotorErrors.append(testResults[0])
            self.crossProductErrors.append(testResults[1])
            self.exactRotorSpins.append(testResults[2])
            self.approximateRotorSpins.append(testResults[3])
            self.approximateCrossProductSpins.append(testResults[4])
        
    def runRotorAndCrossProductTest(self, timestep, iterations):
        exactRotorResults = self.runExactRotorTest(timestep, iterations)
        
        approximateRotorResults = self.runApproximateRotorTest(timestep, iterations)
        rotorMaximumErrors = self.getMaximumErrors(exactRotorResults, approximateRotorResults)
        
        approximateCrossProductResults = self.runApproximateCrossProductTest(timestep, iterations)
        crossProductMaximumErrors = self.getMaximumErrors(exactRotorResults, approximateCrossProductResults)

        return (rotorMaximumErrors, crossProductMaximumErrors, exactRotorResults, approximateRotorResults, approximateCrossProductResults)


    def getMaximumErrors(self, exactArray, approximateArray):
        length = len(exactArray)
        squareDifferences = (exactArray - approximateArray)**2
        maxErrors = [0]
        
        for i in range(len(squareDifferences)):
            currentError = (squareDifferences[i,:].sum())**0.5
            
            if currentError > max(maxErrors):
                maxErrors.append(currentError)
            else:
                maxErrors.append(max(maxErrors))
                
        return maxErrors[2:]
        

    def runApproximateRotorTest(self, timestep, iterations):
        rotorsList, spinsList = [self.initialRotor],[self.initialSpinByRotor]
        spinResultVectors = np.array([[0,0,1]])
        for i in range(iterations):
            currentRotor = rotorsList[-1]
            currentRotor = currentRotor*(1/currentRotor.getRotorNorm())
            rotorsList.append(currentRotor + 0.5*(self.I)*(self.getBfield((i+1)*timestep,returnMulti = True))*currentRotor*timestep)
            spinsList.append(spinsList[0].applyRotorPair(currentRotor))
            spinResultVectors = np.append(spinResultVectors, [[spinsList[-1].getBladeScale(1), spinsList[-1].getBladeScale(2), spinsList[-1].getBladeScale(4)]], axis = 0)
        
        return spinResultVectors

    def runApproximateCrossProductTest(self, timestep, iterations):
        spinResultVectors = np.array([self.initialSpinByCross])
        
        for i in range(iterations):
            spinDot = self.cross(spinResultVectors[i, :], self.getBfield((i+1)*timestep))
            spinDot = np.array(spinDot)
            newSpin = spinDot*timestep + spinResultVectors[-1]
            spinResultVectors = np.append(spinResultVectors, [newSpin], axis = 0)


        return spinResultVectors
            

    def cross(self, vector1, vector2):
        xComponent = vector1[1]*vector2[2]-vector1[2]*vector2[1]
        yComponent = -vector1[0]*vector2[2]+vector1[2]*vector2[0]
        zComponent = vector1[0]*vector2[1]-vector1[1]*vector2[0]
        return [xComponent, yComponent, zComponent]

    def runExactRotorTest(self, timestep, iterations):
        rotorsList, spinsList = [self.initialRotor],[self.initialSpinByRotor]
        spinResultVectors = np.array([[0,0,1]])
        for i in range(iterations):
            currentTime = ((i+1)*timestep)
            
            bivector1 = -0.5*self.omega*self.I*self.e3
            bivector1 = bivector1.getMagnitudeOfOrder(2), bivector1*(1/bivector1.getMagnitudeOfOrder(2))
            currentRotor = math.cos(bivector1[0]*currentTime) + bivector1[1]*math.sin(bivector1[0]*currentTime)
            
            bivector2 = (0.5*(self.omega+self.B0)*self.I*self.e3) + (0.5*self.B1*self.I*self.e1)
            bivector2 = bivector2.getMagnitudeOfOrder(2), bivector2*(1/bivector2.getMagnitudeOfOrder(2))
            currentRotor = currentRotor*( math.cos(bivector2[0]*currentTime) + bivector2[1]*math.sin(bivector2[0]*currentTime) )
                        
            spinsList.append(spinsList[0].applyRotorPair(currentRotor))
            spinResultVectors = np.append(spinResultVectors, [[spinsList[-1].getBladeScale(1), spinsList[-1].getBladeScale(2), spinsList[-1].getBladeScale(4)]], axis = 0)

        return spinResultVectors

    def getBfield(self, currentTime, returnMulti = False):
        currentBField = [self.B1*math.cos(self.omega*currentTime), self.B1*math.sin(self.omega*currentTime), self.B0]
        if returnMulti:
            return multi(bincodes = [1,2,4], quickScales = currentBField)
        return currentBField

    
    def plotAllSpins(self, plt):
        fig = plt.figure(1)
        ax = fig.add_subplot(131, projection='3d', aspect = 'equal')
        ax.plot(self.exactRotorSpins[-1][:, 0], self.exactRotorSpins[-1][:, 1], self.exactRotorSpins[-1][:, 2], 'r')
        ax.set_title("Exact solution")
        ax.set_xlabel('x-component')
        ax.set_ylabel('y-component')
        ax.set_zlabel('z-component')
        ax = fig.add_subplot(132, projection='3d')
        ax.plot(self.approximateRotorSpins[-1][:, 0], self.approximateRotorSpins[-1][:, 1], self.approximateRotorSpins[-1][:, 2], 'r')
        ax.set_title("Rotor approximation")
        ax.set_xlabel('x-component')
        ax.set_ylabel('y-component')
        ax.set_zlabel('z-component')
        ax = fig.add_subplot(133, projection='3d')
        ax.plot(self.approximateCrossProductSpins[-1][:, 0], self.approximateCrossProductSpins[-1][:, 1], self.approximateCrossProductSpins[-1][:, 2], 'r')
        ax.set_title("CP approximation")
        ax.set_xlabel('x-component')
        ax.set_ylabel('y-component')
        ax.set_zlabel('z-component')
        
    def plotErrorGraphs(self, plt):
        plt.figure(2)
        self.plotErrorVsIterations(plt)
        self.plotErrorVsTimestepSize(plt)
        self.plotRotorErrorContour(plt)
        self.plotCrossProductErrorContour(plt)

    def plotErrorVsIterations(self, plt):
        plt.subplot(221)
        plt.plot(np.arange(0, self.iterations, 1), self.rotorErrors[0], "g", label = "Rotor")
        plt.plot(np.arange(0, self.iterations, 1), self.crossProductErrors[0], "b", label = "Cross-product")
        plt.title("Error vs Timesteps elapsed")
        plt.xlabel("Timesteps", fontsize=12)
        plt.ylabel("Error", fontsize=12)
        rotorLegendHandle, crossProductLegendHandle = mlines.Line2D([], [], color='green', label='Rotor'), mlines.Line2D([], [], color='blue', label='Cross-Product')
        plt.legend(handles = [rotorLegendHandle, crossProductLegendHandle], loc = 2, fontsize = 9)
        
    def plotErrorVsTimestepSize(self, plt):
        allErrorsTuple = self.extractErrors()
        plt.subplot(222)
        plt.plot(self.timestepList, allErrorsTuple[0], "g", label = "Rotor")
        plt.plot(self.timestepList, allErrorsTuple[1], "b", label = "Cross-product")
        plt.title("Error vs Timestep size")
        plt.xlabel("Timestep size (s)", fontsize=12)
        plt.ylabel("Error", fontsize=12)
        rotorLegendHandle, crossProductLegendHandle = mlines.Line2D([], [], color='green', label='Rotor'), mlines.Line2D([], [], color='blue', label='Cross-Product')
        plt.xticks(self.timestepList)
        plt.legend(handles = [rotorLegendHandle], loc = 2, fontsize = 9)

    def extractErrors(self):
        rotorErrors = np.array(self.rotorErrors)[:, -1]
            
        crossProductErrors = []
        for crossProductErrorList in self.crossProductErrors:
            crossProductErrors.append(crossProductErrorList[-1])

        return (rotorErrors, crossProductErrors)

    def plotRotorErrorContour(self, plt):
        rotorErrorArray = np.array(self.rotorErrors).T
        plt.subplot(223)
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        handle = plt.contour(self.timestepList, np.arange(0, self.iterations, 1), rotorErrorArray, colors = 'g')
        manual_locations = [(0.00119, 2000), (0.0015, 2000), (0.0018, 2000), (0.0021, 2000), (0.0024, 2000), (0.0027, 2000), (0.003, 2000)]
        plt.clabel(handle, fmt = '%.2E', inline=1, fontsize=10, manual = manual_locations)
        self.plotReciprocalContour(plt)
        plt.title("Rotor error contour of timestep size and iterations")
        plt.xlabel("Size of timestep (s)", fontsize=12)
        plt.ylabel("Iterations", fontsize=12)
        plt.xticks(self.timestepList)

    def plotCrossProductErrorContour(self, plt):
        crossProductErrorArray = np.array(self.crossProductErrors).T
        plt.subplot(224)
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        handle = plt.contour(self.timestepList, np.arange(0, self.iterations, 1), crossProductErrorArray, colors = 'b')
        manual_locations = [(0.002, 3000), (0.0023, 4226), (0.0025, 5378), (0.00265, 6872), (0.00272, 8135), (0.0027887, 9380)]
        plt.clabel(handle, fmt = '%.0E', inline=1, fontsize=10, manual = manual_locations)
        self.plotReciprocalContour(plt)
        plt.title("CP error countour of timestep size and iterations")
        plt.xlabel("Size of timestep (s)", fontsize=12)
        plt.ylabel("Iterations", fontsize=12)
        plt.xticks(self.timestepList)
        
    def plotReciprocalContour(self, plt):
        fullGridPoints = []
        xGridSpacing = np.arange(min(self.timestepList), max(self.timestepList), (max(self.timestepList)-min(self.timestepList))/20 )
        yGridSpacing = np.arange(0, self.iterations, self.iterations/20 )
        for n in yGridSpacing:
            fullGridPoints.append([])
            for i in xGridSpacing:
                fullGridPoints[-1].append(n*i)
        plt.contour(xGridSpacing, yGridSpacing, fullGridPoints, colors = 'k')

        
    def squareWave(self, currentTime, period):
        return (1+(-1)**(int(currentTime/(period/2))))/2


    def createImages(self, plt, filePath):
        fig = plt.figure(1).clf()
        os.chdir(filePath)
        for i in range(len(self.exactRotorSpins[1])):
            fig = plt.figure(1)
            ax = fig.add_subplot(111, projection='3d', aspect='equal')
            ax.plot(self.exactRotorSpins[-1][:i+1,0], self.exactRotorSpins[-1][:i+1,1], self.exactRotorSpins[-1][:i+1,2], 'r')
            ax.set_xlim3d(-1, 1)
            ax.set_ylim3d(-1, 1)
            ax.set_zlim3d(0, 1)
            plt.savefig(str(i) + '.png')
            fig = plt.figure(1).clf()
            print(str(i) + " done.")

    def outputSpin(self, filename, path):
        os.chdir(path)
        fileHandle = open(filename, 'w')
        for n in self.exactRotorSpins[-1]:
            fileHandle.write(str(n)+"\n")
            print(n)
        fileHandle.close()
        
session = varyingBFieldTestSession([0.001, 0.002, 0.003, 0.004], 1000)
session.runAllTestsAndPlotAllGraphs(plt)
