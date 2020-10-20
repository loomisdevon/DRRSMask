# Dual Rotating Radiation Scattering Mask MCNP Code
# Rotate neutron and gamma flux from 0 to 2PI around mask and return signal received from detector 

# Authors: Ivan Novikov and Devon Loomis



import os
import os.path
import math
import fileinput
import shutil
import subprocess
import io
import string
import numpy as np
import matplotlib.pyplot as plt
import time
import math
from decimal import Decimal
import csv
from tqdm import *

########### GLOBALS #############

#Configuration Name (this is name of input file without .txt)
CONFIG_NAME = 'DDRS3_rand2_absorber1Source'
SOURCE_NAME = 'source6'
tMATRIXCONFIG_NAME = 'DDRS3_rand2_absorber'
tMatrixFilename = tMATRIXCONFIG_NAME + "tMatrix.csv"
#####################
# Source Info
R = 100
Phi = 36
Theta = 45
##############
#################################

###################################### creating all MCNP input files for simulation of the increasing relative distance between source and detector #######################
''' 														Parameters:
														file_name: MCNP input template
														y-pos: y position of center of circle of source path
														r: radius of circle of source path
														init: inital angle away from xy-plane
														final angle away from xy-plane
														angle step size
														limit: longest distance between the two
'''



def smoothing(fluxArray, smoothingParameter):
	smoothingArray = []
	for i in range(len(fluxArray)):
		totSum = 0
		numPnts = 0
		if (i - smoothingParameter < 0):
			for j in range(i,0,-1):
				totSum += fluxArray[j]
				numPnts += 1
			for k in range(i,i+smoothingParameter,1):
				totSum += fluxArray[k]
				numPnts += 1
		elif (i + smoothingParameter > len(fluxArray)):
			for j in range(i,len(fluxArray),1):
				totSum += fluxArray[j]
				numPnts += 1
			for k in range(i,i-smoothingParameter,-1):
				totSum += fluxArray[k]
				numPnts += 1
		else:
			for j in range(i,i+smoothingParameter,1):
				totSum += fluxArray[j]
				numPnts += 1
			for k in range(i,i-smoothingParameter,-1):
				totSum += fluxArray[k]
				numPnts += 1
		average = totSum/numPnts
		smoothingArray.append(average)
	return smoothingArray


def createFiles(file_name, phi, rad, init, final, step_size):
	fileList = []
	marker=0
	rad_phi = math.radians(phi)
	for new_theta in range(init, final, step_size):

		text_search = None
		f =open(file_name)
		for line in f:
			words = line
			sdef = words[0:4]
			if (sdef == "SDEF"):
				text_search = words
				break
		f.close()
		
		rad_theta = math.radians(new_theta)


		x_pos = round(rad * np.cos(rad_theta)*np.sin(rad_phi),3)
		y_pos = round(rad * np.sin(rad_theta)*np.sin(rad_phi),3)
		z_pos = round(rad * np.cos(rad_phi),3)
		r_mag = np.sqrt(x_pos**2+y_pos**2+z_pos**2)
		vecx_pos = round(-x_pos/r_mag,3)
		vecy_pos = round(-y_pos/r_mag,3)
		vecz_pos = round(-z_pos/r_mag,3)
		#theta_rad = np.arctan(z_pos/r)
		#vecz_pos = round(-1 * (theta_rad/(np.pi/2)),5)
		#replacement_text = sdef + " ERG = 1.42 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " VEC= " + str(vecx_pos) + " " + str(vecy_pos) + " " + str(vecz_pos) + " DIR=d1 par=n" + "\n"
		replacement_text = sdef + " ERG = 2.5e-8 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " VEC= " + str(vecx_pos) + " " + str(vecy_pos) + " " + str(vecz_pos) + " DIR=d1 WGT 20 par=n" + "\n"
		#replacement_text = sdef + " ERG = 1.42 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " par=n" + "\n"
		read_name = file_name
		write_name = CONFIG_NAME + SOURCE_NAME + "_" + str(new_theta) + ".txt"

		f1 = open(read_name, 'r')
		f2 = open(write_name, 'w')

		for lines in f1:
			f2.write(lines.replace(text_search, replacement_text))

		f1.close()
		f2.close()
		fileList.append(write_name)

	return (fileList)


################################# delete runtpe files after every set of commands and delete all output files and input files after program run #######################
'''																Parameters
															directory: directory containing all files 
															file: KSEF_2 #####################
															remove_all: test to determine whether to delete all files or only runtpe files
'''

def removeFiles(directory, file1, file2, file3, outfile, initfile, t_file, save_one, remove_all):
	dir_name = directory
	for fname in os.listdir(dir_name):
		if (fname != initfile and fname != t_file):
			if fname.startswith("binRun"):
				os.remove(os.path.join(dir_name, fname))
			if (fname.startswith(file1[:-4]) or fname.startswith(outfile[:-4])) and remove_all:
				if (fname != file1):
					os.remove(os.path.join(dir_name, fname))
			if (fname.startswith(file1[:-4]) or fname.startswith(outfile[:-4])) and save_one: 
				if (fname != file1 and fname != file2 and fname != file3):
					os.remove(os.path.join(dir_name, fname))


####################### read MCNP output file, find and return flux value #########################
#######################_file_: MCNP output file name ##################################
'''
def readFlux(_file_):
	flux_ = 0
	error_ = 0
	with open(_file_, 'r') as outfile:
		for line in outfile:
			if ('+                                   *Gamma flux in detector*' in line):
				lines = [outfile.readline() for i in range(9)] #this reads 9 lines after the fc4 comment
				spectrum = [outfile.readline() for i in range(13)] #this reads 13 lines which contain spectrum
	            #each line has an index [0]-[12]
	    		#print(type(spectrum[1]))
	    		#print(spectrum[1])

	    		#print(float(spectrum[1].split()[1])) #this splits spectrum[i] using spaces
	            #each spectrum[i].split() has three new indeces [0]-[2]
	            #float converts each string to float
	            #Neutron energy is in [0]
	            #Neutron counts are in [1]
	            #Error is in [2]
				#tmp = 0.0
				#print (spectrum)
				for j in range(13):
					flux_ += float(spectrum[j].split()[1])
					error_ += float(spectrum[j].split()[2])

            #Fluxin3[i] = tmp

	return flux_, error_
'''
def readFlux(_file_,energyBin, binWrite):
	flux_Arr = []
	error_Arr = []
	flux_ = 0
	error_ = 0

	with open(_file_, 'r') as outfile:
		for line in outfile:
			if ('+                                   *Neutron Flux In Detector*' in line):

				lines = [outfile.readline() for i in range(9)] #this reads 9 lines after the fc4 comment
				spectrum = [outfile.readline() for i in range(energyBin+1)] #this reads 13 lines which contain spectrum
	            #each line has an index [0]-[12]
	    		#print(type(spectrum[1]))
	    		#print(spectrum[1])

	    		#print(float(spectrum[1].split()[1])) #this splits spectrum[i] using spaces
	            #each spectrum[i].split() has three new indeces [0]-[2]
	            #float converts each string to float
	            #Neutron energy is in [0]
	            #Neutron counts are in [1]
	            #Error is in [2]
				#tmp = 0.0
				for j in range(energyBin+1):
					if (binWrite == j and binWrite != 0):
						flux_ = float(spectrum[j].split()[1])
						error_ = float(spectrum[j].split()[1])
					#flux_Arr.append(float(spectrum[j].split()[1]))
					#error_Arr.append(float(spectrum[j].split()[2]))
					#flux_ += float(spectrum[j].split()[1])
					#error_ += float(spectrum[j].split()[2])
				if (binWrite == 0):
					flux_ = float(spectrum[energyBin].split()[1])
					error_ = float(spectrum[energyBin].split()[2])
            #Fluxin3[i] = tmp

	return flux_, error_

def initialize(_file_):
	global intensity, activity, nps, t
	global radius, init_theta, final_theta, step_theta
	global init_phi, final_phi, step_phi
	global packet


	with open(_file_,"r", newline='') as file:
		file.readline()
		intensity = float(file.readline()[12:])
		activity = float(file.readline()[11:])
		nps = float(file.readline()[6:])
		t = float(file.readline()[4:])
		file.readline()
		radius = float(file.readline()[9:])
		init_theta = int(file.readline()[13:])
		final_theta = int(file.readline()[14:])
		step_theta = int(file.readline()[13:])
		init_phi = float(file.readline()[11:])
		final_phi = float(file.readline()[12:])
		step_phi = float(file.readline()[11:])
		file.readline()
		file.readline()
		packet = int(file.readline()[9:])

#**********************MAIN**************************
#dir_ = 'C:\\Users\\devon\\Documents\\DRRSMask\\Working_Version\\MLEM\\'
dir_ = os.path.dirname(os.path.abspath(CONFIG_NAME+SOURCE_NAME)) + "\\"
file_ = CONFIG_NAME + SOURCE_NAME + '.txt'
outFile_ = CONFIG_NAME + SOURCE_NAME + '_out.txt'
file_name_ = dir_ + file_
outFile_name_ = dir_ + outFile_
keepInFile = CONFIG_NAME + SOURCE_NAME + '_0.txt'
keepOutFile = CONFIG_NAME + SOURCE_NAME + '_out0.txt'
init_file= dir_ + 'init.txt'
t_file = CONFIG_NAME + "tMatrix.csv"



intensity, activity, nps, t = 0,0,0,0
radius,init_theta,final_theta,step_theta = 0,0,0,0
init_phi,final_phi,step_phi = 0,0,0
packet = 0


initialize(init_file)

originalThetaCountsArray = []
originalThetaErrorArray = []



init_theta = Theta


final_theta = init_theta+360
transmissionMatrix = []


start = time.time()
removeFiles(dir_, file_, keepInFile, keepOutFile, outFile_, init_file, t_file, False, True) # purge directory of any existing MCNP files from previous run
#files = createFiles(file_name_, z, radius, init_ang, final_ang, step)
files = createFiles(file_name_, Phi, R, init_theta, final_theta, step_theta) # create all MCNP input files

commands = []
outFileList = []
j = init_theta

#create set of commands for subprocess of all input files
for i in range(int((final_theta - init_theta) / step_theta)):
	binFile = "binRun" + str(j) + ".r"	
	outFile = (CONFIG_NAME + SOURCE_NAME + "_out" + str(j) + ".txt")
	commands.append("mcnp6 i=" + files[i] + " o=" +  outFile + " runtpe=" + binFile)
	outFileList.append(outFile)
	j += step_theta


print("Simulating...")

# give subprocess pak amount of parallel programs to execute until all commands are executed 
for x in tqdm(range(0,int((final_theta - init_theta) / step_theta),(packet))):	
	if (x < (len(commands) - packet)):
		commandsub = commands[x:(x+packet)]
	else:
		commandsub = commands[x:]
	processes = [subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, cwd=dir_) for cmd in commandsub]
	removeFiles(dir_, file_, keepInFile, keepOutFile, outFile_, init_file, t_file, False, False) # remove runtpe files

	for p in processes:
		p.wait()

print ("Checkpoint")
theta = init_theta

fluxList, fluxList1, fluxList2, fluxList3, fluxList4, fluxList5, fluxList6, fluxList7, fluxList8 = [], [], [], [], [], [], [], [], []
errorList, errorList1, errorList2, errorList3, errorList4, errorList5, errorList6, errorList7, errorList8  = [], [], [], [], [], [], [], [], []
sourceThetaList = []

#use for neutrons
#energyBinOfInterest = 13

#use for gammas
energyBinOfInterest = 14

############################read and gather flux values and source distances for each output file and add them to lists###################################
for f in outFileList:
	flux, error = readFlux(f, energyBinOfInterest, 0)
	flux1, error1 = readFlux(f, energyBinOfInterest, 1)
	flux2, error2 = readFlux(f, energyBinOfInterest, 2)
	flux3, error3 = readFlux(f, energyBinOfInterest, 3)
	flux4, error4 = readFlux(f, energyBinOfInterest, 4)
	flux5, error5 = readFlux(f, energyBinOfInterest, 5)
	flux6, error6 = readFlux(f, energyBinOfInterest, 6)
	flux7, error7 = readFlux(f, energyBinOfInterest, 7)
	flux8, error8 = readFlux(f, energyBinOfInterest, 8)
	fluxList.append(flux)
	fluxList1.append(flux1)
	fluxList2.append(flux2)
	fluxList3.append(flux3)
	fluxList4.append(flux4)
	fluxList5.append(flux5)
	fluxList6.append(flux6)
	fluxList7.append(flux7)
	fluxList8.append(flux8)
	errorList.append(error)
	errorList1.append(error1)
	errorList2.append(error2)
	errorList3.append(error3)
	errorList4.append(error4)
	errorList5.append(error5)
	errorList6.append(error6)
	errorList7.append(error7)
	errorList8.append(error8)
	rad_theta = math.radians(theta)
	sourceThetaList.append(rad_theta)
	theta += step_theta

removeFiles(dir_, file_, keepInFile, keepOutFile, outFile_, init_file, t_file, True, False)
end = time.time()
print("Runtime: ", round((end - start)/60, 2), " mins")


#rawFluxArray = np.array(fluxList)
#print (rawFluxArray)
fluxArray = np.array(smoothing(fluxList, 8))
#fluxArray = np.array(fluxList)
errorArray = np.array(errorList)
#print (errorArray)


fluxArray1 = np.array(smoothing(fluxList1, 8))
errorArray1 = np.array(errorList1)

fluxArray2 = np.array(smoothing(fluxList2, 8))
errorArray2 = np.array(errorList2)

fluxArray3 = np.array(smoothing(fluxList3, 8))
errorArray3 = np.array(errorList3)

fluxArray4 = np.array(smoothing(fluxList4, 8))
errorArray4 = np.array(errorList4)

fluxArray5 = np.array(smoothing(fluxList5, 8))
errorArray5 = np.array(errorList5)

fluxArray6 = np.array(smoothing(fluxList6, 8))
errorArray6 = np.array(errorList6)

fluxArray7 = np.array(smoothing(fluxList7, 8))
errorArray7 = np.array(errorList7)

fluxArray8 = np.array(smoothing(fluxList8, 8))
errorArray8 = np.array(errorList8)

angleArray = np.array(sourceThetaList)
countsArray = fluxArray * intensity * t
countsArray1 = fluxArray1 * intensity * t
countsArray2 = fluxArray2 * intensity * t
countsArray3 = fluxArray3 * intensity * t
countsArray4 = fluxArray4 * intensity * t
countsArray5 = fluxArray5 * intensity * t
countsArray6 = fluxArray6 * intensity * t
countsArray7 = fluxArray7 * intensity * t
countsArray8 = fluxArray8 * intensity * t



countsSum = np.sum(countsArray)
countsSum1 = np.sum(countsArray1)
countsSum2 = np.sum(countsArray2)
countsSum3 = np.sum(countsArray3)
countsSum4 = np.sum(countsArray4)
countsSum5 = np.sum(countsArray5)
countsSum6 = np.sum(countsArray6)
countsSum7 = np.sum(countsArray7)
countsSum8 = np.sum(countsArray8)
#normalizedCountsArray = countsArray / countsSum
#normalizedCountsArray = np.copy(countsArray)
normalizedCountsArray = np.copy(countsArray)

normalizedCountsErr = np.sqrt((1/countsArray) + (1/countsSum))
normalizedCountsErrorArray = np.multiply(normalizedCountsArray, normalizedCountsErr)


with open(CONFIG_NAME + SOURCE_NAME + "data.csv","w+", newline='') as file:
	writer=csv.writer(file,delimiter=',')
	for a in normalizedCountsArray:
		writer.writerow([a])

with open(CONFIG_NAME + SOURCE_NAME + "background.csv","w+", newline='') as file:
	
	writer=csv.writer(file,delimiter=',')
	for b in normalizedCountsErrorArray:
		writer.writerow([b])

###########################END MAIN###############################