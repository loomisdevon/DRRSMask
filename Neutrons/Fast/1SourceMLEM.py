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
SOURCE_NAME = 'source3'
tMATRIXCONFIG_NAME = 'DDRS3_rand2_absorber1Source'
tMatrixFilename = tMATRIXCONFIG_NAME + "tMatrix.csv"
#####################
# Source Info
z = 60
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


def createFiles(file_name, z_pos, r, init, final, step_size):
	fileList = []
	marker=0
	for new_ang in range(init, final, step_size):

		text_search = None
		f =open(file_name)
		for line in f:
			words = line
			sdef = words[0:4]
			if (sdef == "SDEF"):
				text_search = words
				break
		f.close()
		
		rad_ang = math.radians(new_ang)

		x_pos = round(r * np.cos(rad_ang),1)
		y_pos = round(r * np.sin(rad_ang),1)
		r_mag = np.sqrt(x_pos**2+y_pos**2+z_pos**2)
		vecx_pos = round(-x_pos/r_mag,1)
		vecy_pos = round(-y_pos/r_mag,1)
		vecz_pos = round(-z_pos/r_mag,1)
		theta_rad = np.arctan(z_pos/r)
		#vecx_pos = round(-1 * np.cos(rad_ang),5)
		#vecy_pos = round(-1 * np.sin(rad_ang),5)
		#theta_rad = np.arctan(z_pos/r)
		#vecz_pos = round(-1 * (theta_rad/(np.pi/2)),5)
		#replacement_text = sdef + " ERG = 1.42 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " VEC= " + str(vecx_pos) + " " + str(vecy_pos) + " " + str(vecz_pos) + " DIR=d1 par=n" + "\n"
		replacement_text = sdef + " ERG = 2.00 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " VEC= " + str(vecx_pos) + " " + str(vecy_pos) + " " + str(vecz_pos) + " DIR=d1 WGT 20 par=n" + "\n"
		#replacement_text = sdef + " ERG = 1.42 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " par=n" + "\n"
		read_name = file_name
		write_name = CONFIG_NAME + SOURCE_NAME + "_" + str(new_ang) + ".txt"

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
def readFlux(_file_,energyBin):
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
					flux_Arr.append(float(spectrum[j].split()[1]))
					error_Arr.append(float(spectrum[j].split()[2]))
					#flux_ += float(spectrum[j].split()[1])
					#error_ += float(spectrum[j].split()[2])
				flux_ = float(spectrum[energyBin].split()[1])
				error_ = float(spectrum[energyBin].split()[2])
            #Fluxin3[i] = tmp

	return flux_, error_

def initialize(_file_):
	global intensity, activity, nps, t
	global rho_,init_ang,final_ang,step
	global zMax,zMin,zStep,initialZFlag,initialZ
	global packet


	with open(_file_,"r", newline='') as file:
		file.readline()
		intensity = float(file.readline()[12:])
		activity = float(file.readline()[11:])
		nps = float(file.readline()[6:])
		t = float(file.readline()[4:])
		file.readline()
		rho_ = float(file.readline()[7:])
		init_ang = int(file.readline()[11:])
		final_ang = int(file.readline()[12:])
		step = int(file.readline()[7:])
		zMax = float(file.readline()[7:])
		zMin = float(file.readline()[7:])
		zStep = float(file.readline()[8:])
		file.readline()
		initialZFlag = bool(file.readline()[15:])
		initialZ = float(file.readline()[11:])
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
rho_,init_ang,final_ang,step = 0,0,0,0
zMax,zMin,zStep,initialZFlag,initialZ = 0,0,0,0,0
packet = 0


initialize(init_file)

originalZCountsArray = []
originalZErrorArray = []



init_ang = Theta


final_ang = init_ang+360
transmissionMatrix = []


start = time.time()
removeFiles(dir_, file_, keepInFile, keepOutFile, outFile_, init_file, t_file, False, True) # purge directory of any existing MCNP files from previous run
#files = createFiles(file_name_, z, radius, init_ang, final_ang, step)
files = createFiles(file_name_, z, rho_, init_ang, final_ang, step) # create all MCNP input files

commands = []
outFileList = []
j = init_ang

#create set of commands for subprocess of all input files
for i in range(int((final_ang - init_ang) / step)):
	binFile = "binRun" + str(j) + ".r"	
	outFile = (CONFIG_NAME + SOURCE_NAME + "_out" + str(j) + ".txt")
	commands.append("mcnp6 i=" + files[i] + " o=" +  outFile + " runtpe=" + binFile)
	outFileList.append(outFile)
	j += step


print("Simulating...")

# give subprocess pak amount of parallel programs to execute until all commands are executed 
for x in tqdm(range(0,int((final_ang - init_ang) / step),(packet))):	
	if (x < (len(commands) - packet)):
		commandsub = commands[x:(x+packet)]
	else:
		commandsub = commands[x:]
	processes = [subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, cwd=dir_) for cmd in commandsub]
	removeFiles(dir_, file_, keepInFile, keepOutFile, outFile_, init_file, t_file, False, False) # remove runtpe files

	for p in processes:
		p.wait()

print ("Checkpoint")
ang = init_ang

fluxList = []
errorList = []
sourceAngleList = []

#use for neutrons
#energyBinOfInterest = 13

#use for gammas
energyBinOfInterest = 14

############################read and gather flux values and source distances for each output file and add them to lists###################################
for f in outFileList:
	flux, error = readFlux(f, energyBinOfInterest)
	fluxList.append(flux)
	errorList.append(error)
	rad_ang = math.radians(ang)
	sourceAngleList.append(rad_ang)
	ang += step

removeFiles(dir_, file_, keepInFile, keepOutFile, outFile_, init_file, t_file, True, False)
end = time.time()
print("Runtime: ", round((end - start)/60, 2), " mins")



rawFluxArray = np.array(fluxList)
#print (rawFluxArray)
fluxArray = np.array(smoothing(fluxList, 8))
#fluxArray = np.array(fluxList)
errorArray = np.array(errorList)
#print (errorArray)

angleArray = np.array(sourceAngleList)
countsArray = fluxArray * intensity * t

countsSum = np.sum(countsArray)
#normalizedCountsArray = countsArray / countsSum
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