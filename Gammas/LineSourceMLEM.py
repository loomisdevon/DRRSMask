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
import random as rd
from tqdm import *

########### GLOBALS #############

#Configuration Name (this is name of input file without .txt)
CONFIG_NAME = 'DDRS3_rand2_absorberLineSource'
SOURCE_NAME = 'source3'
tMATRIXCONFIG_NAME = 'DDRS3_rand2_absorber'
tMatrixFilename = tMATRIXCONFIG_NAME + "tMatrix.csv"
############################
# Source Info
z = 0
Theta = 80
zExtension = 30
###############
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


def createFiles(file_name, z_pos, z_ext, r, init, final, step_size):
	fileList = []
	marker=0
	z_ext_pos = z_pos + z_ext
	z_ext_neg = z_pos - z_ext
	for new_ang in range(init, final, step_size):

		text_search = None
		f =open(file_name)
		for line in f:
			words = line
			#sdef = words[0:4]
			startOfWord = words[0:4]
			if (startOfWord == "SDEF"):
				sdef = words
			#if (sdef == "SDEF"):
				#text_search = words
				#break
			elif (startOfWord == "si1 "):
				extent = words
			#elif (startOfWord == "si4 L"):
				#vec2_info = words
			#elif (startOfWord == "si5 L"):
				#vec3_info = words


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
		replacement_sdef = "SDEF ERG 0.661 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " AXS 0 0 1 EXT d1 VEC " + str(vecx_pos) + " " + str(vecy_pos) + " " + str(vecz_pos) + " DIR d2 RAD 0 WGT 20 par=p" + "\n"
		replacement_extent = "si1 " + str(z_ext_neg) + " " + str(z_ext_pos) + "\n"
		#replacement_vec1 = "si3 L " + str(vecx_pos1) + " " + str(vecy_pos1) + " " + str(vecz_pos1)  + "\n"
		#replacement_vec2 = "si4 L " + str(vecx_pos2) + " " + str(vecy_pos2) + " " + str(vecz_pos2)  + "\n"
		#replacement_vec3 = "si5 L " + str(vecx_pos3) + " " + str(vecy_pos3) + " " + str(vecz_pos3)  + "\n"
		#replacement_text = sdef + " ERG = 1.42 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " VEC= " + str(vecx_pos) + " " + str(vecy_pos) + " " + str(vecz_pos) + " DIR=d1 par=n" + "\n"
		#replacement_text = sdef + " ERG = 0.6617 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " VEC= " + str(vecx_pos) + " " + str(vecy_pos) + " " + str(vecz_pos) + " DIR=d1 WGT 20 par=p" + "\n"
		#replacement_text = sdef + " ERG = 1.42 POS " + str(x_pos) + " " + str(y_pos) + " " + str(z_pos) + " par=n" + "\n"
		
		read_name = file_name
		write_name = CONFIG_NAME + SOURCE_NAME + "_" + str(new_ang) + ".txt"

		f1 = open(read_name, 'r')
		f2 = open(write_name, 'w')

		for lines in f1:
			#f2.write(lines.replace(text_search, replacement_text))
			#f2.write(lines.replace(pos_info, replacement_pos).replace(vec1_info, replacement_vec1).replace(vec2_info, replacement_vec2).replace(vec3_info, replacement_vec3))
			f2.write(lines.replace(sdef, replacement_sdef).replace(extent, replacement_extent))

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
			if ('+                                   *Gamma Flux In Detector*' in line):

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


#z = 0
init_ang = Theta
z_extent = zExtension

final_ang = init_ang+360
transmissionMatrix = []


start = time.time()
removeFiles(dir_, file_, keepInFile, keepOutFile, outFile_, init_file, t_file, False, True) # purge directory of any existing MCNP files from previous run
#files = createFiles(file_name_, z, radius, init_ang, final_ang, step)
files = createFiles(file_name_, z, z_extent, rho_, init_ang, final_ang, step) # create all MCNP input files

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




SOURCERUN_NAME = CONFIG_NAME + SOURCE_NAME
dataFilename = SOURCERUN_NAME + "data.csv"
backgroundFilename = SOURCERUN_NAME + "background.csv"

P=49   #Number of Pixels in Z-direction
N=120  # Number of Pixels in Theta-direction
noOfPixels = P*N  #Total Number of Pixels for MLEM of entire surface plane
#noOfPixels = P   #Total Number of Pixels for MLEM of z-position


# Read in all data from csv files saved in same directory as this source code 
transmissionMatrix=list(csv.reader(open(tMatrixFilename)))
data = list(csv.reader(open(dataFilename)))   # signal curve of source location that we are trying to find
background = list(csv.reader(open(backgroundFilename)))   # background (I have made this Error from MCNP?)


################################ This is an attempt at generating MLEM for strictly z-position ###########################
############### This should hopefully work because the generated Transmission Matrix should be handle the "constant theta" ###########


##### Initialize all arrays that will be used (data type is 64 bit floating point)
#mle_n = [rd.random()/100 for i in range(noOfPixels)]
#mle_n = np.array([[rd.random()/100 for x in range(1)] for y in range(N*P)])
mle_n = np.array([[1 for x in range(1)] for y in range(N*P)])
mle_n1 = np.copy(mle_n)
mle_n = np.array(mle_n,dtype=np.float64)
mle_n1 = np.array(mle_n1,dtype=np.float64)
#print (mle_n)
#mle_n = np.array([[1]*noOfPixels], dtype=np.float64).T
#mle_n1 = np.array([[1]*noOfPixels],dtype=np.float64).T
data = np.array(data, dtype=np.float64)
transmissionMatrix = np.array(transmissionMatrix, dtype=np.float64)
background = np.array(background, dtype=np.float64).T
oneDim = np.array([1]*P).T

##### data and background are 2d arrays that are empty except for first index so extract these 1D arrays
#data = data[0]
background = background[0]

##### N dimensional array of ones
identity = np.ones(N).T
identity = np.array([[1 for x in range(1)] for y in range(N)]).T

###### Populate Entire Transmission Matrix ##############
entireTransmissionMatrix = [[0 for x in range(N*P)] for y in range(N)]
for i in range(P):
	signalCurve = transmissionMatrix[:,i]
	for j in range(N):
		colIndex = i*N+j
		shiftedSignalCurve = np.roll(signalCurve,-j)
		for k in range(len(shiftedSignalCurve)):
			entireTransmissionMatrix[k][colIndex] = shiftedSignalCurve[k]
##########################################################

entireTransmissionMatrix = np.array(entireTransmissionMatrix)

onedData = data[:,0]
#print (entireTransmissionMatrix.shape) 
#print (mle_n)
#print (np.dot(entireTransmissionMatrix,mle_n))



n = 0
euclidArray = []
while True:
	condition=False
	n += 1
	
	##### Calculate the mle_(n+1) array from the mle_(n) array
	##### Follow the Presciption of Kowash Pg. 47
	'''
	innerTermSum = np.array([0]*(N*P))
	innerTermSum = np.array([0]*(N*P))
	mleASum = 0
	for a in range(N):
		for b in range(P):
			mleASum += mle_n[:,0][b]*entireTransmissionMatrix[a][b]
		innerTermSumTerm = (np.multiply(entireTransmissionMatrix[a],onedData[a]/mleASum))
		innerTermSum = np.add(innerTermSum,innerTermSumTerm)
	divideTerm = np.array([0]*(N*P))


	for c in range(N):
		divideTerm = np.add(entireTransmissionMatrix[c],divideTerm)
	factor = np.divide(innerTermSum,divideTerm)
	'''
	
	innerTerm = np.divide(data,(np.dot(entireTransmissionMatrix,mle_n)))

	innerWithT = np.dot(entireTransmissionMatrix.T, innerTerm)
	#print (innerWithT)
	#print (np.dot(identity,entireTransmissionMatrix))
	factor = np.array(np.divide(innerWithT,np.dot(identity,entireTransmissionMatrix).T))
	#print (factor.shape)
	
	#print (factor)
	#temp = np.multiply(mle_n[:,0], factor)
	#print (mle_n[:,0])
	temp = np.multiply(mle_n, factor)
	#print (temp.shape)
	#print (temp)
	for t in range(len(mle_n1)):
		value = temp[t]
		mle_n1[t][0]= value
		#mle_n1[t]= value
	
	#print (mle_n)
	#print (mle_n1)
	###### Check the convergence of mle_(n+1)
	###### If mle_(n+1) has not converged, then assign mle_(n+1) to mle_n and loop again
	#print (mle_n)
	#print (mle_n1)
	#print (np.subtract(mle_n[:,0],mle_n1[:,0]))
	epsilon = 0.1
	euclidTotal = 0
	for i in range(len(mle_n1)):
		#print (abs(mle_n1[i][0] - mle_n[i][0]))
		#print (np.average(mle_n1[:,0]))
		# Euclidean Distance
		#if (n == 100):
			#print ((mle_n1[i][0]-mle_n[i][0])**2)
			#if (((mle_n1[i][0]-mle_n[i][0])**2)>10000):
				#print (mle_n1[i][0])
				#print (mle_n[i][0])
				#print (i)
			#print (euclidTotal)
		euclidTotal += (mle_n1[i][0]-mle_n[i][0])**2


		#if (abs(mle_n1[i][0] - mle_n[i][0]) <= np.average(mle_n1[:,0])*0.00005):
			#condition=True
		#else:
			#condition=False
			#continue
	
	#if (condition==True):
	#print (np.sqrt(euclidTotal))
	euclidArray.append((np.sqrt(euclidTotal)))
	if (n%100==0):
		print (n, ': ',np.sqrt(euclidTotal))
	if (np.sqrt(euclidTotal) < 1e-6):
		print ("MLE has converged!\nIterations: ", n)
		break
	elif (n > 5000):
		print ("MLE has not converged in 1000 iterations. Exiting.")
		break
	else:
		mle_n = np.copy(mle_n1)


'''
with open(SOURCERUN_NAME + "Convergence.csv","w+", newline='') as file:
	
	writer=csv.writer(file,delimiter=',')
	for a in euclidArray:
		writer.writerow([a])
'''		
##### Plot MLE along z-axis
##### Should see a maximum at the z that corresponds to the source position 
#z_val = [z for z in range(-100,101,5)]
#plt.plot(z_val,mle_n1[:,0])
#plt.show()

zMax = 120
zStep = 5
noOfZSteps = P
noOfAngleSteps = 120
num = 0
angleArr = []
totalArr = []

#print (mle_n1)

for x in range(noOfZSteps):
	angleArr = []
	for y in range(noOfAngleSteps):
		#angleArr.append(mle_n1[num][0])
		angleArr.append(mle_n1[num][0])
		num += 1
	totalArr.append(angleArr)
for z in range(noOfZSteps):
	print ("Z: ", (5*z-120), " Intensity: ",np.max(totalArr[z]))

#print (np.array(totalArr[][])

fig1, (ax1,ax2) = plt.subplots(1,2,constrained_layout=True,sharey=True)

zContour = np.arange(-zMax,zMax+1,5)
thetaContour = np.arange(0,360,3)

thetaCONT, zCONT = np.meshgrid(thetaContour,zContour)

levels = [1e-3,5e-3,1e-2,2e-2,3e-2]
cmap = plt.cm.get_cmap("hot")
cmap.set_under("magenta")
cmap.set_over("yellow")
#CS = ax2.contour(thetaCONT,zCONT,totalArr,levels,cmap=cmap)
CS = ax2.contour(thetaCONT,zCONT,totalArr,cmap=cmap)
CS.cmap.set_under("gray")
CS.cmap.set_over("yellow")
#plt.clabel(CS, fmt='%1.2e', colors='black', fontsize=4)
fig1.colorbar(CS)
#plt.contour(totalArr)

#plt.imshow(totalArr, cmap=plt.cm.jet, origin='lower', extent=[0,6.28,-zMax,zMax], aspect='auto')
ax1.imshow(totalArr, cmap=plt.cm.jet, origin='lower', extent=[0,360,-zMax,zMax], aspect='auto')
#plt.imshow(totalArr, cmap=plt.cm.jet, origin='lower')
ax1.set_xlabel("Theta (degrees)")
ax2.set_xlabel("Theta (degrees)")
ax1.set_ylabel("Z (cm)")
plt.show()
###########################END MAIN###############################