import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['animation.convert_path'] = 'C:\\Program Files\\ImageMagick-7.0.10-Q16-HDRI\\magick.exe'
import csv
import random as rd
import os
from tqdm import *
from matplotlib.animation import FuncAnimation


def smoothing(fluxArray, smoothingParameter, colVector):
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

	smoothingArrayCol = [[0 for x in range(1)] for y in range(len(smoothingArray))]
	for j in range(len(smoothingArray)):
		smoothingArrayCol[j][0] = smoothingArray[j]
	
	if colVector:
		return smoothingArrayCol
	else:
		return smoothingArray




fig = plt.figure()
ax = plt.axes(xlim=(0, 4), ylim=(-2, 2))
line, = ax.plot([], [], lw=3)




CONFIG_NAME = 'DDRS3_rand2_absorber1Source'
tMatrixFilename = CONFIG_NAME + "tMatrix.csv"

SOURCERUN_NAME = ".//Source10/" + CONFIG_NAME + "source10"
dataFilename = SOURCERUN_NAME + "data.csv"
backgroundFilename = SOURCERUN_NAME + "background.csv"

actualPhi = 90
actualTheta = 180
time = 30.
normalizedTime = time/5.

#P=60   #Number of Pixels in Phi-direction
P=51 
N=360  # Number of Pixels in Theta-direction
noOfPixels = P*N  #Total Number of Pixels for MLEM of entire surface plane
#noOfPixels = P   #Total Number of Pixels for MLEM of z-position

# Read in all data from csv files saved in same directory as this source code 
transmissionMatrix=list(csv.reader(open(tMatrixFilename)))
data = list(csv.reader(open(dataFilename)))   # signal curve of source location that we are trying to find
background = list(csv.reader(open(backgroundFilename)))   # background (I have made this Error from MCNP?)

#print (data)
#data = smoothing(data.T,7).T
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

#CLYC gamma detection efficiency
data = np.multiply(data, 0.1)

#Multiply by irradiation time
data = np.multiply(data, time)

#Multiply by active area of detector
data = np.multiply(data, 20)

#Introduce Poisson Noise from Counts (This scales inversely with irradiation time)
for inoise in range(len(data[:,0])):
	data[inoise][0] = np.random.normal(data[inoise][0],np.sqrt(data[inoise][0]))

#Round because detector counts is discrete
data = np.round(data)

#Normalize Signal Curve
data /= sum(data)[0]



#print(data[:,0])




#print (data)
#print (data/sum(data)[0])
#print (sum(data))
#data = np.array(smoothing(data.T[0],7,1))

#print (data[:,0])
#print (np.array(smoothing(data.T[0],7,1)).T)
transmissionMatrix = np.array(transmissionMatrix, dtype=np.float64)

#CLYC gamma detection efficiency
transmissionMatrix = np.multiply(transmissionMatrix, 0.15)

#Multiply by irradiation time
transmissionMatrix = np.multiply(transmissionMatrix, time)

#Discrete Detector Counts
transmissionMatrix = np.round(transmissionMatrix)


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
for i in tqdm(range(P)):
	signalCurve = transmissionMatrix[:,i]
	for j in range(N):
		colIndex = i*N+j
		#hiftedSignalCurve = smoothing(np.roll(signalCurve,-j),7,0)
		shiftedSignalCurve = np.roll(signalCurve,-j)
		for k in range(len(shiftedSignalCurve)):
			entireTransmissionMatrix[k][colIndex] = shiftedSignalCurve[k]/(sum(shiftedSignalCurve))
			#entireTransmissionMatrix[k][colIndex] = shiftedSignalCurve[k]
##########################################################
#print ("Here")
entireTransmissionMatrix = np.array(entireTransmissionMatrix)

existingTransmissionMatrix = False
'''
for i in range(N):
	with open("smoothedTMatrix.csv","r") as file_input:
			reader=csv.reader(file_input)
			for row in reader:
				if (not row[0]):
					break
				else:
					existingTransmissionMatrix = True
					break
			file_input.seek(0)
			with open("tempSmoothedTMatrix.csv","w") as file_output:
				writer=csv.writer(file_output,lineterminator='\n')
				if (existingTransmissionMatrix):
					appendIndex = 0
					for row in reader:
						row.append(entireTransmissionMatrix[i][appendIndex])
						writer.writerow(row)
						appendIndex += 1
				elif (existingTransmissionMatrix == False):
					for index in entireTransmissionMatrix[:,0]:
						writer.writerow([index])
	os.remove("smoothedTMatrix.csv")
	os.rename("tempSmoothedTMatrix.csv", "smoothedTMatrix.csv")
'''
'''
with open("smoothedTMatrix.csv","w+", newline='') as file:
	writer=csv.writer(file,delimiter=',')
	for a in data[:,0]:
		writer.writerow([a])
'''

onedData = data[:,0]
#print (entireTransmissionMatrix.shape) 
#print (mle_n)
#print (np.dot(entireTransmissionMatrix,mle_n))



n = 0
euclidArray = []



def init():
    line.set_data([], [])
    return line,



def animate(i):
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
	phi = np.linspace(0, 153, 51)
	mle = []
	for phi_rms in range(P):
		for theta_rms in range(N):
			if (phi_rms*3 == rms_phi):
				mle.append(mle_n[i][0])

	line.set_data(phi, mle)
	return line,

	'''	
	if (n%1000==0):
		#for ()
		#plt.plot(mle_n)
		#plt.show()
		print (n, ': ',np.sqrt(euclidTotal))
	if (np.sqrt(euclidTotal) < 1e-6):
		print ("MLE has converged!\nIterations: ", n)
		break
	elif (n > 5000):
		print ("MLE has not converged in 3000 iterations. Exiting.")
		break
	else:
		'''
	mle_n = np.copy(mle_n1)


'''
with open(CONFIG_NAME + SOURCE_NAME + "Convergence.csv","w+", newline='') as file:
	
	writer=csv.writer(file,delimiter=',')
	for a in euclidArray:
		writer.writerow([a])
'''		
##### Plot MLE along z-axis
##### Should see a maximum at the z that corresponds to the source position 
#z_val = [z for z in range(-100,101,5)]
#plt.plot(z_val,mle_n1[:,0])
#plt.show()

#phiMax = 180
phiMax = 153
phiStep = 3
noOfPhiSteps = P
noOfThetaSteps = 360
num = 0
thetaArr = []
totalArr = []
C = 0
#print (mle_n1)
max_phi_ind = 0
max_theta_ind = 0

for x in range(noOfPhiSteps):
	thetaArr = []
	for y in range(noOfThetaSteps):
		#angleArr.append(mle_n1[num][0])
		thetaArr.append(mle_n1[num][0])
		C += mle_n1[num][0]
		if (mle_n1[num][0] > C):
			#C = mle_n1[num][0]
			max_phi_ind = x
			max_theta_ind = y
		num += 1
		
	totalArr.append(thetaArr)
for phi in range(noOfPhiSteps):
	print ("Phi: ", (3*phi), " Intensity: ", np.max(totalArr[phi]))



print ("Max Phi:", max_phi_ind*3, " Max Theta: ", max_theta_ind)
print ("C: ", C)
#print (np.array(totalArr[][])
num_rms = 0
rms_phi = actualPhi - actualPhi % 3
rms_sum = 0
snr_num_sum = 0
snr_den_sum = 0
for phi_rms in range(P):
	for theta_rms in range(N):
		snr_num_sum += mle_n1[num_rms][0]*mle_n1[num_rms][0]
		if (phi_rms*3 == rms_phi and theta_rms == actualTheta):
			rms_sum += (mle_n1[num_rms][0] - C)*(mle_n1[num_rms][0] - C)
			snr_den_sum += (mle_n1[num_rms][0] - C)*(mle_n1[num_rms][0] - C)
		else:
			rms_sum += mle_n1[num_rms][0]*mle_n1[num_rms][0]
			snr_den_sum += mle_n1[num_rms][0]*mle_n1[num_rms][0]
		num_rms += 1

anim = FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)


writer = ImageMagickFileWriter()
anim.save('mle_anim.gif', writer=writer)
'''
#snr = 10*np.log10(snr_num_sum/(snr_den_sum*N*P))
#snr = 10*np.log10(snr_num_sum/(snr_den_sum/(N*P)))
snr = 10*np.log10(snr_num_sum/(snr_den_sum/(N*P)))
print ("SNR: ", snr)
#fig1, (ax1,ax2) = plt.subplots(1,2,constrained_layout=True,sharey=True)

#zContour = np.arange(0,180,3)
#zContour = np.arange(24,138,3)
zContour = np.arange(0,153,3)
thetaContour = np.arange(0,360,1)

thetaCONT, zCONT = np.meshgrid(thetaContour,zContour)


sigma = np.sqrt(rms_sum/(N*P))
print (sigma)
OneSigma = np.max(totalArr) - 1*sigma
TwoSigma = np.max(totalArr) - 2*sigma
ThreeSigma = np.max(totalArr) - 3*sigma
print (OneSigma, " ", TwoSigma, " ", ThreeSigma)

#levels = [1e-3,5e-3,1e-2,2e-2,3e-2]
#levels = [1e-3,1.125e-3,1.25e-3]
levels = [ThreeSigma,TwoSigma,OneSigma,C]
cmap = plt.cm.get_cmap("hot")
cmap.set_under("magenta")
cmap.set_over("yellow")
'''
'''
CS = ax2.contour(thetaCONT,zCONT,totalArr,levels,cmap=cmap)
#CS = ax2.contour(thetaCONT,zCONT,totalArr,cmap=cmap)
CS.cmap.set_under("gray")
CS.cmap.set_over("yellow")
#plt.clabel(CS, fmt='%1.2e', colors='black', fontsize=4)
fig1.colorbar(CS)
#plt.contour(totalArr)
'''

#plt.imshow(totalArr, cmap=plt.cm.jet, origin='lower', extent=[0,6.28,-zMax,zMax], aspect='auto')
'''
ax1.imshow(totalArr, cmap=plt.cm.jet, origin='lower', extent=[0,360,0,153], aspect='auto')
ax1.plot(actualTheta,actualPhi,'kx',markersize=8)
ax1.plot(actualTheta,actualPhi,'ko',markersize=8,markerfacecolor='none')
'''
'''
plt.imshow(totalArr, cmap=plt.cm.jet, origin='lower', extent=[0,360,0,153], aspect='auto')
plt.plot(actualTheta,actualPhi,'kx',markersize=8)
plt.plot(actualTheta,actualPhi,'ko',markersize=8,markerfacecolor='none')
#ax2.plot(actualTheta,actualPhi,'bx',markersize=2)
#plt.imshow(totalArr, cmap=plt.cm.jet, origin='lower')
#ax1.set_xlabel("Theta (degrees)")
plt.xlabel("Theta (degrees)")
#ax2.set_xlabel("Theta (degrees)")
#ax1.set_ylabel("Phi (degrees)")
plt.ylabel("Theta (degrees)")
plt.show()
'''