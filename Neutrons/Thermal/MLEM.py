import numpy as np
import matplotlib.pyplot as plt
import csv
import random as rd


CONFIG_NAME = 'DDRS3_rand2_absorber1Source'
tMatrixFilename = CONFIG_NAME + "tMatrix.csv"

SOURCERUN_NAME = ".//Source2//" + CONFIG_NAME + "source2"
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
	if (n%1000==0):
		print (n, ': ',np.sqrt(euclidTotal))
	if (np.sqrt(euclidTotal) < 1e-7):
		print ("MLE has converged!\nIterations: ", n)
		break
	elif (n > 50000):
		print ("MLE has not converged in 3000 iterations. Exiting.")
		break
	else:
		mle_n = np.copy(mle_n1)



with open(SOURCERUN_NAME + "Convergence.csv","w+", newline='') as file:
	
	writer=csv.writer(file,delimiter=',')
	for a in euclidArray:
		writer.writerow([a])
		
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
	print ("Z: ", (5*z-110), " Intensity: ",np.max(totalArr[z]))

#print (np.array(totalArr[][])

fig1, (ax1,ax2) = plt.subplots(1,2,constrained_layout=True,sharey=True)

zContour = np.arange(-zMax,zMax+1,5)
thetaContour = np.arange(0,360,3)

thetaCONT, zCONT = np.meshgrid(thetaContour,zContour)

#levels = [1e-3,5e-3,1e-2,2e-2,3e-2]
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


'''

#################### This is an attempt at generating MLEM for whole surface plane ###########################
######### This will most likely not work because we have only generated Transmission Matrix for a constant theta ##########



##### Initialize all arrays that will be used (data type is 64 bit floating point)
mle_n = np.array([[1]*noOfPixels], dtype=np.float64).T
mle_n1 = np.array([[1]*noOfPixels],dtype=np.float64).T
data = np.array(data, dtype=np.float64).T
transmissionMatrix = np.array(transmissionMatrix, dtype=np.float64)
background = np.array(background, dtype=np.float64).T
oneDim = np.array([1]*P).T

##### data and background are 2d arrays that are empty except for first index so extract these 1D arrays
data = data[0]
background = background[0]

##### N dimensional array of ones
identity = np.ones(N)


#### This loop runs the MLEM algorithm until convergence

n = 0
while True:
	condition=False
	n += 1
	

	##### Calculate the mle_(n+1) array from the mle_(n) array
	##### Follow the Presciption of Kowash Pg. 47

	#innerTerm = np.divide(data,(np.dot(transmissionMatrix,mle_n)+background))
	innerTerm = np.divide(data,np.dot(transmissionMatrix,mle_n))
	a = np.dot(transmissionMatrix,identity)
	
	new = []
	factor = []
	innerWithT = np.dot(transmissionMatrix.T, innerTerm)
	
	for q in innerWithT:
		new=np.divide(q,a)
		for r in new:
			factor.append(r)

	factor = np.array(factor)
	temp = np.multiply(mle_n[:,0], factor)
	for t in range(len(mle_n1)):
		value = temp[t]
		mle_n1[t][0]= value

	

	###### Check the convergence of mle_(n+1)
	###### If mle_(n+1) has not converged, then assign mle_(n+1) to mle_n and loop again
	for i in range(len(mle_n1)):
		if (abs(mle_n1[i][0] - mle_n[i][0]) <= np.average(mle_n1[:,0])*0.02):
			condition=True
		else:
			condition=False
			continue
	if (condition==True):
		print ("MLE has converged!")
		break
	elif (n > 100):
		print ("MLE has not converged in 100 iterations. Exiting.")
		break
	else:
		mle_n = np.copy(mle_n1)
		


##### Split mle_(n+1), which is a (N*P) X 1 dimensional array into and N X P 2 dimensional array so that we can plot with matplotlib.imshow()
zMax = 120
zStep = 5
noOfZSteps = 25
noOfAngleSteps = 120
num = 0
angleArr = []
totalArr = []
for x in range(noOfZSteps):
	angleArr = []
	for y in range(noOfAngleSteps):
		angleArr.append(mle_n1[num][0])
		num += 1
	totalArr.append(angleArr)



plt.imshow(totalArr, cmap=plt.cm.jet, origin='lower', extent=[0,6.28,-zMax,zMax], aspect='auto')
plt.xlabel("Theta (Radians)")
plt.ylabel("Z (cm)")
plt.show()

'''
