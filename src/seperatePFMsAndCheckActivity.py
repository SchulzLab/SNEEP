#!/usr/bin/env python

import sys, os, random, math

EPSILON = 0.001 
#EPSILON = 0.0001 
ZERO2 = [0.001, 0.002]
ZERO3 = [0.001, 0.002, 0.003]

ENTROPY = 1.9
#def compansatAddEpsilon(zero_A, A, no_zeros, i ,e):
#
#	if zero_A == False and A[i] > e and no_zeros > 0: 
#		A[i] = A[i] - e
#		no_zeros = no_zeros - 1
#	else:
#		zero_A = True
def getFreq(A, sum_):
	
	if A == 0:
		A = 1

	return(A/sum_)


def detemineEntropy(A,C,G,T, infoFile, TF):

		result = []

		#print("A: " + str(A))
		#print("C: " + str(C))
		#print("G: " + str(G))
		#print("T: " + str(T))

		#compute entropy
		motifLength = len(A)
		for counter in range(len(A)):
			sum_ = A[counter] + C[counter] + G[counter] + T[counter]
			a = getFreq(A[counter], sum_)
			c = getFreq(C[counter], sum_)
			g = getFreq(G[counter], sum_)
			t = getFreq(T[counter], sum_)
			entropy = -(a * math.log2(a)  + c * math.log2(c) + g * math.log2(g)  + t * math.log2(t))
			#print("entropy: " + str(entropy))
			if (entropy >= ENTROPY):
				result.append(True)
			else:
				result.append(False)


		r = [] # positions of the bases we want to remove
		counter = 0 # counter for the positions
		start = 0
		#only remove flanking bases 
		for elem in result:
			if elem == True: # and forward == True: # check if for the flanking positions the entriopt is >= the threshold
				r.append(counter) #add position to remove
				start += 1
			if elem == False:
				#forward = False #stop checking 
				break
			counter += 1
		
		counter = len(A)-1
		end = 0
		for elem in reversed(result):
			if elem == True: # and reverse == True: 
				r.append(counter)
				end += 1
			if elem == False:
				break
			counter -= 1

		A_, C_, G_, T_ = [],[],[],[]
		#print("remove: " + str(r))
		#if len(r) >= 1:
		#	print("jap")
		if (motifLength - len(r)) < 7: ## avoid to short motifs because the scale parameter is really bad for those
			#print(TF)
			A_ = A
			C_ = C
			G_ = G
			T_ = T

		else:	
			for i in range(len(A)):
				if not i in r:
					A_.append(A[i])
					C_.append(C[i])
					G_.append(G[i])
					T_.append(T[i])

		#print("A: " + str(A_))
		#print("C: " + str(C_))
		#print("G: " + str(G_))
		#print("T: " + str(T_))
		with open(infoFile, 'a') as o:
			o.write(TF + '\t' + str(len(A)) + '\t' + str(len(r)) + '\t' +  str(start) + '\t' + str(end) + '\t' + str(len(A_)) + '\n')

		return(A_,C_,G_,T_)



def convertCountToFreq(A,C,G,T):

	sum_ = 0.0
	number_of_zeros = 0
	zero_A = False
	zero_C = False
	zero_G = False
	zero_T = False
	multi_factor = 0 #how much epislon did we introduce?

	for i in range(len(A)):
		sum_ = A[i] + C[i] + G[i] + T[i]
#		print("sum: " + str(sum_))
		if A[i] != 0.0 and (A[i]/sum_) > EPSILON:	
			A[i] = A[i]/ sum_ 
		else:
			if A[i] == 0.0:
				A[i] = 0.0
			else:
				A[i] = A[i]/ sum_ 
			number_of_zeros = number_of_zeros + 1
			zero_A = True

		if C[i] != 0.0 and (C[i]/sum_) > EPSILON:
			C[i] = C[i]/ sum_ 
		else:
			if C[i] == 0.0:
				C[i] = 0.0
			else:
				C[i] = C[i]/ sum_ 
			number_of_zeros = number_of_zeros + 1
			zero_C = True

		if G[i] != 0.0 and (G[i]/sum_) > EPSILON:
			G[i] = G[i]/ sum_
		else:
			if G[i] == 0.0:
				G[i] = 0.0
			else:
				G[i] = G[i]/ sum_ 

			number_of_zeros = number_of_zeros + 1
			zero_G = True

		if T[i] != 0.0 and (T[i]/sum_) > EPSILON:
			T[i] = T[i]/ sum_
		else:
			if T[i] == 0.0:
				T[i] = 0.0
			else:
				T[i] = T[i]/ sum_ 

			number_of_zeros = number_of_zeros + 1
			zero_T = True
	
#		print(number_of_zeros)	 
		#add slightly varying epsilon to zeros to avoid 0 in the diffBind score
		if number_of_zeros == 1: # add epsilon to the position	
			multi_factor = EPSILON
			if zero_A == True:
				A[i] = A[i] + EPSILON
			if zero_C == True:
				C[i] = C[i]+ EPSILON
			if zero_G == True:
				G[i] = G[i] +  EPSILON
			if zero_T == True:
				T[i] = T[i] + EPSILON

		c = 0
		random.shuffle(ZERO2) #shuffle list random -> avoid assining same value to same letter 
		if number_of_zeros == 2: # add slightly variny epsioln to the positions _> sum is still 2* epsilon	
			multi_factor = 0.003
			if zero_A == True:
				A[i] = A[i] + ZERO2[c]
				c+= 1
			if zero_C == True:
				C[i] = C[i] + ZERO2[c]
				c+= 1
			if zero_G == True:
				G[i] = G[i] + ZERO2[c]
				c+= 1
			if zero_T == True:
				T[i] = T[i] + ZERO2[c]
				c+= 1
	
		c = 0
		random.shuffle(ZERO3)
		if number_of_zeros == 3: # add slightly varing epsilon to the position	-> sum is still 3* epsilon
			multi_factor = 0.006
			if zero_A == True:
				A[i] = A[i] + ZERO3[c]
				c+= 1
			if zero_C == True:
				C[i] = C[i] + ZERO3[c]
				c+= 1
			if zero_G == True:
				G[i] = G[i] + ZERO3[c]
				c+= 1
			if zero_T == True:
				T[i] = T[i] + ZERO3[c]
				c+= 1
		e = (multi_factor)/ (4-number_of_zeros)
		no_zeros = 4 - number_of_zeros
		while  no_zeros > 0:
			#if number_of_zeros == 3:
			#	e = (number_of_zeros * EPSILON)
			#else:
			if zero_A == False and (A[i]- e >=0.001) and no_zeros > 0: 
				A[i] = A[i] - e
				no_zeros = no_zeros - 1
			else:
				zero_A = True

			if zero_C == False and (C[i] - e >= 0.001)  and no_zeros> 0: 
				C[i] = C[i] - e	
				no_zeros = no_zeros - 1
			else: 
				zero_C = True
	
			if zero_G == False and (G[i] - e >= 0.001) and no_zeros > 0: 
				G[i] = G[i] - e
				no_zeros = no_zeros - 1
			else:
				zero_G = True

			if zero_T == False and (T[i] - e >= 0.001) and no_zeros > 0: 
				T[i] = T[i] - e		
				no_zeros = no_zeros - 1
			else:
				zero_T = True

		number_of_zeros = 0	
		zero_A = False
		zero_C = False
		zero_G = False
		zero_T = False
		multi_factor = 0 #how much epislon did we introduce?
#		print("A: " + str(A))
#		print("C: " + str(C))
#		print("G: " + str(G))
#		print("T: " + str(T))
	

def check(A,C,G,T):
	check = True
	for i in range(len(A)):
		if  A[i] +  C[i] +  G[i] +  T[i] != 1.0:
			check = False
			print(  "i: " + str(i)+ " " + str(A[i] + C[i] + G[i] +  T[i]))

def write(A):
	
	for elem in A[:-1]:
		output.write(str(elem) + '\t')
	output.write(str(A[-1]) + '\n')
	return

if (len(sys.argv) < 7):
	print("python3 seperatePFMsAndCheckActivity.py activeTFsFile, TransfacPFMFile, outputPFMDir, ENSG_GeneName, activityThreshold, infoFile")
else:
	activeTFsFile = sys.argv[1]
	TRANSFACFile = sys.argv[2]
	outputDir = sys.argv[3]
	ENSG_GeneFile = sys.argv[4]
	threshold = float(sys.argv[5])
	infoFile = sys.argv[6]

	random.seed(123)
	with open(infoFile, 'w') as o:
		o.write("TF\toriginalLen\tbasesWithHighEntropy\tstart\tend\tadaptedLen\n")
	#create mapping ensemble id -> gene name and via vice
	ensemble_name = {}
	name_ensemble = {}
	with open(ENSG_GeneFile, 'r') as input_:
		input_.readline()
		for line in input_:
			line = line.strip().split('\t')
			ensemble_name[line[0]] = line[1]
			name_ensemble[line[1]] = line[0]
		
	#store activity of all TFs
	ensembleIdsTFs = ensemble_name.keys()
	activeTFs = {}	
	with open(activeTFsFile, 'r') as input_:
		input_.readline()
		for line in input_:
			line = line.strip().split('\t')
			if line[0] in ensembleIdsTFs:
				activeTFs[line[0]] = float(line[1])
#	print(len(ensembleIdsTFs))
#	print(len(activeTFs.keys()))

	A = []
	C = []
	G = []
	T = []
	id_ = False
	active = activeTFs.keys()
	with open(TRANSFACFile, 'r') as input_:
		for line in input_:
			line = line.strip()
			#if line[0] == "D":
			if line[0] == "I":
				TF = ""
				TF2 = ""
				line = line.split(" ")
				#TF_filename = line[2]
				TF_filename = line[1]
				#need to distinguish TFs since there are variants and combined binding sites
				#helper = line[2]
				helper = line[1]
				pos = helper.find("(") 
				pos2 = helper.find("::")
				if pos !=  -1:
					if pos2 != -1: 
						helper = helper.split("::")
						TF = helper[0].strip()
						TF2 = helper[1].strip()
					else:
						TF = helper[:pos]
				elif pos2 != -1:
					helper = helper.split("::")
					TF = helper[0].strip()
					TF2 = helper[1].strip()
				else: 
					TF = helper.strip()
				#TODO nachsehen warum fast keien TFs ausgeschlossen werden	
				if TF in name_ensemble.keys() and name_ensemble[TF] in activeTFs and activeTFs[name_ensemble[TF]] > threshold:
					#check if it is a combined TF like ARNT::STAT2
					if TF2 != "":
						if TF2 in name_ensemble.keys() and name_ensemble[TF2] in activeTFs and activeTFs[name_ensemble[TF2]] > threshold:
							id_ = True
						if TF not in name_ensemble.keys():
							print("TF2  without ensemble id: " + TF2)
					else:
						id_ = True
				if TF not in name_ensemble.keys():
					print("TF  without ensemble id: " + TF)		
			if id_ == True:
				#if line[0] == "P" or line[0] == "D" or line[0] == "I" or line[0] == "A" or line[0] == 'X' or line[0] == 'C':
				if line[0][0] == "P" or line[0][0] == "D" or line[0][0] == "I" or line[0][0] == "A" or line[0][0] == 'X' or line[0][0] == 'C' or line[0][0] == "N" or line[0][0] == "B":
					continue
				elif line[0] == "/":
					id_ = False
					A,C,G,T = detemineEntropy(A,C,G,T, infoFile, TF_filename)
					convertCountToFreq(A,C,G,T)
	#				if check(A,C,G,T) == False:
	#					print ("AAAAAHHHH: " + name)
	#					print("A: " + str(A))
	#					print("C: " + str(C))
	#					print("G: " + str(G))
	#					print("T: " + str(T))
					with open(outputDir + '/' + TF_filename + ".txt", 'w') as output: 
						write(A)
						write(C)
						write(G)
						write(T)
					A = []
					C = []
					G = []
					T = []
				else: #read matrix
					line = line.split('\t')
					A.append(round(float(line[1]))) # necessary since jaspar sometime give float values ?!
					C.append(round(float(line[2])))
					G.append(round(float(line[3])))
					T.append(round(float(line[4])))
