import sys, os

EPSILON = 0.001 

#def compansatAddEpsilon(zero_A, A, no_zeros, i ,e):
#
#	if zero_A == False and A[i] > e and no_zeros > 0: 
#		A[i] = A[i] - e
#		no_zeros = no_zeros - 1
#	else:
#		zero_A = True


def convertCountToFreq(A,C,G,T):
	
	sum_ = 0.0
	number_of_zeros = 0
	zero_A = False
	zero_C = False
	zero_G = False
	zero_T = False

	for i in range(len(A)):
		sum_ = A[i] + C[i] + G[i] + T[i]
#		print("sum: " + str(sum_))
		if A[i] != 0.0 and (A[i]/sum_) > 2*EPSILON:	
			A[i] = A[i]/ sum_ 
		else:
			A[i] = EPSILON
			number_of_zeros = number_of_zeros + 1
			zero_A = True

		if C[i] != 0.0 and (C[i]/sum_) > 2*EPSILON:
			C[i] = C[i]/ sum_ 
		else:
			C[i] = EPSILON
			number_of_zeros = number_of_zeros + 1
			zero_C = True

		if G[i] != 0.0 and (G[i]/sum_) > 2*EPSILON:
			G[i] = G[i]/ sum_
		else:
			G[i] = EPSILON
			number_of_zeros = number_of_zeros + 1
			zero_G = True

		if T[i] != 0.0 and (T[i]/sum_) > 2*EPSILON:
			T[i] = T[i]/ sum_
		else:
			T[i] = EPSILON
			number_of_zeros = number_of_zeros + 1
			zero_T = True
	
#		print(number_of_zeros)	 
		e = (number_of_zeros * EPSILON)/ (4-number_of_zeros)
		no_zeros = 4 - number_of_zeros
		while  no_zeros > 0:
			#if number_of_zeros == 3:
			#	e = (number_of_zeros * EPSILON)
			#else:
			if zero_A == False and A[i] > e and no_zeros > 0: 
				A[i] = A[i] - e
				no_zeros = no_zeros - 1
			else:
				zero_A = True

			if zero_C == False and C[i] > e and no_zeros> 0: 
				C[i] = C[i] - e	
				no_zeros = no_zeros - 1
			else: 
				zero_C = True
	
			if zero_G == False and G[i] > e and no_zeros > 0: 
				G[i] = G[i] - e
				no_zeros = no_zeros - 1
			else:
				zero_G = True

			if zero_T == False and T[i] > e and no_zeros > 0: 
				T[i] = T[i] - e		
				no_zeros = no_zeros - 1
			else:
				zero_T = True

		number_of_zeros = 0	
		zero_A = False
		zero_C = False
		zero_G = False
		zero_T = False
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

if (len(sys.argv) < 6):
	print("python3 seperatePFMsAndCheckActivity.py activeTFsFile, TransfacPFMFile, outputPFMDir, ENSG_GeneName, activityThreshold")
else:
	activeTFsFile = sys.argv[1]
	TRANSFACFile = sys.argv[2]
	outputDir = sys.argv[3]
	ENSG_GeneFile = sys.argv[4]
	threshold = float(sys.argv[5])

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
					A.append(float(line[1]))
					C.append(float(line[2]))
					G.append(float(line[3]))
					T.append(float(line[4]))
