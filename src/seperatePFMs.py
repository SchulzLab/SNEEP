import sys, os,random, math

EPSILON = 0.001 
ZERO2 = [0.001, 0.002]
ZERO3 = [0.001, 0.002, 0.003]

#ENTROPY = 2.0
ENTROPY = 1.9
ENTROPY_2 = 1.8

#EPSILON = 0.0001 
#ZERO2 = [0.0001, 0.0002]
#ZERO3 = [0.0001, 0.0002, 0.0003]
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
	
		c_entropy = 0

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
				if (entropy >= ENTROPY_2):
					c_entropy += 1


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
		if (motifLength - len(r)) < 7:
			print(TF)
			A_ = A
			C_ = C
			G_ = G
			T_ = T
			# adapte c_entropy, add also bases with entropy higher than 1.9
			c_entropy += len(r)

		else:	

		#print("remove: " + str(r))
		#if len(r) >= 1:
		#	print("jap")
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
			o.write(TF + '\t' + str(len(A)) + '\t' + str(len(r)) + '\t' +  str(start) + '\t' + str(end) + '\t' + str(len(A_)) + '\t' + str(c_entropy) + '\n')

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

if (len(sys.argv) < 4):
	print("python3 seperatePFMsAndCheckActivity.py TransfacPFMFile, outputPFMDir, infoFile")
else:
	TRANSFACFile = sys.argv[1]
	outputDir = sys.argv[2]
	infoFile = sys.argv[3]

	random.seed(123)

	with open(infoFile, 'w') as o:
		o.write("TF\toriginalLen\tbasesWithHighEntropy\tstart\tend\tadaptedLen\t#basesEntropyHigher1.8\n")
	A = []
	C = []
	G = []
	T = []
	id_ = False
#	active = activeTFs.keys()
	TF_filename = ""
	with open(TRANSFACFile, 'r') as input_:
		for line in input_:
			line = line.strip()
			if line[0] == "I": #for JASPAR2020
				TF = ""
				TF2 = ""
				line = line.split()
				TF_filename = line[1]
				#print(TF_filename)
				#need to distinguish TFs since there are variants and combined binding sites
				helper = line[1]
				pos = helper.find("(") 
				pos2 = helper.find("::")
				if pos !=  -1:
					if pos2 != -1: 
						helper = helper.split("::")
						TF = helper[0]
						TF2 = helper[1]
					else:
						TF = helper[:pos]
				elif pos2 != -1:
					helper = helper.split("::")
					TF = helper[0]
					TF2 = helper[1]
				else: 
					TF = helper
			elif line[0] == "P" or line[0] == "A" or line[0] == "C" or line[0] =="X" or line[0] == "D" or line[0] == "N" or line[0] == "B":
				continue
			elif line[0] == "/":
				A,C,G,T = detemineEntropy(A,C,G,T, infoFile, TF_filename)
				convertCountToFreq(A,C,G,T)
			#	check(A,C,G,T)
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
				#print(line)
				A.append(round(float(line[1])))
				C.append(round(float(line[2])))
				G.append(round(float(line[3])))
				T.append(round(float(line[4])))
