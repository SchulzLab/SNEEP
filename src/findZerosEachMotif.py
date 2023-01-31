import scipy.optimize as sp
#import RegscorePy as reg
import sys, os, random, math
import numpy as np
import warnings
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from  scipy import stats

from functools import partial

#suppress warnings
warnings.filterwarnings('ignore')

STEPSIZE = 0.05 # oder 0.1  or depending on motif length
ROUNDER = 2
PERCENTAGE = 0.25
VALUE = 0.01

def laplace_max_abs(x, l ,k):	

	result = k * (1-np.exp(-(abs(x))/l))**(k-1) * (1/l)*np.exp(-abs(x)/l)
	return result

#laplace_abs_max MLE
def laplace_max_abs_MLE(l,x,k):

	n = len(x)
	result = -n/l + ((1/l**2) * np.sum(abs(x))) - ( ((k-1)/l**2) * np.sum( abs(x) / (np.exp(abs(x)/l) -1) ))
	return result



def lapace_max_abs_cdf( x, l , k):
	result = (1 - np.exp(- abs(x) / l))**k
	return result


#laplace MLE 
def laplace_MLE(l,x):

	n = len(x)
	result = -n/l + (1/l**2) * np.sum(abs(x))

	return result

def MSE_func(l, obs, step_sizes, motifLength):  ## input: current scale estimated with newton, tail values histogram, number stepsizes, motifLength

	estimate_values = []
	for i in step_sizes: 
		estimate_values.append(laplace_max_abs(i,l, motifLength))

	MSE = (np.sum(np.array(obs) - np.array(estimate_values))**2) # normally sum is multipled by 1/numberSamples, but this is contant in our case
	return MSE

#def MSE_func(obs,pred):
#
#	return (np.sum(np.array(obs - pred))**2)

def readScores(f, motif, length,  counter):

	values = []
	with open(f, 'r') as i:
		i.readline() #skip header
		for line in i:
			line = line.strip().split('\t')
			values.append(abs(float(line[0])))
			if counter == 0:	
				print("here")
				motif = line[2]
				length = float(line[1]) * 2
			counter +=1
	print(motif)
	print(length)
	return np.array(values)
	#return values


if (len(sys.argv) < 4):
	print("python3 findZeros.py  maxDiffBind_length_x,  outputFile, outputFile distribution (pdf)")
else:
	motif = ""
	counter = 0
	maxScores = []
	motifLength = 0

	with open(sys.argv[1], 'r') as i:
		i.readline() #skip header
		for line in i:
			line = line.strip().split('\t')
			maxScores.append(abs(float(line[0])))
			if counter == 0:
				motif = line[2]
				motifLength = float(line[1]) * 2
				counter +=1

	counter = len(maxScores) # number of kmers
	print(counter)
	random.seed(13)

	## original
	sol2 = sp.newton( lambda l:laplace_max_abs_MLE(l, np.array(maxScores), motifLength), 0.001, maxiter=100)
	print("newton scale: " + str(sol2))

	## refine fitting by fitting scale to tail to do so compute residual /mean squared error 1/n sum(observed- predicted)^2 

	## Step 1: get values from histogram (bin size 0.001) and compute cutoff for 0.25% of the data

	## get counts and freq for stepsize 0.001
	step_size = list(np.arange(0, max(maxScores)+ STEPSIZE, STEPSIZE)) # get bin intervals (0 - STEPSIZE)(2*STEPSIZE - 3* STEPSIZE).... plus STEPSIZE, otherwise last bin is missing
	#print(step_size)
	counts_per_bin  = stats.binned_statistic(maxScores, [], statistic = 'count', bins = step_size)[0] # get counts per bin
	freq_all_bins = counts_per_bin / counter / STEPSIZE # compute frequeny per bin additional / STEPSIZE otherwise integral is not 1
	#print(list(counts_per_bin))

	## determine maxDiffBinding cutoff for X percentage of the data point -> right tail 
	sorted_max_scores = sorted(list(maxScores), reverse=True)
	#print(sorted_max_scores)
	## consider as right tail X% of the data points 
	tail_cutoff = round(sorted_max_scores[round(counter * PERCENTAGE)], ROUNDER)
	print("tailcutoff: " + str(tail_cutoff))
	number_bins_tail_cutoff = len(step_size) - int(tail_cutoff / STEPSIZE) # needs to substracttailcutoff from len of step-size otherwise left tail is considered
	print("number bins for right tail: " + str(number_bins_tail_cutoff))

	freq_per_bin = freq_all_bins[-number_bins_tail_cutoff:]
	#print(freq_per_bin)
	
	## Step 2: for scale evaluate laplace max abs in steps of 0.001
	estimate_values = []
	for i in step_size[-number_bins_tail_cutoff:]: 
		estimate_values.append(laplace_max_abs(i,sol2, motifLength))

	## Step 3: compute MSE for tail based on previously defined cutoff,
	#print(np.sum(np.array(freq_per_bin) - np.array(estimate_values))**2) # normally sum is multipled by 1/numberSamples, but this is contant in our case
	MSE_start = MSE_func(sol2,freq_per_bin, step_size[-number_bins_tail_cutoff:], motifLength)
	print("MSE for newton scale: " + str(MSE_start) + " scale: " + str(sol2))

	## Step 4: find local minimum 
	## do we need to increase or decrease the scale?
	increased_scale = MSE_func(sol2+VALUE,freq_per_bin, step_size[-number_bins_tail_cutoff:], motifLength)
	optimale_scale = 0
	previous = 0
	if (MSE_start > increased_scale): # increasing scale reduce error
	
		print("increase")
		previous = MSE_start
		print("previous: " + str(previous) + " increased: " + str(increased_scale))
		c = 1
		while (previous > increased_scale) == True:
			previous = increased_scale
			c += 1
			increased_scale = MSE_func(sol2+(c*VALUE),freq_per_bin, step_size[-number_bins_tail_cutoff:], motifLength)
			print("previous: " + str(previous) + " increased: " + str(increased_scale))
		optimale_scale = sol2 + ((c-1)*VALUE)
		print("best scale: " + str(optimale_scale) + ' c: ' + str(c-1) + " bestMSE: " + str(previous))

	else:
		decreased_scale = MSE_func(sol2-VALUE,freq_per_bin, step_size[-number_bins_tail_cutoff:], motifLength)
		if (MSE_start > decreased_scale): # decreasing scale reduce error
			print("decrease")
			previous = MSE_start
			print("previous: " + str(previous) + " decreased: " + str(decreased_scale))
			c = 1
			while (previous > decreased_scale) == True:
				previous = decreased_scale
				c += 1
				decreased_scale = MSE_func(sol2-(c*VALUE),freq_per_bin, step_size[-number_bins_tail_cutoff:], motifLength)
				print("previous: " + str(previous) + " decreased: " + str(decreased_scale))
			optimale_scale = sol2 - ((c-1)*VALUE)
			print("best scale: " + str(optimale_scale) + ' c: ' + str(c-1) + " bestMSE: " + str(previous))
	

		else: # scale is already optimale
			print("best already found")
			optimale_scale = sol2
			previous = MSE_start
			
	## Step 5: plot laplace max abs laplace for newton scale value and updated one 

	with open(sys.argv[2], 'a') as o:
		o.write(motif + '\t' + str(sol2) + '\t' + str(MSE_start) + '\t' + str(optimale_scale) + '\t' + str(previous) + '\t' +  str(motifLength/2) + '\t' + str(counter)+ "\n" )

	#plot histogram *empirical data) and laplace abs max distribution with computed scale value
	m= np.linspace(0,5, 100)
	helper = []
	helper2 = []
	for i in m:
		helper.append(laplace_max_abs(i,sol2, motifLength))
		helper2.append(laplace_max_abs(i,optimale_scale, motifLength))

	plt.plot(m, helper, 'm', label = "scale")
	plt.plot(m, helper2, 'g', label = "optimized scale")
	plt.hist(abs(np.array(maxScores)), density= True, bins = 100)


	plt.ylabel('Probability')
	plt.xlabel('maximal differential binding scores')
	plt.legend(loc="upper left")
	plt.title(motif + " "  + str(motifLength/2) + " optimized scale: " + str(round(optimale_scale, 2)) + " MSE: " + str(round(previous, 3)))

	plt.savefig( sys.argv[3] + "/" +  motif+  ".pdf")
