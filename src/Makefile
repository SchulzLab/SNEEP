
#sources
SRC =  	Matrix_new.cpp \
	differentialBindingAffinity.cpp \
	determineSignificantHitsWithoutDifferentialBindingAffinity.cpp \
	pvalue_copy.cpp  \
	callBashCommand.cpp \
	HandleInOutput.cpp \
	differentialBindingAffinity_multipleSNPs.cpp \
	differentialBindingAffinity_multipleSNPsSEMs.cpp \
	sampleRandomRsIDs.cpp

platform=$(shell uname -s)

ifeq ($(platform), Darwin)

CC = clang++
CFLAGS =  -std=c++11  -Xpreprocessor -fopenmp
#CFLAGS =  -std=c++11  -stdlib=libc++ 
LDFLAGS = -stdlib=libc++
all: $(SRC) Makefile
	$(CC) Matrix_new.cpp $(CFLAGS) $(LDFLAGS) -o Marix
	$(CC) differentialBindingAffinity.cpp $(CFLAGS) $(LDFLAGS) -o differentialBindingAffinity
	$(CC) determineSignificantHitsWithoutDifferentialBindingAffinity.cpp $(CFLAGS) $(LDFLAGS) -o determineSignificantHitsWithoutDifferentialBindingAffinity
	$(CC) pvalue_copy.cpp $(CFLAGS) $(LDFLAGS) -o pvalue
	$(CC) callBashCommand.cpp $(CFLAGS) $(LDFLAGS) -o callBashCommand
	$(CC) HandleInOutput.cpp $(CFLAGS) $(LDFLAGS) -o HandleInOutput
	$(CC) differentialBindingAffinity_multipleSNPs.cpp $(CFLAGS) $(LDFLAGS) -o differentialBindingAffinity_multipleSNPs
	$(CC) differentialBindingAffinity_multipleSNPsSEMs.cpp $(CFLAGS) $(LDFLAGS) -o differentialBindingAffinity_multipleSNPsSEMs
	$(CC) sampleRandomRsIDs2.cpp $(CFLAGS) $(LDFLAGS) -o randomRsIds_2.0
	
else

CC = g++ 
CFLAGS = -std=c++11 -fopenmp 
CFLAGS_ = -std=c++11 -fopenmp 
LDFLAGS = -static
LFLAG = -I/TL/opt/include
lFLAG = -lcurl
IFLAG = -I json/include/ # bedtools2/src/utils/driver/ bedtools2/src/utils/version/

#for bedtools (copy pased form the bedtools 2 makefile)
#SRC_DIR= bedtools2/src
#HTSDIR = $(SRC_DIR)/utils/htslib


#SUBDIRS = $(SRC_DIR)/annotateBed \
		  $(SRC_DIR)/bamToBed \
		  $(SRC_DIR)/bamToFastq \
		  $(SRC_DIR)/bedToBam \
		  $(SRC_DIR)/bedpeToBam \
		  $(SRC_DIR)/bedToIgv \
		  $(SRC_DIR)/bed12ToBed6 \
		  $(SRC_DIR)/closestFile \
		  $(SRC_DIR)/clusterBed \
		  $(SRC_DIR)/complementFile \
		  $(SRC_DIR)/coverageFile \
		  $(SRC_DIR)/expand \
		  $(SRC_DIR)/fastaFromBed \
		  $(SRC_DIR)/flankBed \
		  $(SRC_DIR)/genomeCoverageBed \
		  $(SRC_DIR)/getOverlap \
		  $(SRC_DIR)/groupBy \
		  $(SRC_DIR)/intersectFile \
		  $(SRC_DIR)/fisher \
		  $(SRC_DIR)/jaccard \
		  $(SRC_DIR)/linksBed \
		  $(SRC_DIR)/maskFastaFromBed \
		  $(SRC_DIR)/mapFile \
		  $(SRC_DIR)/mergeFile \
		  $(SRC_DIR)/multiBamCov \
		  $(SRC_DIR)/multiIntersectBed \
		  $(SRC_DIR)/nucBed \
		  $(SRC_DIR)/pairToBed \
		  $(SRC_DIR)/pairToPair \
		  $(SRC_DIR)/randomBed \
		  $(SRC_DIR)/regressTest \
		  $(SRC_DIR)/reldist \
		  $(SRC_DIR)/sampleFile \
		  $(SRC_DIR)/shiftBed \
		  $(SRC_DIR)/shuffleBed \
		  $(SRC_DIR)/slopBed \
		  $(SRC_DIR)/sortBed \
		  $(SRC_DIR)/spacingFile \
		  $(SRC_DIR)/split \
		  $(SRC_DIR)/subtractFile \
		  $(SRC_DIR)/summaryFile \
		  $(SRC_DIR)/tagBam \
		  $(SRC_DIR)/unionBedGraphs \
		  $(SRC_DIR)/windowBed \
		  $(SRC_DIR)/windowMaker

#UTIL_SUBDIRS =	$(SRC_DIR)/utils/FileRecordTools \
				$(SRC_DIR)/utils/FileRecordTools/FileReaders \
				$(SRC_DIR)/utils/FileRecordTools/Records \
				$(SRC_DIR)/utils/bedFile \
				$(SRC_DIR)/utils/BinTree \
				$(SRC_DIR)/utils/version \
				$(SRC_DIR)/utils/bedGraphFile \
				$(SRC_DIR)/utils/chromsweep \
				$(SRC_DIR)/utils/Contexts \
				$(SRC_DIR)/utils/general \
				$(SRC_DIR)/utils/gzstream \
				$(SRC_DIR)/utils/fileType \
				$(SRC_DIR)/utils/bedFilePE \
				$(SRC_DIR)/utils/KeyListOps \
				$(SRC_DIR)/utils/NewChromsweep \
				$(SRC_DIR)/utils/sequenceUtilities \
				$(SRC_DIR)/utils/tabFile \
				$(SRC_DIR)/utils/BamTools-Ancillary \
				$(SRC_DIR)/utils/BlockedIntervals \
				$(SRC_DIR)/utils/Fasta \
				$(SRC_DIR)/utils/VectorOps \
				$(SRC_DIR)/utils/GenomeFile \
				$(SRC_DIR)/utils/RecordOutputMgr \
				$(SRC_DIR)/utils/ToolBase \
				$(SRC_DIR)/utils/driver

#INCLUDES_BEDTOOLS =		$(addprefix -I,$(SUBDIRS) $(UTIL_SUBDIRS)) \
				-I$(SRC_DIR)/utils/BamTools/include \
				-I$(HTSDIR) \
				-I$(SRC_DIR)/utils/lineFileUtilities \
				-I$(SRC_DIR)/utils/Point \
				-I$(SRC_DIR)/utils/stringUtilities \
				-lz



all: $(SRC) Makefile
	$(CC) Matrix_new.cpp $(CFLAGS) $(LDFLAGS) -o Matrix
	$(CC) differentialBindingAffinity.cpp $(CFLAGS) $(LDFLAGS) -o differentialBindingAffinity
	$(CC) determineSignificantHitsWithoutDifferentialBindingAffinity.cpp $(CFLAGS) $(LDFLAGS) -o determineSignificantHitsWithoutDifferentialBindingAffinity
	$(CC) pvalue_copy.cpp $(CFLAGS) $(LDFLAGS) -o pvalue
	$(CC) callBashCommand.cpp $(CFLAGS) $(LDFLAGS) -o callBashCommand
	$(CC) HandleInOutput.cpp $(CFLAGS) $(LDFLAGS) -o HandleInOutput
	#$(CC) buffering.cpp $(CFLAGS) $(LDFLAGS) -o buffering $(INCLUDES_BEDTOOLS)
	#$(CC) sampleRandomRsIDs.cpp $(CFLAGS_)     -o randomRsIds $(LFLAG) $(lFLAG) $(IFLAG) $(INCLUDES_BEDTOOLS)
	$(CC) sampleRandomRsIDs2.cpp $(CFLAGS) -o randomRsIds_2.0 
#	$(CC) differentialBindingAffinity_multipleSNPs.cpp $(CFLAGS) -o differentialBindingAffinity_multipleSNPs $(LFLAG) $(lFLAG) $(IFLAG) $(INCLUDES_BEDTOOLS)
	$(CC) differentialBindingAffinity_multipleSNPs.cpp $(CFLAGS_) -o differentialBindingAffinity_multipleSNPs
#	$(CC) differentialBindingAffinity_multipleSNPsSEMs.cpp $(CFLAGS) $(LDFLAGS) -o differentialBindingAffinity_multipleSNPsSEMs
	$(CC) getSNPInfo.cpp $(CFLAGS_)  -o getSNPInfo $(LFLAG) $(lFLAG) $(IFLAG)

endif

clean:
	-/bin/rm -f core *.o *~
