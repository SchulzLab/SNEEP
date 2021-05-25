
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
	$(CC) pvalue_copy.cpp $(CFLAGS) $(LDFLAGS) -o pvalue
	$(CC) callBashCommand.cpp $(CFLAGS) $(LDFLAGS) -o callBashCommand
	$(CC) HandleInOutput.cpp $(CFLAGS) $(LDFLAGS) -o HandleInOutput
	$(CC) differentialBindingAffinity_multipleSNPs.cpp $(CFLAGS) $(LDFLAGS) -o differentialBindingAffinity_multipleSNPs
	$(CC) sampleRandomRsIDs2.cpp $(CFLAGS) $(LDFLAGS) -o randomRsIds_2.0
	
else

CC = g++ 
CFLAGS = -std=c++11 -fopenmp 
CFLAGS_ = -std=c++11 -fopenmp 
LDFLAGS = -static
LFLAG = -I/TL/opt/include
lFLAG = -lcurl

all: $(SRC) Makefile
	$(CC) Matrix_new.cpp $(CFLAGS) $(LDFLAGS) -o Matrix
	$(CC) pvalue_copy.cpp $(CFLAGS) $(LDFLAGS) -o pvalue
	$(CC) callBashCommand.cpp $(CFLAGS) $(LDFLAGS) -o callBashCommand
	$(CC) HandleInOutput.cpp $(CFLAGS) $(LDFLAGS) -o HandleInOutput
	$(CC) sampleRandomRsIDs2.cpp $(CFLAGS) -o randomRsIds_2.0 
	$(CC) differentialBindingAffinity_multipleSNPs.cpp $(CFLAGS_) -o differentialBindingAffinity_multipleSNPs

endif

clean:
	-/bin/rm -f core *.o *~
