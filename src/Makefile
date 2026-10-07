# builds SNEEP (differentialBindingAffinity_multipleSNPs); the C++ code is header-only, i.e. one compilation unit
# getSNPInfo.cpp (preprocessing of the dbSNP file, see docs flag -k) is built separately if needed:
#	g++ -std=c++11 -O2 -fopenmp getSNPInfo.cpp -o getSNPInfo

HEADERS = Matrix_new.hpp pvalue_copy.hpp callBashCommand.hpp HandleInOutput.hpp sampleRandomRsIDs2.hpp stringUtils.hpp

platform=$(shell uname -s)

.PHONY: all clean

all: differentialBindingAffinity_multipleSNPs

ifeq ($(platform), Darwin)

CC = clang++
CFLAGS = -std=c++11 -O2 -Xpreprocessor -fopenmp -lomp
LDFLAGS = -stdlib=libc++ "-L/usr/local/opt/libomp/lib"
CPPFLAGS = "-I/usr/local/opt/libomp/include"

differentialBindingAffinity_multipleSNPs: differentialBindingAffinity_multipleSNPs.cpp $(HEADERS) Makefile
	$(CC) differentialBindingAffinity_multipleSNPs.cpp -lc $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o differentialBindingAffinity_multipleSNPs

else

CC = g++
CFLAGS = -std=c++11 -O2 -fopenmp

differentialBindingAffinity_multipleSNPs: differentialBindingAffinity_multipleSNPs.cpp $(HEADERS) Makefile
	$(CC) differentialBindingAffinity_multipleSNPs.cpp $(CFLAGS) -o differentialBindingAffinity_multipleSNPs

endif

clean:
	-/bin/rm -f core *.o *~ differentialBindingAffinity_multipleSNPs
