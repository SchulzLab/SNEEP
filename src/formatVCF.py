#!/usr/bin/env python

import sys, os

if (len(sys.argv) < 3):
    print("python3 formatVCF.py inputFile, outputFile")
else:
    inputFile = sys.argv[1]
    outputFile = sys.argv[2]

    with open(inputFile, 'r') as i, open(outputFile, 'w') as o:
        for line in i:
            if line[0] == "#":
                continue
            else:
                line = line.strip().split('\t')
                chr_ = line[0]
                pos = line[1]
                rsid = line[2]
                if rsid ==".": # in bed-like file not given info is given as - (not . as in vcf file)
                    rsid = "-"
                ref = line[3]
                alt = line[4]
                info = line[7]
                maf = "-1"
                stop = False
                if chr_ == "." or pos == "." or ref == "." or alt == ".": # important info missing -> skip this snp
                    continue
                elif "," in ref and "," in alt: # if for the reference and the alternative allele more than one option is provided -> skip the snp
                    continue
                else:
                    ## is maf given?
                    if "MAF" in info: 
                        info = info.split("MAF=")
                        maf = info[1]
                        #print(maf)
                        if ";" in maf:
                            maf = maf.split(';')[0]
                         #   print(maf)
                    ## is there more than one alternative allele given?
                    ref = ref.split(',')
                    ## is there more than one reference allele given?
                    alt = alt.split(',')

                    ## write output to new file
                    if len(ref) > 1: # more than one reference allele but only one alternative allele
                        for elem in ref:
                            o.write("chr" + chr_ + '\t' + str(int(pos)-1) + '\t' + pos + '\t' +  elem + '\t'  + alt[0] + '\t'+  rsid + '\t' + maf + '\n')

                    elif len(alt) > 1: # more than one reference allele but only one alternative allele
                        for elem in alt:
                            o.write("chr" + chr_ + '\t' + str(int(pos)-1) + '\t' + pos + '\t'+  ref[0] + '\t' + elem + '\t'+  rsid + '\t' + maf + '\n' )
                    else:
                            o.write("chr" + chr_ + '\t' + str(int(pos)-1) + '\t' + pos + '\t' +  ref[0] + '\t'  + alt[0] + '\t'+  rsid + '\t' + maf + '\n') 










