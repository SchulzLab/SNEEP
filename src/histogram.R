library(ggplot2)

#input files
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1] # input file
outputDir <- args[2] # output file

data <- read.delim(input_file, header = TRUE, sep = "\t") #read data

ggplot(data, aes(x=MAF)) + 
	geom_histogram(binwidth = 0.01, color="darkblue", fill="lightblue", boundary = 0.0) +
	labs(x = "Minor Allele Frequency (MAF)", y = "count") +
	theme_minimal() +
	theme(text = element_text(size = 18), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) 

ggsave(paste(outputDir, sep = ""))
