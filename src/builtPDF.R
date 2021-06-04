library(rmarkdown) #for RMarkdown
library(tinytex) #latex package
library(knitr)
#install.packages("kableExtra")
options(kableExtra.latex.load_packages = FALSE)
library(kableExtra)
#install.packages("bookdown")
library(bookdown) #for cross references within a document
#install.packages("gprofiler2")

#pandoc must be installed 
#tinytex::tlmgr_install('pdfcrop') #install pdfcrop hat Florian gemacht

#tinytex::install_tinytex() #muss man einmal ausführen
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1] # input file
numberRounds <-args[2] #number of randmly sampled background rounds
minTFCount <- args[3] #min number TF counts wchih needs to be exceeded for sig TF hits
sourceDir <- args[4]

result = paste0(input_dir, "/result.txt")
print(result)
TFcount = paste0(input_dir, "/TF_count.txt")
#TFcount = paste0(input_dir, "/sampling/TF_count.txt")
print(TFcount)
info = paste0(input_dir, "/info.txt" )
print(info)

#toc = F -> avoid showing table of content ( which is automatically included when using bookdown)
rmarkdown::render(paste0(sourceDir,'visualization.Rmd'), output_format = bookdown::pdf_document2(toc = F), params = list( rounds = numberRounds, resultFile = result, TFcountFile = TFcount, infoFile = info, minTFCount = minTFCount), output_file = paste0(input_dir, "/summaryReport.pdf")) 
