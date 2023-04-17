###############################################################################
# Author: Francesco Vitali                                                    #
# Affiliation: Research Centre for Agriculture and Environment, Council       #
# for Agricultural Research and Economics (CREA-AA)                           #
# e-mail: francesco.vitali@crea.gov.it                                        #
# twitter: @svito_fi                                                          #
# github: FrancescoVit                                                        #
# Date: 15/04/2023                                                            #
#                                                                             #
# Licence: MIT licence                                                        #
# Copyright (c) 2023 Excalibur H2020 Project                                  #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a copy#
# of this software and associated documentation files (the "Software"), to    #
# deal in the Software without restriction, including without limitation the  #
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or #
# sell copies of the Software, and to permit persons to whom the Software is  #
# furnished to do so, subject to the following conditions:                    #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     # 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
###############################################################################
###############################################################################
## Script name: Biolog export converter 
##
## Purpose of script: Takes in multiple export from biolog software (in the one 
## line header format) and join them in an unique csv file that is compatible 
## with Ductape. In this script, all plates data are reduced to a desired 
## experiment length in hours, also knowing desired read point per hours (2 or 4)
##
## TO DO: include the possibility to let user decide at which time/hour to crop 
## all run. Eventually if user do not know, include possibility to crop to 
## autmatically crop to the shortest. Also insert possibility to no crop
##
################################################################################

# Set User parameters
analysisDir <- "/Path/to/folder"  
# Change it accordingly to the full path of the folder containing the exported files

# Get the file in the folder
filesList <- list.files(path= analysisDir, 
                        pattern="*.csv", 
                        full.names=TRUE, 
                        recursive=FALSE, 
                        all.files = F)

# Set variables
ExpLength = 73 # length in hours of the experiment, choose the minimum common to all plates
MeasurePerHour = 2 # number of measurement points per hour

# load libraries
library(tidyverse)

#
# Main loop for reading file and creating a dataframe for sample sheet
#

# build output df

df_out = data.frame()
df_newline = data.frame(matrix(nrow = 2, ncol = 97))

for (i in 1:length(filesList)) {
  
  biolog_file <- filesList[i] # get filenames
  modern_export = read.table(file = biolog_file, sep = ";", header = F) 
  
  # adjust reading time across different plates
  data = modern_export[1:nrow(modern_export),-c(1:8)]
  row_to_keep = data[1,] 
  
  data %>% 
    filter(!grepl(x = data$V9, ".25|.75")) -> data 
  # here the filtering to 2 measurements per hours is hardcoded by removing
  # occurences at .25 and .75 hours, to do here use the MeasurePerHour variable 
  # to automatically filter
  
  # extract data till desired duration in hours
  data = data[1:(ExpLength*MeasurePerHour+2),]
  
  # create a df for each loop to store template
  df_loop = data.frame(matrix(nrow = nrow(data) + 10, ncol = 97)) 
  
  # put correct data in place
  df_loop[11:nrow(df_loop),] =  data
  df_loop[1,2] = modern_export[2,1]
  df_loop[2,2] = modern_export[2,2]
  df_loop[3,2] = modern_export[2,3]
  df_loop[4,2] = modern_export[2,4]
  df_loop[5,2] = modern_export[2,7]
  df_loop[6,2] = modern_export[2,6]
  df_loop[7,2] = modern_export[2,6]
  df_loop[8,2] = modern_export[2,6]
  df_loop[9,2] = modern_export[2,8]
  
  # set names in first column
  df_loop[1,1] <- "Data File"
  df_loop[2,1] <- "Set up Time"
  df_loop[3,1] <- "Position"
  df_loop[4,1] <- "Plate Type"
  df_loop[5,1] <- "Strain Type"
  df_loop[6,1] <- "Sample Number"
  df_loop[7,1] <- "Strain Name"
  df_loop[8,1] <- "Strain Number"
  df_loop[9,1] <- "Other"
  
  # append to output dataframe, adding two empty lines 
  # between header and plate data
  df_out = rbind(df_out, df_loop, df_newline)
  }

# remove NA
df_out[is.na(df_out)] <- ""

# write output to csv file
write.table(x = df_out, 
            file = "./Elaborated_ductape_input_phenome.csv", 
            # change here accordin to strain name, must be the same that will be 
            # used as organism name in the DuctApe pipeline (in command dape add)
            sep = ",",
            dec = ".",
            row.names=F,
            col.names = F)

