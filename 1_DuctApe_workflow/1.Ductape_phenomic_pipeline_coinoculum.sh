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
# Title of the script: DuctApe analysis workflow			      #
# prepared as supporting material for the workshop "EXCALIBUR Training Series:#
# Addressing microbial metabolic profile by means of Phenotype Microarray     #
# technology (BIOLOG)                                                         #
#                                                                             #
# Purpose of the script: Demonstrate the use of DuctApe to analyze phenotype  #
# microarray Phenomic_data						      #
###############################################################################

#!/bin/bash

# Remember to activate the ductape conda environment before runnin the script with
# conda activate ductape

# Set the strain to use as comparison in ring plot
comparisonStrain=CoInoculum 

# Setup the project
dape init

# Add the different strain. Name of the strain must match the name
# of the phenomic file from BIOLOG plate in the "Phenomic_data" folder
dape add StrainA -c red
dape add StrainB -c blue
dape add CoInoculum -c violet

# Add phenomic data
dphenome add-dir ./Phenomic_data

# Negative control well substraction
dphenome zero

# Start main phenome module
dphenome start -g
# -g option is used in this tutorial to avoid Kegg map fetching, remove this 
# option if metabolic network or mapping has to be performed by DuctApe

# Plot results
# Adding option -s produce svg plots
dphenome plot
dphenome rings
mv ActivityRing.png ActivityRingNormal.png
# Rename file or get overwrited in diff mode
dphenome rings -o $comparisonStrain
mv ActivityRing.png ActivityRingDiff.png

# Save other results
dphenome export
dphenome stats
