# Deviations-HWE-CCR5-Delta32
This repository contains the scripts and results in Wei and Nielsen, 2019 "Deviations from Hardy Weinberg Equilibrium at CCR5-∆32 in Large Sequencing Data Sets"

For analyzing gnomAD data:
filterGnomadVCFData.m, where its input files gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf and gnomad.exomes.r2.1.1.sites.vcf can be downloaded from the gnomAD website.

mergeAndAnalyzeGnomadWES_WGS.m will take the output files from filterGnomadVCFData.m as input files

outputBi_gnomad_individualPop.m take the output file from mergeAndAnalyzeGnomadWES_WGS.m as input, and its output files Bi_gnomad_individualPopControlSNPs.txt, and Bi_gnomad_individualPopDelta32.txt are also provided. 



For Cox analysis:
It uses individual level UK Biobank data thus we cannot provide the input and output.

writeFileForCoxModelRemoveKinMaxPower.m This script will remove kinship while maximazing the number of minor allele homozygotes. 

testDelta32AndTagSNPCox.R and its output Delta32AndTagSNPCoxResult.txt are provided

testParentDeathMotherFatherTogether.R and its output Delta32ResultsCoxForParentDeath.txt and its output are provided



For simulation:
parentMortalitySimulation.R and its output file SimulationParentMortality20Percent.txt are provided. 


For figures:
plotFigure1.m and its input files ukbArrayControlSNPsCounts.txt, ukbWESControlSNPsCounts.txt, Bi_gnomad_individualPopControlSNPs.txt are provided

plotFigure2.m
