rm(list = ls())
gc(reset=TRUE)
library(survival)
(.packages())

#Unfortunately, we cannot share files for Cox analysis, because it contains individual level data that needs UK Biobank approval, but we provide the output of this script
dataDelta32 <- read.delim("dataForCoxWithSexCenterAndRemoveKin3Delta32.txt", header=TRUE, sep="\t")
dataDelta32$Center <- factor(dataDelta32$Center)

testDelta32<- (coxph(Surv(dataDelta32$ageStart,dataDelta32$ageEnd,dataDelta32$deathEven)~ dataDelta32$Sex + dataDelta32$SNP
                     +dataDelta32$pc1+dataDelta32$pc2+dataDelta32$pc3+dataDelta32$pc4+dataDelta32$pc5 +dataDelta32$pc6+dataDelta32$pc7+dataDelta32$pc8+dataDelta32$pc9+dataDelta32$pc10
                     +dataDelta32$pc11+dataDelta32$pc12+dataDelta32$pc13+dataDelta32$pc14+dataDelta32$pc15 +dataDelta32$pc16+dataDelta32$pc17+dataDelta32$pc18+dataDelta32$pc19+dataDelta32$pc20
                     +dataDelta32$pc21+dataDelta32$pc22+dataDelta32$pc23+dataDelta32$pc24+dataDelta32$pc25 +dataDelta32$pc26+dataDelta32$pc27+dataDelta32$pc28+dataDelta32$pc29+dataDelta32$pc30
                     +dataDelta32$pc31+dataDelta32$pc32+dataDelta32$pc33+dataDelta32$pc34+dataDelta32$pc35 +dataDelta32$pc36+dataDelta32$pc37+dataDelta32$pc38+dataDelta32$pc9+dataDelta32$pc40
                     +dataDelta32$Center))

sink("Delta32AndTagSNPCoxResult.txt")
print("Delta32")
print((summary(testDelta32)))


dataRS113010081 <- read.delim("dataForCoxWithSexCenterAndRemoveKin3rs113010081.txt", header=TRUE, sep="\t")
dataRS113010081$Center <- factor(dataRS113010081$Center)

testRS113010081<- (coxph(Surv(dataRS113010081$ageStart,dataRS113010081$ageEnd,dataRS113010081$deathEven)~ dataRS113010081$Sex + dataRS113010081$SNP
                         +dataRS113010081$pc1+dataRS113010081$pc2+dataRS113010081$pc3+dataRS113010081$pc4+dataRS113010081$pc5 +dataRS113010081$pc6+dataRS113010081$pc7+dataRS113010081$pc8+dataRS113010081$pc9+dataRS113010081$pc10
                         +dataRS113010081$pc11+dataRS113010081$pc12+dataRS113010081$pc13+dataRS113010081$pc14+dataRS113010081$pc15 +dataRS113010081$pc16+dataRS113010081$pc17+dataRS113010081$pc18+dataRS113010081$pc19+dataRS113010081$pc20
                         +dataRS113010081$pc21+dataRS113010081$pc22+dataRS113010081$pc23+dataRS113010081$pc24+dataRS113010081$pc25 +dataRS113010081$pc26+dataRS113010081$pc27+dataRS113010081$pc28+dataRS113010081$pc29+dataRS113010081$pc30
                         +dataRS113010081$pc31+dataRS113010081$pc32+dataRS113010081$pc33+dataRS113010081$pc34+dataRS113010081$pc35 +dataRS113010081$pc36+dataRS113010081$pc37+dataRS113010081$pc38+dataRS113010081$pc9+dataRS113010081$pc40
                         +dataRS113010081$Center))

print("dataForCoxWithSexCenterAndRemoveKin3rs113010081")
print(summary(testRS113010081))

sink()
