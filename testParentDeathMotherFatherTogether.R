rm(list = ls())
gc(reset=TRUE)
library(survival)
(.packages())

#Unfortunately, we cannot share files for Cox analysis, because it contains individual level data that needs UK Biobank approval, but we provide the output of this script

dataDelta32 <- read.delim("dataForCoxParentDeathNoKin.txt", header=TRUE, sep="\t")

dataDelta32 <- dataDelta32[dataDelta32$deathEven ==1,]
dataDelta32$Center <- factor(dataDelta32$Center)
testDelta32AddParentNoKin <- (coxph(Surv(dataDelta32$ageEnd,dataDelta32$deathEven)~ dataDelta32$Delta32Additive + dataDelta32$parent
                                        +dataDelta32$pc1+dataDelta32$pc2+dataDelta32$pc3+dataDelta32$pc4+dataDelta32$pc5 +dataDelta32$pc6+dataDelta32$pc7+dataDelta32$pc8+dataDelta32$pc9+dataDelta32$pc10
                                        +dataDelta32$pc11+dataDelta32$pc12+dataDelta32$pc13+dataDelta32$pc14+dataDelta32$pc15 +dataDelta32$pc16+dataDelta32$pc17+dataDelta32$pc18+dataDelta32$pc19+dataDelta32$pc20
                                        +dataDelta32$pc21+dataDelta32$pc22+dataDelta32$pc23+dataDelta32$pc24+dataDelta32$pc25 +dataDelta32$pc26+dataDelta32$pc27+dataDelta32$pc28+dataDelta32$pc29+dataDelta32$pc30
                                        +dataDelta32$pc31+dataDelta32$pc32+dataDelta32$pc33+dataDelta32$pc34+dataDelta32$pc35 +dataDelta32$pc36+dataDelta32$pc37+dataDelta32$pc38+dataDelta32$pc9+dataDelta32$pc40
                                    +dataDelta32$Center)) 

sink("Delta32ResultsCoxForParentDeath.txt")

print("dataForCoxParentDeathNoKin.txt")
print(summary(testDelta32AddParentNoKin))
sink()
