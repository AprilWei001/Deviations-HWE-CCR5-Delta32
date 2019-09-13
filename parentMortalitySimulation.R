rm(list = ls())
gc(reset=TRUE)
library(survival)
(.packages())

ratio = 1.2;
hazard <- read.delim("uk8082MaleFemaleForSimulation.txt", header=FALSE, sep="\t")
numParent = 225000;
sizeOfRep = 1000;

Results1 = rep(0,1,sizeOfRep);
Results2 = rep(0,1,sizeOfRep);
Results3 = rep(0,1,sizeOfRep);
Results4 = rep(0,1,sizeOfRep);
Results5 = rep(0,1,sizeOfRep);
Results6 = rep(0,1,sizeOfRep);
maf = 0.1159;

for (i in 1:sizeOfRep) {
  fathers = runif(numParent,0,1);
  mothers = runif(numParent,0,1);
  indF2 = (fathers <= maf^2); 
  indF1 = (fathers > maf^2) & (fathers <= maf*(2-maf));
  indF0 = (fathers > maf*(2-maf));
  fathers[indF1] = 1;
  fathers[indF0] = 0;
  fathers[indF2] = 2;
  indM2 = (mothers <= maf^2);
  indM1 = (mothers > maf^2) & (mothers <= maf*(2-maf));
  indM0 = (mothers > maf*(2-maf));
  mothers[indM1] = 1;
  mothers[indM0] = 0;
  mothers[indM2] = 2;
  hapFromFather = fathers * runif(numParent,0,1);
  hapFromFather[(hapFromFather <= 0.5) & (!indF2)] = 0;
  hapFromFather[(hapFromFather > 0.5 )] = 1;
  hapFromFather[indF2] = 1;
  hapFromMother = mothers * runif(numParent,0,1);
  hapFromMother[(hapFromMother <= 0.5) & (!indM2)] =0;
  hapFromMother[hapFromMother > 0.5] = 1;
  hapFromMother[indM2] = 1;
  children = hapFromFather+hapFromMother;
  childrenGeno = c(children,children);
  parentSex = c(rep(1,numParent),rep(0,numParent));
  motherAgeDeath = rep(0,numParent);
  fatherAgeDeath = rep(0,numParent);
  for (j in 1:100) {
    indSurMother = (motherAgeDeath >= j-1);
    indSurFather = (fatherAgeDeath >= j-1);
    survMother = runif(numParent,0,1);
    if (j == 1) {
      indSurM2 = (survMother > ((hazard[1,2] + hazard[2,2])*ratio)) & indSurMother & indM2;
      indSurMOther = (survMother > (hazard[1,2] + hazard[2,2])) & indSurMother & (! indM2);
    }
    else{
      indSurM2 = (survMother > hazard[j+1,2]*ratio) & indSurMother & indM2;
      indSurMOther = (survMother > hazard[j+1,2]) & indSurMother & (! indM2);
    }
    motherAgeDeath[indSurM2] = motherAgeDeath[indSurM2] + 1;
    motherAgeDeath[indSurMOther] = motherAgeDeath[indSurMOther] + 1;
    survFather =  runif(numParent,0,1);
    if (j == 1) { # correct for the fact that the first year survival is seperated in to the first week after birth and 0-1 year.
      indSurF2 = (survFather > ((hazard[1,1] + hazard[2,1])*ratio)) & indSurFather & indF2;
      indSurFOther = (survMother > (hazard[1,1] + hazard[2,1])) & indSurFather & (!indF2);
    }
    else {
      indSurF2 = (survFather > hazard[j+1,1] * ratio) & indSurFather & indF2;
      indSurFOther = (survFather > hazard[j+1,1]) & indSurFather & (!indF2);
    }
    fatherAgeDeath[indSurF2] = fatherAgeDeath[indSurF2] + 1;
    fatherAgeDeath[indSurFOther] = fatherAgeDeath[indSurFOther] + 1;
  }
  parentDeathAge = c(fatherAgeDeath,motherAgeDeath);
  indNotDie = (parentDeathAge == j);
  parentDeathAge[indNotDie] = parentDeathAge[indNotDie] + 1;
  parentDeathEvent = rep(1, 1, 2*numParent);
  coxResult <- coxph(Surv(parentDeathAge, parentDeathEvent)~ childrenGeno + parentSex);
  indDie =  (parentDeathAge < j);
  maf1 = sum(children)/numParent/2;
  coxResult2 <- coxph(Surv(parentDeathAge[indDie], parentDeathEvent[indDie])~ childrenGeno[indDie] + parentSex[indDie]);
  maf2 = sum(childrenGeno[indDie])/sum(parentDeathEvent[indDie])/2;
  Results1[i] = exp(coxResult$coefficients[1]/maf1 * 2) -1;
  Results2[i] = (exp(coxResult$coefficients[1])-1)/maf1 * 2;
  Results3[i] = (exp(2*coxResult$coefficients[1])-1)/maf1;
  Results4[i] = exp(coxResult2$coefficients[1]/maf2 * 2) -1;
  Results5[i] = (exp(coxResult2$coefficients[1])-1)/maf2 * 2;
  Results6[i] = (exp(2*coxResult2$coefficients[1])-1)/maf2;
}

sink("SimulationParentMortality20Percent.txt")
print("  exp(coxResult$coefficients[1]/maf1 * 2) -1")
print(Results1)
print(" (exp(coxResult$coefficients[1])-1)/maf1 * 2")
print(Results2)
print("(exp(2*coxResult$coefficients[1])-1)/maf1")
print(Results3)
print("  exp(coxResult2$coefficients[1]/maf2 * 2) -1")
print(Results4)
print("(exp(coxResult2$coefficients[1])-1)/maf2 * 2")
print(Results5)
print(" (exp(2*coxResult2$coefficients[1])-1)/maf2")
print(Results6)
sink()
