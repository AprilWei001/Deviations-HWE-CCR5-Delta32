clear all
close all
clc

%Unfortunately, all the files used here contains individual level data, thus we cannot share them. The final output result from survival analysis is available to download
PCA= importdata('../../../../Biobank-data/PCAresultWBFiveAgeGroup.fam');
pcaID = PCA(:,1);
PCA = PCA(:,4:end);

%add code to remove the non-aa individuals in the pair, unless both are aa.
AgeFile = importdata('../AgeAndAncestryInformation_version2.txt');
FID = AgeFile(:,1);

[diffID indDiff] = setdiff(FID,pcaID);
AgeFile(indDiff,:)=[];
FID(indDiff)= [];
British = AgeFile(:,3);
indBritish = find(British == 1);


Sex = importdata('../../../../Biobank-data/whiteBritishMaleFemale.fam');
Sex = Sex(:,3)-1;
Sex(indDiff) = [];

Center = importdata('../../../../Biobank-data/CenterID.fam');
Center = Center(:,3);
Center(indDiff) = [];

daysAtRecruit = AgeFile(:,4);
indBirth = find(daysAtRecruit>0);
AgeAtRecruit = AgeFile(:,5);

%-----deal with the approximation of using 15th of month as bd------------
AgeFile0 = importdata('../AgeAndAncestryInformation.txt');
AgeFile0(indDiff,:)= [];
AgeRecruitOriginal = AgeFile0(:,5);
ind = find(AgeAtRecruit < AgeRecruitOriginal);
AgeAtRecruit(ind) = AgeRecruitOriginal(ind);

indProbmatic = find(AgeAtRecruit > AgeRecruitOriginal+2); %17 are problematic
%These are the ones should be removed perhaps, because age exit is also
%inaccurate
%among them one dead, and two are British
%the two British did not die

ind = find(AgeAtRecruit > AgeRecruitOriginal+1);
AgeAtRecruit(ind) = AgeRecruitOriginal(ind)+1;
%----deal with the approximation of using 15th of month as bd------------

indAgeRecruit = find(AgeAtRecruit>0);

AgeExit = (datenum(2016,02,16) - daysAtRecruit)/365.25 + AgeAtRecruit;
AgeAtDeath = AgeFile(:,6);

clear AgeFile

%-----------------Here to switch imput output genotype file ----------
 fid = fopen('CCRgenesVCF.txt');
 rowDelta32 = 12;
 for rowpass = 1:rowDelta32+7;
     tline = fgetl(fid);
 end
 fidOut = fopen('dataForCoxWithSexCenterAndRemoveKin3Delta32.txt','w');


%fid = fopen('del32Linkrs113010081VCF.txt');
%rowThis = 1;
%for rowpass = 1:rowThis+7;
%     tline = fgetl(fid);
%end
%fidOut = fopen('dataForCoxWithSexCenterAndRemoveKin3rs113010081.txt','w');
%------------------------------------------

outGeno = strsplit(tline);
outGeno(1:9)
outGeno = outGeno(10:end);
outGeno(indDiff) = [];
Geno = zeros(length(outGeno),1);
indaa = find(strcmp(outGeno,'1/1')==1);
Geno(indaa) = 1;
%indAa = find(strcmp(outGeno,'0/1')~=1);

indNotMissing = find(strcmp(outGeno,'./.')~=1);
fclose(fid); 


indBritishRecruit = intersect(indBritish,indBirth);
indBritishRecruit = setdiff(indBritishRecruit,intersect(indProbmatic,indBritishRecruit));
indBritishRecruit = intersect(indBritishRecruit,indNotMissing);
outGeno = outGeno(indBritishRecruit);
Geno = Geno(indBritishRecruit);
AgeAtRecruit = AgeAtRecruit(indBritishRecruit);
AgeExit = AgeExit(indBritishRecruit);
AgeAtDeath = AgeAtDeath(indBritishRecruit);
row = length(indBritishRecruit);
Event = zeros(row,1);
indDeath = find(AgeAtDeath>0);
Event(indDeath) = 1;
AgeExit(indDeath) = AgeAtDeath(indDeath);
PCA = PCA(indBritishRecruit,:);
Sex = Sex(indBritishRecruit);
FID = FID(indBritishRecruit);
Center = Center(indBritishRecruit);

%----------Control Relatives----------------
Rel = importdata('ukb3367_rel_s488364.dat');
Rel = Rel.data;

uniqueRelID = unique([Rel(:,1);Rel(:,2)]);
countUnique = zeros(length(uniqueRelID),1);
combineRel = [Rel(:,1);Rel(:,2)];
parfor i = 1:length(uniqueRelID);
    countUnique(i) = length(find(combineRel == uniqueRelID(i)));
end

FIDExcess = uniqueRelID(find(countUnique>10));
indKing = find(Rel(:,end)> 0);%0.0884) %0.0884 if removing up to second degree;
RelID = Rel(indKing,1:2);
removeDueToKin = [];

for i = 1:length(indKing);
    ind1 = find(FID == RelID(i,1));
    ind2 = find(FID == RelID(i,2));
    indEx1 = find(FIDExcess == RelID(i,1));
    indEx2 = find(FIDExcess == RelID(i,2));
    if length(ind1)*length(ind2) == 1;
        if length(indEx1)* length(indEx2) == 0;
            if Geno(ind1) == 1;
                removeDueToKin = [removeDueToKin;ind2];
            else
                removeDueToKin = [removeDueToKin;ind1];
            end
        end
    end
end
length(removeDueToKin)

for i = 1:length(FIDExcess);
    ind = find(FID == FIDExcess(i));
    if length(ind) == 1;
        removeDueToKin = [removeDueToKin;ind];
    end
end


Result = [FID,Sex,Center,AgeAtRecruit,AgeExit,Event,Geno,PCA];
Result(removeDueToKin,:) = [];
AgeExit(removeDueToKin,:) = [];
AgeAtRecruit(removeDueToKin,:) = [];
Event(removeDueToKin) = [];
outGeno(removeDueToKin) = [];

ind = find(AgeExit <= AgeAtRecruit);
Result(ind,:) = [];
Event(ind) = [];
outGeno(ind) = [];


fprintf(fidOut,'FID\tSex\tCenter\tageStart\tageEnd\tdeathEven\tSNP\t');

for i = 1:40;
    fprintf(fidOut,strcat('pc',num2str(i),'\t'));
end
fprintf(fidOut,'\n');

row = size(Result,1);

for i = 1:row;
    fprintf(fidOut,'%d\t',Result(i,:));
    fprintf(fidOut,'\n');
end
fclose(fidOut);

