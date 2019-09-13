clear all
close all
clc

WES = importdata('gnomad_WES_fin_swe_est_nwe_onf.txt');
indDelta32WES = find(contains(WES.textdata(:,4),'TACAGTCAGTATCAATTCTGGAAGAATTTCCAG'));
Delta32WES = WES.data(indDelta32WES,:);
ind = find(strcmp(WES.textdata(:,1), 'X')~=1);
rsNumWES = WES.textdata(ind,3);
WES = WES.data(ind,:);
WGS = importdata('gnomad_WGS_exome_fin_est_nwe_onf.txt');
indDelta32WGS = find(contains(WGS.textdata(:,4),'TACAGTCAGTATCAATTCTGGAAGAATTTCCAG'));
Delta32WGS = WGS.data(indDelta32WGS,:);
ind = find(strcmp(WGS.textdata(:,1), 'X')~=1);
rsNumWGS = WGS.textdata(ind,3);
WGS = WGS.data(ind,:);

[sharedRS indWES indWGS] = intersect(rsNumWES,rsNumWGS);
DataWES = WES(indWES,:);
DataWGS = WGS(indWGS,:);

DataMerged = DataWES(:,sort([[1:4:20],[2:4:20],[4:4:20]]));
DataMerged(:,1:3) = DataMerged(:,1:3)+DataWGS(:,[1,2,4]);
DataMerged(:,7:end) = DataMerged(:,7:end)+DataWGS(:,sort([[5:4:16],[6:4:16],[8:4:16]]));
DataMerged = [DataMerged,DataMerged(:,1:3)+DataMerged(:,4:6)+DataMerged(:,7:9)+DataMerged(:,10:12)+DataMerged(:,13:15)];
pPop = DataMerged(:,1:3:end)./DataMerged(:,2:3:end);

DataDelta32 = Delta32WES(:,sort([[1:4:20],[2:4:20],[4:4:20]]));
DataDelta32(:,1:3) = DataDelta32(:,1:3)+Delta32WGS(:,[1,2,4]);
DataDelta32(:,7:end) = DataDelta32(:,7:end) + Delta32WGS(:,sort([[5:4:16],[6:4:16],[8:4:16]]));
DataDelta32 = [DataDelta32,sum(DataDelta32(:,1:3:end)),sum(DataDelta32(:,2:3:end)),sum(DataDelta32(:,3:3:end))];
pDelta32 = DataDelta32(:,1:3:end)./DataDelta32(:,2:3:end);

ResultDelta32 = [];
CI95 = [];
pVal = [];
numSNPs = [];
SampleSize = [];
medianHWE = [];
fidout = fopen('Bi_gnomad_individualPopControlSNPs.txt','w');
for i = 1:6;
    sampleSize = DataDelta32(3*(i-1)+2)*0.9;
    ind = intersect(find(pPop(:,i)> pDelta32(i)-0.01),find(pPop(:,i)< pDelta32+0.01));
    Data = DataMerged(ind,3*(i-1)+1:3*i);
    ind = find(Data(:,2) > sampleSize);
    Data = Data(ind,:);
    HWEData = Data(:,3)./((Data(:,1)./Data(:,2)).^2.*Data(:,2)/2);
    ind = intersect(find(HWEData ~= 0),find(isnan(HWEData)~=1));
    HWEData = sort(HWEData(ind));
    HWE = DataDelta32(3*i)/((DataDelta32(3*(i-1)+1)/DataDelta32(3*(i-1)+2))^2*DataDelta32(3*(i-1)+2)/2);
    for j = 1:length(HWEData);
        fprintf(fidout,'%d\t%d\n',[i,HWEData(j)]);
    end
end
fclose(fidout);
