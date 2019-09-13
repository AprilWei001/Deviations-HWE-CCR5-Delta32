clear all
close all
clc

SampleSize =[12524  13064  2408  25309  16520 69825  41059  395714];
Range = [0.88 1.17 ;   0.82 1.14;   0.69 1.32 ; 0.90 1.21 ;   0.87 1.18;  0.95 1.20 ; 0.92 1.69; 0.93 1.04];
Med =[1.015   0.995  0.999  1.004  1.017   1.023 1.025   0.992];
DataDelta32 = [228   167    25   302   171   893  537 4349;
    2903        2832         505        5103        3112       14455  8537  83044;
    9393       10065        1878       19904       13237       54477  31985 308321];


MAF = (DataDelta32(1,:)+DataDelta32(2,:)/2)./SampleSize;
DataExpected = [MAF.^2.*SampleSize;2*MAF.*(1-MAF).*SampleSize;(1-MAF).^2.*SampleSize];
ChiSquare = sum((DataDelta32-DataExpected).^2./DataExpected,1);

              
Delta32 =	[1.012  0.871  0.782  0.939  0.947  0.946  0.9555  0.818];
MeanDev = [1.0806  0.9935 0.9999 1.0436 1.0176 1.0759 1.0785 0.9926];
%Pval = [ 0.482   0.057   0.069   0.093   0.135   0.030 0.078 ];
distDelta32 = Delta32;
%---first convert Pval to Zscore, then distDelta32/Zscore(pVal)*1.96 would
%be the CI for the result

Cat = [{'FIN'},{'SWE'},{'EST'},{'NWE'},{'ONF'},{'ALL'},{'UKB-WES'}];

CI = [];
Pval =[];
bootNumber = 10000;
for i = 1:7;
    genoOri = [zeros(DataDelta32(3,i),1);ones(DataDelta32(2,i),1); 2*ones(DataDelta32(1,i),1)];
    N = SampleSize(i);
    HWE = zeros(bootNumber,1);
    parfor j = 1:bootNumber;
        indRand = randsample(N,N,'true');
        bootGeno = genoOri(indRand);
        ind1 = find(bootGeno == 1);
        ind2 = find(bootGeno == 2);
        pBoot = (length(ind2)+length(ind1)/2)/N;
        HWE(j) = length(ind2)/(pBoot^2*N);
    end
    Pval = [Pval;length(find(HWE >= 1))/bootNumber];
    HWE = sort(HWE);
    %Pval = [Pval;length(find(HWE >= 0))/bootNumber];
    CI = [CI; HWE(249), HWE(9751)];
end

Pval
subplot('Position',[0.1 0.1 0.85 0.85])
Color =[0.3 0.3 0.3; 
    1 0.2 1;
    0.3 0.3 1
    
    0 0.5  .5;
    0.5 0 0.5;
    0 1 1; 
    0.3 .9 0.3;
    .8 0.2 0.2;];
plot([0 10],[1 1],'--','LineWidth',2,'Color',[0.5 0.5 0.5]);
hold on
for i = 1:7;
    plot([i i],CI(i,:),'LineWidth',4,'Color',Color(i,:));
    plot([i],[distDelta32(i)],'.','Color',Color(i,:),'LineWidth',3,'MarkerSize',30);
    plot([i],Med(i),'p','Color','r','LineWidth',2,'MarkerSize',20);
    %plot([i-0.05,i],[distDelta32(i)-0.01,distDelta32(i)],'-','Color','r','LineWidth',1); 
    %plot([i-0.05,i],[distDelta32(i)+0.01,distDelta32(i)],'-','Color','r','LineWidth',1);
end

set(gca,'XTick',[1:7])
xlim([0.5 7.5])
ylabel('Deviations from HWE')%,'Interpreter','latex');
box on
set(gca,'XTickLabel',Cat)
set(gca,'YTick', [0.4:0.2:1.2])
set(gca,'YTickLabel',[{'0.4'},{'0.6'},{'0.8'},{'1.0'},{'1.2'}]);
%ylim([0.6 1.8])
set(gca,'FontSize',25)
set(gcf,'PaperPosition',[0 0 18 6])
saveas(1,'Figure 2.png');
