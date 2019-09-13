clear all
close all
clc

A = importdata('ukbWESControlSNPsCounts.txt');
A = A.data;
ind = find(A(:,4)> 36953); %41059*0.9
% ind = intersect(find(A(:,5)>= 0.1170-0.0025),find(A(:,5)<= 0.1170+0.0025));
A = A(ind,:);

hold on


%xlim([-0.5 0.5])

gnomad = importdata('Bi_gnomad_individualPopControlSNPs.txt');
delta32Gnomad = importdata('Bi_gnomad_individualPopDelta32.txt');

Array = importdata('ukbArrayControlSNPsCounts.txt');
SNP = Array.textdata;
Array = Array.data;
indDelta32 = find(strcmp(SNP,'rs62625034')==1)
pArray = (Array(:,1)+Array(:,2)/2)./sum(Array,2);
ArrayResult = Array(:,1)./(sum(Array,2).*pArray.*pArray);
arrayHWEDelta32 = ArrayResult(indDelta32);
sum(Array(indDelta32,:))
ind = find(sum(Array,2)>0.9*sum(Array(indDelta32,:)));
ArrayResult = ArrayResult(ind,:);

Cat = [{'FIN'},{'SWE'},{'EST'},{'NWE'},{'ONF'},{'ALL'},{'UKB-WES'},{'UKB-Array'}];
gnomad = [gnomad;[ones(length(A(:,6)),1)*7,A(:,6)+1];[ones(length(ArrayResult),1)*8,ArrayResult]];
delta32Gnomad = [delta32Gnomad;7, 0.9555;8,arrayHWEDelta32];
Frac = [];
CatAll = [];
subplot(2,1,2)
hold on
gnColor = [1. 0 0; 1 0 1; 0 1  1 ; 0 1 0 ;  0 0 1; 0 0 0];
gnLin = [{'--'},{'--'},{'--'},{'--'},{'--'},{'-'}];
for i = 1:6
    ind = find(gnomad(:,1)==i);
    A = gnomad(ind,2);
    h = cdfplot(A)
    set(h,'Linewidth',1.5,'LineStyle',gnLin{i},'Color',gnColor(i,:));
    Frac = [Frac,length(find(A >1))/length(find(A ~=1))];
    CatAll = [CatAll;repmat(Cat(i),length(ind),1)];
end
xlim([0 2])
xlabel('\it{B}_{i}')
ylabel('CDF')
set(gca,'XTick',[0:0.5:2])
title('')
box on
legend(Cat(1:6),'Location','best');
text(-0.2, 1, 'b','FontSize',40);
plot([0,2],[0.025,0.025],'--','HandleVisibility','off','LineWidth',1.5,'Color',[0.5 0.5 0.5]);
plot([0,2],[0.975,0.975],'--','HandleVisibility','off','LineWidth',1.5,'Color',[0.5 0.5 0.5]);
plot([1 1],[0 1],'--','HandleVisibility','off','LineWidth',1.5,'Color',[0.5 0.5 0.5]);
set(gca,'YTick',[0:0.5:1])
set(gca,'FontSize',25)


subplot(2,1,1)

hold on
ukbColor = [1 0 0 ;0 0 1];
for i = 7:8
    ind = find(gnomad(:,1)==i);
    A = gnomad(ind,2);
    h = cdfplot(A);
    ind2 = intersect(find(A > 0.5),find(A < 2));
    B = A(ind2);
    B = sort(B);
    set(h,'Linewidth',1.5,'LineStyle','-','Color',ukbColor(i-6,:));
    Frac = [Frac,length(find(A >1))/length(find(A ~=1))];
    CatAll = [CatAll;repmat(Cat(i),length(ind),1)];
end
text(-0.2, 1, 'a','FontSize',40);
xlim([0 2])
xlabel('\it{B}_{i}')
ylabel('CDF')
set(gca,'XTick',[0:0.5:2])
title('')
box on
legend(Cat(7:8),'Location','best');
Frac
plot([0,2],[0.025,0.025],'--','HandleVisibility','off','LineWidth',1.5,'Color',[0.5 0.5 0.5]);
plot([0,2],[0.975,0.975],'--','HandleVisibility','off','LineWidth',1.5,'Color',[0.5 0.5 0.5]);
plot([1 1],[0 1],'--','HandleVisibility','off','LineWidth',1.5,'Color',[0.5 0.5 0.5]);
set(gca,'YTick',[0:0.5:1])
set(gca,'FontSize',25)
set(gcf,'PaperPosition',[0 0 12 12])
saveas(1,'Figure 1');

%violinplot(gnomad(:,2),CatAll);
%boxplot(gnomad(:,2),CatAll);
%plot([0 8],[1,1],'LineWidth',2)
