clear all
close all
clc

fid = fopen('gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf');
fout= fopen('gnomad_WGS_exome_fin_est_nwe_onf.txt','w');

%change the imput output file below to process the WES or NGS data from gnomAD.
%fid = fopen('gnomad.exomes.r2.1.1.sites.vcf');
%fout= fopen('gnomad_WES_fin_swe_est_nwe_onf.txt','w');

tline = fgetl(fid);

while strcmp(tline(1),'#') == 1
    tline = fgetl(fid);
end
matchFields = [{'AC_fin'},{'AN_fin'},{'AF_fin'},{'nhomalt_fin'},...
    {'AC_nfe_swe'},{'AN_nfe_swe'},{'AF_nfe_swe'},{'nhomalt_nfe_swe'},...
    {'AC_nfe_est'},{'AN_nfe_est'},{'AF_nfe_est'},{'nhomalt_nfe_est'},...
    {'AC_nfe_nwe'},{'AN_nfe_nwe'},{'AF_nfe_nwe'},{'nhomalt_nfe_nwe'},...
    {'AC_nfe_onf'},{'AN_nfe_onf'},{'AF_nfe_onf'},{'nhomalt_nfe_onf'},...
    ];



while ischar(tline)
    pos1 = strfind(tline,';AC_fin=');
    if  (tline(pos1+8) ~= '0') 
       % pos2 = strfind(tline,';AC_nfe_swe=');%swe does not exist in WGS
    %if (tline(pos2+8) ~= '0') 
        pos3 = strfind(tline,';AC_nfe_est=');
    if    (tline(pos3+8) ~= '0') 
        pos4 = strfind(tline,';AC_nfe_nwe=');
    if    (tline(pos4+8) ~= '0') 
        pos5 = strfind(tline,';AC_nfe_onf=');
    if    (tline(pos5+8) ~= '0')
        SNPinfo = strsplit(tline(1:200),';');
        SNPheader = SNPinfo(1);
        Fin_info = strsplit(tline(pos1+1:pos1+100),';');
        Fin_info = Fin_info(1:4);
       % SWE_info = strsplit(tline(pos2+1:pos2+100),';');
       % SWE_info = SWE_info(1:4);
        EST_info = strsplit(tline(pos3+1:pos3+100),';');
        EST_info = EST_info(1:4);
        NWE_info = strsplit(tline(pos4+1:pos4+100),';');
        NWE_info = NWE_info(1:4);
        ONF_info = strsplit(tline(pos5+1:pos5+100),';');
        ONF_info = ONF_info(1:4);
        SNPfin = strjoin(cellfun(@(s) strcat(s{2},'\t'),cellfun(@(s) strsplit(s,'='),...
             Fin_info(1:4),'UniformOutput',false),'UniformOutput',false),'');
        %SNPswe = strjoin(cellfun(@(s) strcat(s{2},'\t'),cellfun(@(s) strsplit(s,'='),...
        %     SWE_info(1:4),'UniformOutput',false),'UniformOutput',false),'');
        SNPest = strjoin(cellfun(@(s) strcat(s{2},'\t'),cellfun(@(s) strsplit(s,'='),...
             EST_info(1:4),'UniformOutput',false),'UniformOutput',false),'');
        SNPnwe = strjoin(cellfun(@(s) strcat(s{2},'\t'),cellfun(@(s) strsplit(s,'='),...
             NWE_info(1:4),'UniformOutput',false),'UniformOutput',false),'');
        SNPonf = strjoin(cellfun(@(s) strcat(s{2},'\t'),cellfun(@(s) strsplit(s,'='),...
             ONF_info(1:4),'UniformOutput',false),'UniformOutput',false),'');
        fprintf(fout,SNPheader{1});
        fprintf(fout,'\t');
        fprintf(fout,SNPfin(1:end));
       % fprintf(fout,SNPswe(1:end));
        fprintf(fout,SNPest(1:end));
        fprintf(fout,SNPnwe(1:end));
        fprintf(fout,SNPonf(1:end-2));
        fprintf(fout,'\n');
    %end
    end
    end
    end
    end
    tline = fgetl(fid);
end

fclose(fout);
fclose(fid);
%fout_ageEff = fopen('gnomad_delta32.txt','w');
%     if contains(tline,'TACAGTCAGTATCAATTCTGGAAGAATTTCCAG')
%         fprintf(fout_ageEff,tline);
%         fprintf(fout_ageEff,'\n');
%     end
%fclose(fout_ageEff);
% 
% matchFields = [{'AC_fin'},{'AN_fin'},{'nhomalt_fin'},...
%     {'AC_nfe'},{'AN_nfe'},{'nhomalt_nfe'},...
%     {'age_hist_het_bin_freq'},{'age_hist_het_n_smaller'},{'age_hist_het_n_larger'},...
%     {'age_hist_hom_bin_freq'},{'age_hist_hom_n_smaller'},{'age_hist_hom_n_larger'}];
% 
% Fields = strsplit(tline,';');
% 
% for i = 2:length(Fields);
%     field = strsplit(char(Fields(i)),'=');
%     Fields{i} = field{1};
% end
% lenFields = length(matchFields);
% indFields = [];
% for i = 1:lenFields
%     ind = find(strcmp(Fields,matchFields{i})==1);
%     if length(ind)==1
%     matchFields{i}
%     end
%     indFields = [indFields;ind];
% end
