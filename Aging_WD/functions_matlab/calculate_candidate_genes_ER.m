function [p,ER] = calculate_candidate_genes_ER(SPC,Cgene)
% This function performs permutation test for significance of up/down
% weighting of candidate gene lists in a SPC
% INPUTs:
        % SPC:results from the spatial correlation, containing correlation
        % coefficient for each ranked gene
        % Cgene:the candidate disease risk genes
% ABSOLUTE: ABSOLUTE=true for candiate gene set and ABSOLUTE=false for oligo gene set

disp('  Running  enrichment ratio analysis for candidate gene')

SPCweight=SPC.rho;
SPCweight = zscore(SPCweight);
SPCgenes=SPC.gene;

CANDind = [];
for i = 1:length(Cgene) 
    for j = 1:length(SPCgenes)
        if strcmp(Cgene(i),SPCgenes(j))
            CANDind = [CANDind,j];
            break
        end
    end
end

CAND_SPCwe = SPCweight(CANDind);

mean_CAND_SPCw = mean(CAND_SPCwe);

perm_SPCw=[];
for r=1:10000
    permo=randperm(length(SPCgenes));
    tem_SPCw=SPCweight(permo(1:length(Cgene)));
    
    mean_tem_SPCw=mean(tem_SPCw);
    
    perm_SPCw=[perm_SPCw;mean_tem_SPCw];
end
perm_SPCw1 = abs(perm_SPCw);
n = length(find(perm_SPCw1 > abs(mean_CAND_SPCw)));
p = n/length(perm_SPCw1);
% p = p/2;
ER = (mean_CAND_SPCw - mean(perm_SPCw)) / std(perm_SPCw);