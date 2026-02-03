function sc_Res = spatialCorr(ge_EX,mri_paras,GeneSymbol)
% this function is to calculate spatial correlation between gene expression
% and mri parameters, and then rank the gene symbol accorrding to
% correlation coefficient of each gene
% INPUTs:
%       ge_EX: gene expression data
%       mri_paras: mri parameters, such as gray volume
%       GeneSybol: gene symbol
% OUTPUTs:
%       sc_Res: save all the spatial correlation results

sc_Res = struct();
% Spearman correlation analysis
[rho,pval] = corr(ge_EX,mri_paras,"type","Spearman");
rho_PEM = [];
parfor j = 1:10000
    order=randperm(size(mri_paras,1));
    mri_paras_perm = mri_paras(order);
    [rho_perm,pval_perm] = corr(ge_EX,mri_paras_perm,"type","Spearman");
    rho_PEM = [rho_PEM,rho_perm];
end
count = sum(abs(rho_PEM) > abs(rho),2);
pval_PEM = count ./ 10000;
[rho, rho_idx]= sort(rho,'descend');
pval_PEM = pval_PEM(rho_idx);
fdr_PEM = mafdr(pval_PEM,'BHFDR',true);
sort_geneSym = GeneSymbol(rho_idx);
    
sc_Res.rho = rho;
sc_Res.fdr = fdr_PEM;
sc_Res.gene = sort_geneSym;