%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_statistical_analysis.m
%%%
%%% MATLAB script to perform statsitical analysis using inhome matlab code. 
%%% 
%%% Original: University of Science and Technology of China, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data
% brain age prediction results
BA = load('BAResults.mat');

% load brain mask
gm_mask_file = fullfile('Data','Train_VBM','GM_mask_s.nii.gz');
gm_mask_nift = load_untouch_nii(gm_mask_file);
gm_mask_img = gm_mask_nift.img;

% yeo 7 cortical network
yeo_file = fullfile('masks','Yeo_7Network_MNI152_2mm.nii.gz');
yeo_nifti = load_untouch_nii(yeo_file);
yeo_img = yeo_nifti.img;

% subcortex
subcort_file = fullfile('masks','subcort_atlas_2mm.nii.gz');
subcort_nifti = load_untouch_nii(subcort_file);
subcort_img = subcort_nifti.img;

% cerebellum
cerebellum_file = fullfile('masks','CRB_mask_2mm.nii.gz');
cerebellum_nifti = load_untouch_nii(cerebellum_file);
cerebellum_img = cerebellum_nifti.img;

%% Step1: correlated clincial information with PAD_corrected

% load clnical information for WD
group = {'LWD','NWD'}; % LWD, NWD
corr_pad_clinc = zeros(length(group),2);
for i = 1:length(group)
    % load clinical information
    T = readtable(fullfile('ClinicInfo','WD_Sort_Clin_Info.xlsx'),'Sheet',group{i});
    clinic = T.urinary; % UWDRS-N Diseaseduration Ceruloplasmin urinary
    pad = BA.Results.(group{i}).PAD_corrected;
    [r,p] = corr(clinic,pad);

    corr_pad_clinc(i,1) = r;
    corr_pad_clinc(i,2) = p;
end

%% Step2: calculate group differences of gray volume between WD groups and HC
output_path = './Results/GV_Diff';
w_groups = {'CWD','HWD','LWD','NWD'};
h_group = 'HC';

gm_indx = find(gm_mask_img);
% load gray matter volume of HC
H_GM_MTX = load([h_group,'_data.mat']);
H_GM_MTX = H_GM_MTX.([h_group,'_data']).gv;

for i = 1:length(w_groups)
    % load WD gray matter volume
    W_GM_MTX = load([w_groups{i},'_data.mat']);
    W_GM_MTX = W_GM_MTX.([w_groups{i},'_data']).gv;

    sigd_gvt_3d = zeros(size(gm_mask_img));

    % two sample t test
    [h,p,ci,stats] = ttest2(W_GM_MTX,H_GM_MTX);
    fdr = mafdr(p,'BHFDR',true);
    sig_idx = find(fdr <= 0.05);
    t_value = zeros(length(gm_indx),1);
    t_value(sig_idx) = stats.tstat(sig_idx);
    sigd_gvt_3d(gm_indx) = t_value;

    % map back to brain
    nifti_data = gm_mask_nift;
    nifti_data.hdr.dime.datatype=16;
    nifti_data.hdr.dime.bitpix=32;
    
    nifti_data.img = sigd_gvt_3d;
    save_untouch_nii(nifti_data,fullfile(output_path,[w_groups{i},'_gv_tvalue.nii']));

end

%% Step3: correlated PAD with gray volume of WD
output_path = './Results/GV_PAD_Corr';

group = {'CWD','HWD','LWD','NWD'};
gm_indx = find(gm_mask_img);
for i = 1:length(group)
    % pad corrected
    PAD_MTX = BA.Results.(group{i}).PAD_corrected;
    
    % load gray matter volume
    GM_MTX = load([group{i},'_data.mat']);
    GM_MTX = GM_MTX.([group{i},'_data']).gv;
    
    % Pearson's correlation
    [R_gv,P_gv] = corr(GM_MTX,PAD_MTX);

    R_3d = zeros(size(gm_mask_img));
    R_3d(gm_indx) = R_gv;

    % FDR correct
    fdr_gv = mafdr(P_gv,'BHFDR',true);
    
    % extract statistcal significance of R and FDR
    sig_indx = find(fdr_gv <= 0.05);

    sigR_MTX = zeros(length(R_gv),1);
    sigfdr_MTX = zeros(length(R_gv),1);

    sigR_MTX(sig_indx) = R_gv(sig_indx);
    sigfdr_MTX(sig_indx) = 1-fdr_gv(sig_indx);

    % map back to brain
    sigR_3d = zeros(size(gm_mask_img));
    sigfdr_3d = zeros(size(gm_mask_img));

    sigR_3d(gm_indx) = sigR_MTX;
    sigfdr_3d(gm_indx) = sigfdr_MTX;

    nifti_data = gm_mask_nift;
    nifti_data.hdr.dime.datatype=16;
    nifti_data.hdr.dime.bitpix=32;
    
    nifti_data.img = R_3d;
    save_untouch_nii(nifti_data,fullfile(output_path,[group{i},'_R.nii']));

    nifti_data.img = sigR_3d;
    save_untouch_nii(nifti_data,fullfile(output_path,[group{i},'_sigR.nii']));

    nifti_data.img = sigfdr_3d;
    save_untouch_nii(nifti_data,fullfile(output_path,[group{i},'_sigfdr.nii']));
end

%% Step4: calculate spatial differences of coefficients of PAD-GMV in each network

% load data
NWD_file = fullfile('Results','GV_PAD_Corr','NWD_sigR.nii');
NWD_nifti = load_untouch_nii(NWD_file);
NWD_img = NWD_nifti.img;

LWD_file = fullfile('Results','GV_PAD_Corr','LWD_sigR.nii');
LWD_nifti = load_untouch_nii(LWD_file);
LWD_img = LWD_nifti.img;
[qnet,rnet] = stat_diffofsp_distr(yeo_img,subcort_img,cerebellum_img,LWD_img,NWD_img);

%% Step5: gene overlapped analysis

% set output path
outputDIR = './Results/Spatial_Corr';

% load data
load(fullfile('Results','Spatial_Corr','cort_spat_corr.mat'))
load(fullfile('Results','Spatial_Corr','subcort_spat_corr.mat'))
load(fullfile('GeneprocessedData','cellTypeSpecific_GE.mat'))

% cortex
G1_ranked_GE = cort_spat_corr.LWD.gene;
G2_ranked_GE = cort_spat_corr.NWD.gene;
% subcortex
S_G1_ranked_GE = subcort_spat_corr.LWD.gene;
S_G2_ranked_GE = subcort_spat_corr.NWD.gene;
% cell type-specific gene
spec_GE = specific_gene_symbol;
g_name = uniq_g_name;
% permutation times
times = 10000;

% Calculate the top 5% gene expression overlap betwenn LWD and NWD
% Cortex
sig_gene_ov_LWD_NWD = sigGE_ov_bewGroup(G1_ranked_GE,G2_ranked_GE,times);
save(fullfile(outputDIR,'cort_sig_gene_ov_LWD_NWD.mat'),"sig_gene_ov_LWD_NWD")

% Subcortex
sig_gene_ov_LWD_NWD = sigGE_ov_bewGroup(S_G1_ranked_GE,S_G2_ranked_GE,times);
save(fullfile(outputDIR,'subcort_sig_gene_ov_LWD_NWD.mat'),"sig_gene_ov_LWD_NWD")

% Calculate the differences of the overlapped cell type specific gene with
% significant gene between LWD and NWD (just for cortex)
cort_ovlp_SP_SIG_GE_results = specGE_ov_rankGE(G1_ranked_GE,G2_ranked_GE,spec_GE,g_name,times);
save(fullfile(outputDIR,'cort_ovlp_SP_SIG_GE_results.mat'),"cort_ovlp_SP_SIG_GE_results")

%% Step6: enrichment analysis for disease risk gene
% The candidate gene include Wilson's disease (WD), Mitochondrial disorders (MD),
% Parkinson's disease (PD), Alzheimer's Disease (AD), Huntinton's disease (HD),
% and dystonia.

% set output path
outputDIR = './Results/Spatial_Corr';

% load disease gene symbols
T = readtable(fullfile('GeneprocessedData','Dis_Gene.xlsx'));
headers = T.Properties.VariableNames;
headers = headers(2:end);

% load rank gene and their weight of spatial correlation
load(fullfile('Results','Spatial_Corr','cort_spat_corr.mat')) % cortex
load(fullfile('Results','Spatial_Corr','subcort_spat_corr.mat')) % subcortex

% enrichment ratio analysis
groups = {'LWD','NWD'};
ER_cort = struct();
ER_subcort = struct();
for i = 1:length(groups)
    cort_SPC = cort_spat_corr.(groups{i});
    subcort_SPC = subcort_spat_corr.(groups{i});

    for j = 1:length(headers)
        Cgene = T.(headers{j});
        % cortex
        [cort_p(j),cort_ER(j)] = calculate_candidate_genes_ER(cort_SPC,Cgene);
        
        % subcortex
        [subcort_p(j),subcort_ER(j)] = calculate_candidate_genes_ER(subcort_SPC,Cgene);
    end
    cort_fdr = mafdr(cort_p,'BHFDR',true);
    ER_cort.(groups{i}).fdr = cort_fdr;
    ER_cort.(groups{i}).ER = cort_ER;
    
    subcort_fdr = mafdr(subcort_p,'BHFDR',true);
    ER_subcort.(groups{i}).fdr = subcort_fdr;
    ER_subcort.(groups{i}).ER = subcort_ER;
end
save(fullfile(outputDIR,'ER_cort.mat'),"ER_cort")
save(fullfile(outputDIR,'ER_subcort.mat'),"ER_subcort")

%% Step7: Associated PAD with GV in LWD and NWD respectively based on SVM weight

% load SVM weitht
svm_w_path = fullfile('Results','SVM','SVM_weights.nii');
svm_w_nii = load_untouch_nii(svm_w_path);
svm_w_img = svm_w_nii.img;

% load BAresults
load BAResults.mat

% load gv data
load LWD_data.mat
load NWD_data.mat

% load GM mask
gm_tem_path = fullfile('Data','Train_VBM','GM_mask_s.nii.gz');
gm_tem_nii = load_untouch_nii(gm_tem_path);
gm_tem_img = gm_tem_nii.img;
gm_ind = find(gm_tem_img);

p_ind = find(svm_w_img > 0);  p_in_gm_ind = whereind(p_ind,gm_ind);  % NWD 
n_ind = find(svm_w_img < 0);  n_in_gm_ind = whereind(n_ind,gm_ind);  % LWD

gm_p = NWD_data.gv(:,p_in_gm_ind); mean_gm_p = mean(gm_p,2);
gm_n = LWD_data.gv(:,n_in_gm_ind); mean_gm_n = mean(gm_n,2);

mean_gm_pn = mean(NWD_data.gv(:,n_in_gm_ind),2);
mean_gm_np = mean(LWD_data.gv(:,p_in_gm_ind),2);

pad_p = Results.NWD.PAD_corrected;
pad_n = Results.LWD.PAD_corrected;

[rho_p,pval_p] = corr(mean_gm_p,pad_p); [rho_pn,pval_pn] = corr(mean_gm_pn,pad_p);
[rho_n,pval_n] = corr(mean_gm_n,pad_n); [rho_np,pval_np] = corr(mean_gm_np,pad_n);

% load clinical information
group = {'LWD','NWD'}; % LWD, NWD
T = readtable(fullfile('ClinicInfo','WD_Sort_Clin_Info.xlsx'),'Sheet',group{1});
clinic = T.UWDRS_N;
[rho_n_clin,pval_n_clin] = corr(clinic,mean_gm_n);
[rho_np_clin,pval_np_clin] = corr(clinic,mean_gm_np);

T = readtable(fullfile('ClinicInfo','WD_Sort_Clin_Info.xlsx'),'Sheet',group{2});
clinic = T.UWDRS_N;
[rho_p_clin,pval_p_clin] = corr(clinic,mean_gm_p);
[rho_pn_clin,pval_pn_clin] = corr(clinic,mean_gm_pn);