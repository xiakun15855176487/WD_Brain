%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_SPC_MRIandGE_analysis.m
%%%
%%% MATLAB script to perform ﻿spatial correlation analysis 
%%% for investigting spatial relationships between MRI parameters and gene 
%%%  expression, which implements by using inhome matlab code. 
%%%
%%% Original: University of Science and Technology of China, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare data

% load surfer atlas
surf_atlas_path = fullfile('masks','MMP_MNI152_2mm.nii.gz');
surf_nifti = load_untouch_nii(surf_atlas_path);
surf_img = surf_nifti.img;
surf_subreg = unique(surf_img);
surf_subreg = surf_subreg(2:151); % select left hemesphere

% load subcortical atlas
sub_atlas_path = fullfile('masks','Tian_Subcortex_S4_3T.nii');
sub_atlas_nifti = load_untouch_nii(sub_atlas_path);
sub_atlas_img = sub_atlas_nifti.img;
sub_subreg = unique(sub_atlas_img);
sub_subreg = sub_subreg(2:end); % left subregion

% Project MRI parameters (association between PAD and gray volume) to surf atlas 
group = {'LWD','NWD'};
lc_mriP_group = struct();
sub_mriP_group = struct();
for i = 1:length(group)
    % load MRI parameter
    mri_para_path = fullfile('Results','GV_PAD_Corr',[group{i},'_R.nii']);
    mri_para_nifti = load_untouch_nii(mri_para_path);
    mri_para_img = mri_para_nifti.img;
    mmp_mriP = zeros(length(surf_subreg),1);
    for j = 1:length(surf_subreg)
        sub_idx = find(surf_img == surf_subreg(j));
        sub_mriP = mri_para_img(sub_idx);
        non_idx = find(sub_mriP);
        sub_mriP = sub_mriP(non_idx);
        mmp_mriP(j) = mean(sub_mriP);
    end
    lc_mriP_group.(group{i}) = mmp_mriP;

    subc_mriP = zeros(length(sub_subreg),1);
    for k = 1:length(sub_subreg)
        s_sub_idx = find(sub_atlas_img == sub_subreg(k));
        s_sub_mriP = mri_para_img(s_sub_idx);
        s_non_idx = find(s_sub_mriP);
        s_sub_mriP = s_sub_mriP(s_non_idx);
        subc_mriP(k) = mean(s_sub_mriP);
    end
    sub_mriP_group.(group{i}) = subc_mriP;
end

% load gene expression data
l_cge_path = fullfile('GeneprocessedData','100DS360scaledRobustSigmoidNSGRNAseqQC1Lcortex_ROI_NOdistCorrEuclidean.mat');
l_cge_data = load(l_cge_path); % cortex
l_cgene_Exp = l_cge_data.parcelExpression(:,2:end);
ind_gene = find(~isnan(l_cgene_Exp(:,1)));

sge_path = fullfile('GeneprocessedData','SParcelExpression.mat');
sge_data = load(sge_path); % subcortex
sge_Exp = sge_data.SParcelExpression(:,2:end);
ind_sgene = find(~isnan(sge_Exp(:,1)));

% load gene symptoms
ge_sym_path = fullfile('GeneprocessedData','GeneSymbol.mat');
load(ge_sym_path)

%% spatial correlation analysis
% set output path
outputDIR = './Results/Spatial_Corr';

cort_spat_corr = struct();
subcort_spat_corr = struct();
for i = 1:length(group)
    % cortex
    cmri_paras = lc_mriP_group.(group{i});
    l_cgene_Exp_Na = l_cgene_Exp(ind_gene,:);
    cmri_paras = cmri_paras(ind_gene);
    cort_sc_Res = spatialCorr(l_cgene_Exp_Na,cmri_paras,GeneSymbol);
    cort_spat_corr.(group{i}) = cort_sc_Res;

    % subcortex
    smri_paras = sub_mriP_group.(group{i});
    sind = find(~isnan(smri_paras(:,1)));
    sind_mg = intersect(sind,ind_sgene);

    smri_paras_Na = smri_paras(sind_mg);
    sge_Exp_Na = sge_Exp(sind_mg,:);
    subcort_sc_Res = spatialCorr(sge_Exp_Na,smri_paras_Na,GeneSymbol);
    subcort_spat_corr.(group{i}) = subcort_sc_Res;
end

save(fullfile(outputDIR,'cort_spat_corr_validation.mat'),"cort_spat_corr")
save(fullfile(outputDIR,'subcort_spat_corr_validation.mat'),"subcort_spat_corr")



