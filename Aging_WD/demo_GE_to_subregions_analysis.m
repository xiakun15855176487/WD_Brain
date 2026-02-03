%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_GE_to_subregions_analysis.m
%%%
%%% MATLAB script to how map gene expression data to subcortical regions
%%% using inhome matlab code. 
%%%
%%% Original: University of Science and Technology of China, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load gene data
ge_data_path = fullfile('GeneprocessedData/','100DS82scaledRobustSigmoidNSGRNAseqQC1wholeBrain_ROI_NOdistCorrEuclidean.mat');
ge_data = load(ge_data_path);

%% draw sample roi
% subcoritcal samples's ID correspondace to gene expression data
% 35    Left-Thalamus-Proper
% 36    Left-Caudate
% 37    Left-Putamen
% 38    Left-Pallidum
% 39    Left-Hippocampus
% 40    Left-Amygdala
% 76    Right-Thalamus-Proper
% 77    Right-Caudate
% 78    Right-Putamen
% 79    Right-Pallidum
% 80    Right-Hippocampus
% 81    Right-Amygdala
% 82    Right-Accumbens-area
sample_id = [35,36,37,38,39,40,41,76,77,78,79,80,81,82];
sample_coordinates = ge_data.SampleCoordinates;
output_dir = './sub_ROIs';
draw_sample_roi(sample_id,sample_coordinates,output_dir)

%% extract gene expression data of subcortical regions
SampleGeneExpression = ge_data.SampleGeneExpression;
sample_roi_path = './Subcort_ROIs';
SUB_atls_file = fullfile('masks','Tian_Subcortex_S4_3T.nii');
[SParcelExpression,roi_locatedin_atls] = SUB_geneExpression(SampleGeneExpression,sample_roi_path,SUB_atls_file);
% [subParcelExpression,roi_locatedin_atls] = corr_sampled_gene(SampleGeneExpression,sample_roi_path,SUB_atls_file);

save(fullfile("GeneprocessedData","BNSParcelExpression.mat"),"SParcelExpression")
save(fullfile("GeneprocessedData","Sroi_locatedin_atls.mat"),"roi_locatedin_atls")