%% Some explain for the matlab scripts
% In the matlab scripts, varibles of LWD and NWD are corresponding to 
% PWD and DWD in the article respectively.

%% add inhome matlab function path
addpath ./functions_matlab
addpath ./functions_matlab/spider_plot/

%% Set FSL and afni enviroment
setenv('FSLDIR', '/Users/devlinhu/fsl'); 
setenv('PATH', [getenv('PATH') ':' getenv('FSLDIR') '/bin']); % fsl
setenv('PATH', [getenv('PATH') ':/Users/devlinhu/abin']); % afni
setenv('PATH', [getenv('PATH') ':/usr/local/bin']); % docker

%% Brain age prediction
demo_BA_prediction_analysis

%% Correlated PAD with gray matter volume (GMV) of WD----PAD-GMV
demo_statistical_analysis % step3 in this functional script

%% Spatial correlation analysis between gene expression and PAD-GMV
% map the gene data to subcoritcal atlas,here using Tian_Sub atlas
% demo_GE_to_subregions_analysis

demo_SPC_MRIandGE_analysis

%% Statictical analysis
demo_statistical_analysis

%% Visualization of results
demo_visualization
