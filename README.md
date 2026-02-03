Brain Age Prediction Pipeline
This is a MATLAB analysis pipeline for performing brain age prediction, statistical analysis, and gene expression association studies. The pipeline integrates multimodal neuroimaging data to analyze associations between brain structure changes and age, and explores their spatial correlations with gene expression.

🔧 Main Functional Modules
1. Environment Setup & Dependencies
Configure FSL (FMRIB Software Library) and AFNI (Analysis of Functional NeuroImages) environments

Add custom MATLAB function paths

Support Docker integration

2. Brain Age Prediction Analysis
Use machine learning models to predict brain age

Calculate Predicted Age Difference (PAD)

3. Structure-Age Association Analysis
Statistical correlation analysis between gray matter volume (GMV) and PAD

Voxel-based or region-based analysis methods

Multiple comparison correction

4. Gene Expression Association Analysis
Map gene expression data to subcortical regions (using Tian subregion atlas)

Spatial correlation analysis: MRI features and gene expression patterns

Statistical significance testing

5. Visualization Tools
Generate result visualizations

Display multivariate data using spider plots

Create statistical result plots

