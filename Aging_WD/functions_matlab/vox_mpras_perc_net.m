function vox_pec = vox_mpras_perc_net(yeo_img,subcort_img,cerebellum_img,W_img)
% this function to calcuate the proportion of voxels of significant 
% mri parameters to the voxels of brain networks, including yeo cortex
% network, subcortex and cerebellum
% INPUTs:
%       yeo_img: 3D matrix of yeo network, including 7 cortical network
%       subcort_img, 3D matrix of subcortical mask
%       cerebellum_img: 3D matrix of cerebellum mask
%       W_img: 3D matrix of mri parameters
% OUTPUTs:
%       vox_pec: the proportion of voxels of significant mri parameters to 
%                the voxels of each brain networks

% yeo network
yeo_net = unique(yeo_img);
yeo_net = yeo_net(2:end);
vox_pec = zeros(length(yeo_net)+2,1); % for yeo_net, subcortex, cerebellum
for i = 1:length(yeo_net)
    net_idx = find(yeo_img == yeo_net(i));

    W_net = W_img(net_idx);
    nonz_idx = find(W_net);
    vox_pec(i) = length(nonz_idx) / length(net_idx);
end

% subcortex
subc_idx = find(subcort_img);

W_subcort = W_img(subc_idx);
nonz_idx_s = find(W_subcort);
vox_pec(i+1) = length(nonz_idx_s) / length(subc_idx);

% cerebellum
cereb_idx = find(cerebellum_img);

W_cereb = W_img(cereb_idx);
nonz_idx_c = find(W_cereb);
vox_pec(i+2) = length(nonz_idx_c) / length(cereb_idx);

vox_pec = vox_pec .* 100;
