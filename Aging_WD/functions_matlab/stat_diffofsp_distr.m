function [qnet,rnet] = stat_diffofsp_distr(yeo_img,subcort_img,cerebellum_img,W1_img,W2_img)
% this function to calculate the spatial distribution differences of  
% mri parameters between two groups within each brain network, including 
% yeo cortex network, subcortex and cerebellum
% INPUTs:
%       yeo_img: 3D matrix of yeo network, including 7 cortical network
%       subcort_img, 3D matrix of subcortical mask
%       cerebellum_img: 3D matrix of cerebellum mask
%       W1_img: 3D matrix of mri parameters for group1
%       W2_img: 3D matrix of mri parameters for group2
% OUTPUTs:
%       qnet: q value of each network---FDR-corrected
%       rnet: mean value of each network

% yeo network
yeo_net = unique(yeo_img);
yeo_net = yeo_net(2:end);
pnet = zeros(length(yeo_net)+2,1); % for yeo_net, subcortex, cerebellum
rnet = zeros(length(yeo_net)+2,2); 
for i = 1:length(yeo_net)
    net_idx = find(yeo_img == yeo_net(i));

    W_net1 = W1_img(net_idx);
    nonz_idx1 = find(W_net1);
    W_net1 = W_net1(nonz_idx1);

    W_net2 = W2_img(net_idx);
    nonz_idx2 = find(W_net2);
    W_net2 = W_net2(nonz_idx2);

    [h,p,ci,stats] = ttest2(W_net1,W_net2);
    pnet(i) = p;
    rnet(i,1) = mean(W_net1);
    rnet(i,2) = mean(W_net2);
end

% subcortex
subc_idx = find(subcort_img);

W_subcort1 = W1_img(subc_idx);
nonz_idx_s1 = find(W_subcort1);
W_subcort1 = W_subcort1(nonz_idx_s1);

W_subcort2 = W2_img(subc_idx);
nonz_idx_s2 = find(W_subcort2);
W_subcort2 = W_subcort2(nonz_idx_s2);

[h,p,ci,stats] = ttest2(W_subcort1,W_subcort2);
pnet(i+1) = p;
rnet(i+1,1) = mean(W_subcort1);
rnet(i+1,2) = mean(W_subcort2);

% cerebellum
cereb_idx = find(cerebellum_img);

W_cereb1 = W1_img(cereb_idx);
nonz_idx_c1 = find(W_cereb1);
W_cereb1 = W_cereb1(nonz_idx_c1);

W_cereb2 = W2_img(cereb_idx);
nonz_idx_c2 = find(W_cereb2);
W_cereb2 = W_cereb2(nonz_idx_c2);

[h,p,ci,stats] = ttest2(W_cereb1,W_cereb2);
pnet(i+2) = p;
rnet(i+2,1) = mean(W_cereb1);
rnet(i+2,2) = mean(W_cereb2);


qnet = mafdr(pnet,'BHFDR',true);
