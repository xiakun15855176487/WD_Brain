function sig_gene_ov = sigGE_ov_bewGroup(G1_ranked_GE,G2_ranked_GE,times)
% this function is to calculate significant gene overlapped between two
% groups.The ranked gene from spatial correlation analysis
% INPUTs:
%        G1_ranked_GE: ranked gene of group-1 from spatial correlation
%        analysis（group 1 is LWD）
%        G2_ranked_GE: ranked gene of group-2 from spatial correlation
%        analysis (group 2 is NWD)
%        times: permutation times
% OUTPUTs:
%         sig_gene_ov: save all results of gene overlapped results

G1_sig_GE_Pos = G1_ranked_GE(1:250);
G2_sig_GE_Pos = G2_ranked_GE(1:250);

G1_sig_GE_Neg = G1_ranked_GE(10027-249:end);
G2_sig_GE_Neg = G2_ranked_GE(10027-249:end);

genesymbolOv_Pos = genesymbolOverlap(G1_sig_GE_Pos,G2_sig_GE_Pos);
genesymbolOv_Neg = genesymbolOverlap(G1_sig_GE_Neg,G2_sig_GE_Neg);
num_ov_Pos = length(genesymbolOv_Pos);
num_ov_Neg = length(genesymbolOv_Neg);

% Permutaion
num_ov_Pos_Perm = zeros(1,times);
num_ov_Neg_Perm = zeros(1,times);
parfor i = 1:times
    fprintf('The %d times of permutation\n',i);
    ge_order1 = randperm(length(G1_ranked_GE));
    ge_order2 = randperm(length(G2_ranked_GE));

    G1_ranked_GE_Perm = G1_ranked_GE(ge_order1);
    G2_ranked_GE_Perm = G2_ranked_GE(ge_order2);

    G1_sig_GE_Pos_Perm = G1_ranked_GE_Perm(1:250);
    G2_sig_GE_Pos_Perm = G2_ranked_GE_Perm(1:250);

    G1_sig_GE_Neg_Perm = G1_ranked_GE_Perm(10027-249:end);
    G2_sig_GE_Neg_Perm = G2_ranked_GE_Perm(10027-249:end);

    genesymbolOv_Neg_Perm = genesymbolOverlap(G1_sig_GE_Neg_Perm,G2_sig_GE_Neg_Perm);
    genesymbolOv_Pos_Perm = genesymbolOverlap(G1_sig_GE_Pos_Perm,G2_sig_GE_Pos_Perm);
    
    num_ov_Pos_Perm(i) = length(genesymbolOv_Pos_Perm);
    num_ov_Neg_Perm(i) = length(genesymbolOv_Neg_Perm);

end
P_Pos = length(find(num_ov_Pos_Perm > num_ov_Pos)) / times;
P_Neg = length(find(num_ov_Neg_Perm > num_ov_Neg)) / times;

sig_gene_ov = struct();
sig_gene_ov.P_Pos = P_Pos;
sig_gene_ov.P_Neg = P_Neg;
sig_gene_ov.num_ov_Pos = num_ov_Pos;
sig_gene_ov.num_ov_Neg = num_ov_Neg;