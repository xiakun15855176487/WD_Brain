function ovlp_GE_results = specGE_ov_rankGE(G1_ranked_GE,G2_ranked_GE,spec_GE,g_name,times)
% this function is to calculate gene overlapped between cell-type specific
% gene and ranked gene from spatial correlation analysis, and calculate the
% the differnces of the number of overlapped gene between two groups
% INPUTs:
%        G1_ranked_GE: ranked gene of group-1 from spatial correlation
%        analysis（group 1 is LWD）
%        G2_ranked_GE: ranked gene of group-2 from spatial correlation
%        analysis (group 2 is NWD)
%        spec_GE: cell-type specific gene
%        times: permutation times
% OUTPUTs:
%         ovlp_GE_results: save all results of gene overlapped results

% group 1
sig_GE_Pos1 = G1_ranked_GE(1:250);
sig_GE_Neg1 = G1_ranked_GE(10027-249:end);
spec_ovlp_GE_Pos1 = cell(length(spec_GE),1);
spec_ovlp_GE_Neg1 = cell(length(spec_GE),1);

% group 2
sig_GE_Pos2 = G2_ranked_GE(1:250);
sig_GE_Neg2 = G2_ranked_GE(10027-249:end);
spec_ovlp_GE_Pos2 = cell(length(spec_GE),1);
spec_ovlp_GE_Neg2 = cell(length(spec_GE),1);

P_pos = zeros(1,length(spec_GE));
P_neg = zeros(1,length(spec_GE));
for i = 1:length(spec_GE)
    fprintf('The %dth cell typer-specifc gene (%s)\n',i,g_name(i))
    spec_genesym = spec_GE{i};
    spec_genesym = unique(spec_genesym); % remove the repeated gene symbols
    
    % group 1
    spec_ovlp_GE_Pos1{i} = genesymbolOverlap(spec_genesym,sig_GE_Pos1); % Positive ranked gene
    spec_ovlp_GE_Neg1{i} = genesymbolOverlap(spec_genesym,sig_GE_Neg1); % Negative ranked gene 

    % group 2
    spec_ovlp_GE_Pos2{i} = genesymbolOverlap(spec_genesym,sig_GE_Pos2); % Positive ranked gene
    spec_ovlp_GE_Neg2{i} = genesymbolOverlap(spec_genesym,sig_GE_Neg2); % Negative ranked gene 
    
    % differnce of the number of overlapped gene between two groups
    diff_num_Pos = abs(length(spec_ovlp_GE_Pos1{i})-length(spec_ovlp_GE_Pos2{i}));
    diff_num_Neg = abs(length(spec_ovlp_GE_Neg1{i})-length(spec_ovlp_GE_Neg2{i}));

    % permutation
    dif_num_Pos_per = zeros(1,times);
    dif_num_Neg_per = zeros(1,times);
    parfor j = 1:times
        fprintf('The %d times of permutation\n',j);
        ge_order1 = randperm(length(G1_ranked_GE));
        ge_order2 = randperm(length(G1_ranked_GE));

        % group 1
        G1_ranked_GE_Perm = G1_ranked_GE(ge_order1);
        sig_GE_Pos1_Perm = G1_ranked_GE_Perm(1:250);
        sig_GE_Neg1_Perm = G1_ranked_GE_Perm(10027-249:end);
        ovlp_GE_Pos1_Perm = genesymbolOverlap(spec_genesym,sig_GE_Pos1_Perm);
        ovlp_GE_Neg1_Perm = genesymbolOverlap(spec_genesym,sig_GE_Neg1_Perm);

        % group2
        G2_ranked_GE_Perm = G1_ranked_GE(ge_order2);
        sig_GE_Pos2_Perm = G2_ranked_GE_Perm(1:250);
        sig_GE_Neg2_Perm = G2_ranked_GE_Perm(10027-249:end);
        ovlp_GE_Pos2_Perm = genesymbolOverlap(spec_genesym,sig_GE_Pos2_Perm);
        ovlp_GE_Neg2_Perm = genesymbolOverlap(spec_genesym,sig_GE_Neg2_Perm);

        dif_num_Pos_per(j) = abs(length(ovlp_GE_Pos1_Perm)-length(ovlp_GE_Pos2_Perm));
        dif_num_Neg_per(j) = abs(length(ovlp_GE_Neg1_Perm)-length(ovlp_GE_Neg2_Perm));
    end
    P_pos(i) = length(find(dif_num_Pos_per > diff_num_Pos)) / times;
    P_neg(i) = length(find(dif_num_Neg_per > diff_num_Neg)) / times;
end

% save results
ovlp_GE_results = struct();
ovlp_GE_results.cellTypeName = g_name;
ovlp_GE_results.OVLP_Pos1 = spec_ovlp_GE_Pos1;
ovlp_GE_results.OVLP_Neg1 = spec_ovlp_GE_Neg1;
ovlp_GE_results.OVLP_Pos2 = spec_ovlp_GE_Pos2;
ovlp_GE_results.OVLP_Neg2 = spec_ovlp_GE_Neg2;
ovlp_GE_results.P_pos = P_pos;
ovlp_GE_results.P_neg = P_neg;