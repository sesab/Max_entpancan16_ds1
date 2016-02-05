load('pancan16_ds1.mat'); % raw observation data

data = mat_1; % mat_1: real observation. mat_2, mat_3, mat_4: random matrices (frequency preserving permutation)

% Set parameters and fix data
data = data + 1; % add 1 as simbol 0 is not supported (0,1 -> 1,2)
q=max(max(data));
[M,N] = size(data);
pseudocount_weight = 0.000; % relative weight of pseudo count   
theta = 0; % threshold for sequence id in reweighting
p=N; % total number of eigenvalues to accept
eig_min=1; % accept all eigenvalues <= eig_min
eig_max=1; % accept all eigenvalues > eig_min

% DCA
[F_apc, accepted, invC, C, D, Gamma, Lambda] = HopfieldPottsDCA(data, N, M, q, pseudocount_weight, theta, p, eig_min, eig_max);
