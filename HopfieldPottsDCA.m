function [F_apc, accepted, invC, C, D, Gamma, Lambda] = HopfieldPottsDCA(data, N, M, q, pseudocount_weight, theta, p, min_eig, max_eig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Hopfield-Potts approach to Direct-Coupling Aanalysis 
%
% Copyright for this implementation: 
%             2013 - Martin Weigt
%                    martin.weigt@upmc.fr
% 
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
% Any publication resulting from applications of DCA should cite:
%
%     S Cocco, R Monasson and M Weigt (2013)
%     From principal component analysis to direct-coupling analysis
%     of coevolution in proteins: Low-eigenvalue modes are needed
%     for structure prediction, PLoS Computational Biology XX, XXXX
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Pij_true,Pi_true, Meff]=Compute_True_Frequencies(data,M,N,q, theta); % read binary alteration matrix and build pairwise matrix. Applies theta correction
    fprintf('### N = %d M = %d Meff = %.2f q = %d theta = %.2f pcw = %.2f p = %d\n', N,M,Meff,q, theta, pseudocount_weight, p); % print diagnostic stats
    [Pij,Pi] = with_pc(Pij_true,Pi_true,pseudocount_weight,N,q); % apply pseudocount correction
    %Pij = Pij_true; % unmask to remove pseudocount correction 
    %Pi = Pi_true; % unmask to remove pseudocount correction

    [C,C_self] = Compute_C(Pij,Pi,N,q); % calculate covariance matrix C

    D = sqrtm(C_self);

    Gamma = (D\C)/D;     % Pearson correlation matrix Gamma

    [V,Lambda] = eig(Gamma);  

    Vtilde = D\V;    % Patterns up to prefactor sqrt(abs(1-1/lambda))

    Lambda1 = eye( N*(q-1));

    % log-likelihood contribution of patterns 

    ll = diag(Lambda) - ones( N*(q-1),1 ) - log( diag( Lambda ) );

    [~,b] = sort( ll, 'descend');

    % consider only p patterns of highest ll
    accepted = 0;
    for i=1:p
        if (Lambda( b(i), b(i) ) <= min_eig) || (Lambda( b(i), b(i) ) >= max_eig)
            Lambda1( b(i), b(i) ) = Lambda( b(i), b(i) );
            accepted = accepted + 1;
        end
    end
    Lambda2 = ( Lambda1 - eye(N*(q-1)) )/Lambda1;
    invC = -Vtilde * Lambda2 * Vtilde';
    F_apc = old_calc_norm_apc( invC, N, q); % Calculate final scores - sampling-corrected Frobenius norm of the couplings
end


function [Pij_true,Pi_true,Meff] = Compute_True_Frequencies(align,M,N,q,theta)
% computes reweighted frequency counts
    W = ones(1,M);
    if( theta > 0.0 )
        W = (1./(1+sum(squareform(pdist(align,'hamm')<theta))));
    end
    Meff=sum(W);

    Pij_true = zeros(N,N,q,q);
    Pi_true = zeros(N,q);

    for j=1:M
        for i=1:N
            Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
        end
    end
    Pi_true = Pi_true/Meff;

    for l=1:M
        for i=1:N-1
            for j=i+1:N
                Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
                Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
            end
        end
    end
    Pij_true = Pij_true/Meff;

    scra = eye(q,q);
    for i=1:N
        for alpha=1:q
            for beta=1:q
                Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
            end
        end
    end
end


function [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q)
% adds pseudocount

    Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(N,N,q,q);
    Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(N,q);

    scra = eye(q);

    for i=1:N
        for alpha = 1:q
            for beta = 1:q
               Pij(i,i,alpha,beta) =  (1.-pseudocount_weight)*Pij_true(i,i,alpha,beta) + pseudocount_weight/q*scra(alpha,beta);
            end
        end
    end 
end


function [C,C_self] = Compute_C(Pij,Pi,N,q)
C=zeros(N*(q-1),N*(q-1));
C_self=zeros(N*(q-1),N*(q-1));
for i=1:N
    for j=1:N
        for alpha=1:q-1
            for beta=1:q-1
                 C(mapkey(i,alpha,q),mapkey(j,beta,q)) = Pij(i,j,alpha,beta) - Pi(i,alpha)*Pi(j,beta);
                 if ( i == j )
                     C_self(mapkey(i,alpha,q),mapkey(j,beta,q)) = C(mapkey(i,alpha,q),mapkey(j,beta,q));
                 end
            end
        end
    end
end
end


function A=mapkey(i,alpha,q)
    A = (q-1)*(i-1)+alpha;
end


function W=ReturnW(C,i,j,q)
    W = zeros(q,q);
    W(1:q-1,1:q-1) = C(mapkey(i,1:q-1,q),mapkey(j,1:q-1,q)) ;
end


function F_apc = old_calc_norm_apc( invC, N, q)
    F = zeros(N);
    F_apc = zeros(N);

    for i = 1:(N-1)
        for j = (i+1):N
            J_mf = ReturnW(invC,i,j,q);
            J_j = mean(J_mf,1);
            J_i = mean(J_mf,2);
            J = mean(J_i);
         
            for a = 1:q
                for b = 1:q
                    J_mf(a,b) = J_mf(a,b) - J_i(a) - J_j(b) + J;
                end
            end
            
            F(i,j) = norm( J_mf, 'fro' );
            F(j,i) = F(i,j);
            
        end
    end

    F_i = mean(F);
    F_av = mean(F_i);

    for i = 1:(N-1)
        for j = (i+1):N
            F_apc(i,j) = F(i,j) - F_i(i)*F_i(j)/F_av;
            F_apc(j,i) = F_apc(i,j);
        end
    end
end
