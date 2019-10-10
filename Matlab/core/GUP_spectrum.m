% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function [U,Lambda] = GUP_spectrum(L,stype)

% Calculates the eigendecomposition of L

% INPUT:    
% L            : A symmetric positive definite matrix
% stype        : 'ascend' or 'descend' ordering of eigenvalues
%
% Output:    
% U            : The eigenvectors of L
% lambda       : The (ordered) eigenvalues of L

[U,Lambda,~]=svd(L);
    
[Lambda,inds] = sort(diag(Lambda),stype);
U = U(:,inds);
    
signs=sign(U(1,:));
signs(signs==0)=1;
U = U*diag(signs);
end