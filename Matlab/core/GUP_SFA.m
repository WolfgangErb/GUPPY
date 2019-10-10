% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function [S,V,sigma] = GUP_SFA(U,f,g,SFO,theta)

% Calculates the space frequency operator S (or R) and its eigendecomposition 
% from the filters f, g and the Fourier transform U. 

% INPUT:    
% U            : The Fourier transform matrix of the graph G
% f            : The spatial filter
% g            : The frequency filter
% SFO          : 'S' or 'R' space-frequency operator
% theta        : rotation value of 'R' operator
%
% Output:    
% S            : The space-frequency operator S (or R)
% V            : The eigenvectors of S
% sigma        : The (ordered) eigenvalues of S

 if ~exist('SFO','var')
     SFO = 'S';
 end
 
 if ~exist('theta','var')
     theta = pi/4;
 end

C = U*diag(sqrt(g))*U';
M = diag(f);

 if SFO == 'S'
     S = C*M*C;  
 elseif SFO == 'R'
     S = cos(theta)*M + sin(theta)*C;
 end
 
[V,sigma] = GUP_spectrum(S,'descend');
     
end