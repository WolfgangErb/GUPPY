% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function [rerange, imrange] = GUP_drawshape(G,f,g,plotpar)

% Calculates the numerical range of the pair (M_f,C_g) and, if desired,
% the location of the eigenvectors of the space-frequency operator inside
% the numerical range

% INPUT:    
% G            : The graph G
% f            : The spatial filter
% g            : The frequency filter
% plotpar      : The following parameters are relevant:
%                P           : Approximation order of numerical range
%                KK          : width for bandlimiting
%                MM          : size of nodes
%                MMplus      : size of outside ring 
%                fontsize    : fontsize
%                title       : 'title of plot'
%                SFO         : 'S' or 'R' operator
%                theta       : rotation value of 'R' operator
%                color       : color of numerical range
% OUTPUT:
% [rerange,imrange]          : boundary values of the numerical range

C = G.U*diag(g)*G.U';
M = diag(f);

if isfield(plotpar,'SFO')
   [~,V,~] = GUP_SFA(G.U,f,g,plotpar.SFO,plotpar.theta);
   alpha = diag(V(:,1:plotpar.KK)'*M*V(:,1:plotpar.KK));
   beta = diag(V(:,1:plotpar.KK)'*C*V(:,1:plotpar.KK));
end

[rerange, imrange] = GUP_numrange(M,C,plotpar.P);
fill(rerange,imrange,plotpar.color);

if isfield(plotpar,'SFO')
   hold on
   scatter(alpha,beta,plotpar.MM,'k','filled');
   scatter(alpha(1),beta(1),plotpar.MMplus,'k');
   hold off
end

axis equal;
axis([-0.05 1.05 -0.05 1.05]);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', plotpar.fontsize);
if isfield(plotpar,'title')
   title(plotpar.title);
end

end



