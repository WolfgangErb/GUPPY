% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

clear all
close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'sensor1';

%All possibilities
%G.type = 'bunny';
%G.type = 'guppy';
%G.type = 'sensor1;
%G.type = 'sensor2';
%G.type = 'star';
%G.type = 'rand';

%Generate graph
[G.nodes,G.edges,G.A] = GUP_gengraph(G.type);

%Calculate the graph Laplacian
G.N = length(G.nodes(:,1));
G.deg = sum(G.A,1);
isD = diag(1./sqrt(G.deg));
G.L = eye(G.N) - isD*G.A*isD;

% calculate GFT
[G.U,G.lambda] = GUP_spectrum(G.L,'ascend');

%Choose filter parameters
filtpar.idxcen = 51;          %Index of center node (spatial filter)
filtpar.fcen = G.nodes(51,:); %Center node (spatial filter)
filtpar.frad = 0.15;          %Radius (spatial filter)
filtpar.gK = 35;              %Range of spectrum (spectral filter)
filtpar.alpha = 1/2;          %Parameter for modified filter
filtpar.beta = 2;             %Parameter for modified filter

%-------------------------------------------------------------------

[f,g] = GUP_genfilter(G,'modified-distance-projection',filtpar);

C = G.U*diag(g)*G.U';
M = diag(f);
color = [1,155/255,0];
P = 7;

[rerange1, imrange1, rerange2, imrange2] = GUP_numrange(M,C,P);

[xbound, ybound] = GUP_numrange(M,C);

%-------------------------------------------------------------------

figure
fill(rerange2,imrange2,[1,250/255,175/255]);
hold on
fill(xbound,ybound,[1,150/255,0/255]);
hold on
fill(rerange1,imrange1,[1,200/255,105/255]);
hold on 
scatter(rerange1,imrange1,30,'MarkerFaceColor',[1,200/255,105/255],'MarkerEdgeColor','k','Linewidth',1);
hold on 
scatter(rerange2,imrange2,30,'MarkerFaceColor',[1,250/255,175/255],'MarkerEdgeColor','k','Linewidth',1);
hold off
axis equal;
axis([-0.25 1.25 -0.25 1.25]);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 14)
title('Approximating the Numerical Range')
set(h,'Xtick',[0,0.5,1])
set(h,'Ytick',[0,0.5,1])

