% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

clear all
close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'bunny';

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
filtpar.frad = 0.015;         %Radius (spatial filter)
filtpar.gK = 200;             %Range of spectrum (spectral filter)
filtpar.alpha = 1;            %Parameter for modified filter
filtpar.beta = 2;             %Parameter for modified filter

%Choose plot parameters
plotpar.P = 400;
plotpar.KK = G.N;
plotpar.MM = 5;
plotpar.MMplus = 10;
plotpar.fontsize = 10;

%-------------------------------------------------------------------


% Plot
figure('Units', 'pixels', ...
'Position', [0 50 1200 300]);

subplot(2,8,[1,2,9,10]),
[f,g] = GUP_genfilter(G,'slepian',filtpar);
plotpar.title = 'P-P filter';
plotpar.color = [1,250/255,150/255];
[rerange1,imrange1] = GUP_drawshape(G,f,g,plotpar);

subplot(2,8,[3,4,11,12]), 
[f,g] = GUP_genfilter(G,'distance-projection',filtpar);
plotpar.title = 'D-P filter';
plotpar.color = [1,200/255,100/255];
[rerange2,imrange2] = GUP_drawshape(G,f,g,plotpar);

subplot(2,8,[5,6,13,14]),
[f,g] = GUP_genfilter(G,'modified-distance-projection',filtpar);
plotpar.title = 'Modified D-P filter';
plotpar.color = [1,150/255,50/255];
[rerange3,imrange3] = GUP_drawshape(G,f,g,plotpar);

subplot(2,8,[7,8,15,16]), 
[f,g] = GUP_genfilter(G,'distance-laplace',filtpar);
plotpar.title = 'D-L filter';
plotpar.color = [1,100/255,0/255];
[rerange4,imrange4] = GUP_drawshape(G,f,g,plotpar);

%-------------------------------------------------------------------

figure
fill(rerange1,imrange1,[1,250/255,150/255]);
hold on
fill(rerange2,imrange2,[1,200/255,100/255]);
hold on
fill(rerange3,imrange3,[1,150/255,50/255]);
hold on
fill(rerange4,imrange4,[1,100/255,0/255]);
hold off
axis equal;
axis([-0.05 1.05 -0.05 1.05]);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 10)
title('Uncertainty region of all 4 filters')