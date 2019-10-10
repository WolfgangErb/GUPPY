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
filtpar.alpha = 1/2;          %Parameter for modified filter
filtpar.beta = 2;             %Parameter for modified filter

%Choose plotting parameter
plotpar.MM = 2;               %size of dots
plotpar.ub = 0.02;            %upper boundary
plotpar.lb = 0.02;            %left boundary

% draw graph
figure
GUP_drawgraph(G.nodes,G.edges,plotpar);
title('The Stanford bunny')

% create and plot second eigenfunction of Laplacian 
kk = 2;

plotpar.uaxis = max(G.U(:,kk));
plotpar.laxis = min(G.U(:,kk));

figure
GUP_drawsignal(G.nodes,G.edges,G.U(:,kk),plotpar);
title('Second eigenvector of Laplacian');
 
% create modified-distance-projection pair (f,g)
[f,g] = GUP_genfilter(G,'M-D-P',filtpar);

% plot spatial filter f
plotpar.uaxis = 1;
plotpar.laxis = 0;

figure
GUP_drawsignal(G.nodes,G.edges,f,plotpar);
title('The spatial filter f (modified distance filter)')

% plot spectral filter g
figure
GUP_drawspectrum(G,G.U*g);
title('The spectral filter g (modified projection filter)')


% Plot uncertainty region of the pair (f,g)
plotpar.P = 200;
plotpar.KK = G.N;
plotpar.fontsize = 12;

figure('Units', 'pixels', ...
'Position', [0 50 600 600]);

plotpar.title = 'Uncertainty region of the filter pair (f,g) (M-D-P filter pair)';
plotpar.color = [1,150/255,50/255];
GUP_drawshape(G,f,g,plotpar);