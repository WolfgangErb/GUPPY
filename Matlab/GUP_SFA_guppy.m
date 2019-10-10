% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

clear all
close all

% Paths
addpath(genpath('./core/'))
addpath(genpath('./data/'))

%Choose graph
G.type = 'guppy';

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
filtpar.idxcen = 49;          %Index of center node (spatial filter)
filtpar.fcen = G.nodes(49,:); %Center node (spatial filter)
filtpar.frad = 60;            %Radius (spatial filter)
filtpar.gK = 200;             %Range of spectrum (spectral filter)
filtpar.alpha = 1/2;          %Parameter for modified filter
filtpar.beta = 2;             %Parameter for modified filter

%Choose plotting parameter
plotpar.MM = 0.5;             %size of dots
plotpar.lb = 20;              %left boundary
plotpar.ub = 20;              %upper boundary
plotpar.uaxis = 0.5;          %range of values for colorbar
plotpar.fontsize = 12;        %fontsize
plotpar.colorbar = 'n';       %put 'y' if you want a colorbar

plotspect.part = 1/2;         %part of the spectrum to be plotted
plotspect.yhight = 0.69;      %max hight of the spectrum in plot
plotspect.fontsize = 12;      %fontsize

%Eigenfunction of space-frequency operator to plot
kk = 1;
SFO = 'S';
theta = 0.9*pi/2;

%Create all space-frequency decompositions

[f,g] = GUP_genfilter(G,'slepian',filtpar);
[~,V1,sigma1] = GUP_SFA(G.U,f,g,SFO,theta);

[f,g] = GUP_genfilter(G,'distance-projection',filtpar);
[~,V2,sigma2] = GUP_SFA(G.U,f,g,SFO,theta);

[f,g] = GUP_genfilter(G,'modified-distance-projection',filtpar);
[~,V3,sigma3] = GUP_SFA(G.U,f,g,SFO,theta);

[f,g] = GUP_genfilter(G,'distance-laplace',filtpar);
[~,V4,sigma4] = GUP_SFA(G.U,f,g,SFO,theta);

%-------------------------------------------------------------------

%Plot - spectral localization of largest eigenvector

figure('Units', 'pixels', ...
'Position', [0 50 1200 300]);

subplot(2,8,[1,2,9,10]),

GUP_drawspectrum(G,V1(:,kk),plotspect);
title('P-P-filter')

subplot(2,8,[3,4,11,12]), 

GUP_drawspectrum(G,V2(:,kk),plotspect);
title('D-P-filter')

subplot(2,8,[5,6,13,14]),

GUP_drawspectrum(G,V3(:,kk),plotspect);
title('M-D-P-filter')

subplot(2,8,[7,8,15,16]), 

GUP_drawspectrum(G,V4(:,kk),plotspect);
title('D-L-filter')

hold off


%-------------------------------------------------------------------

%Plot - spatial localization of largest eigenvector

figure('Units', 'pixels', ...
'Position', [0 50 1200 300]);

subplot(2,8,[1,2,9,10]),

GUP_drawsignal(G.nodes,G.edges,abs(V1(:,kk)),plotpar);
set(gca,'XTick',[], 'YTick', [])
title('P-P-filter')

subplot(2,8,[3,4,11,12]), 

GUP_drawsignal(G.nodes,G.edges,abs(V2(:,kk)),plotpar)
set(gca,'XTick',[], 'YTick', [])
title('D-P-filter')

subplot(2,8,[5,6,13,14]),

GUP_drawsignal(G.nodes,G.edges,abs(V3(:,kk)),plotpar)
set(gca,'XTick',[], 'YTick', [])
title('M-D-P-filter')

subplot(2,8,[7,8,15,16]), 

GUP_drawsignal(G.nodes,G.edges,abs(V4(:,kk)),plotpar)
set(gca,'XTick',[], 'YTick', [])
title('D-L-filter')

h = colorbar('vert');
set(h,'Position',[0.92 0.2 0.014 0.635])

hold off

% Plot spectra of space-frequency operators

figure
plot(1:G.N,sigma1,'color',[0,0,0],'marker','.','MarkerSize',10);
hold on
plot(1:G.N,sigma2,'color',[1,100/255,0/255],'marker','.','MarkerSize',10);
hold on
plot(1:G.N,sigma3,'color',[1,150/255,50/255],'marker','.','MarkerSize',10);
hold on
plot(1:G.N,sigma4,'color',[1,200/255,100/255],'marker','.','MarkerSize',10);
hold on
plot(1:G.N,sigma1,'color',[0,0,0],'marker','.','MarkerSize',10);
hold off

title('Spectrum of space-frequency operator');
axis square;
axis([0 300 -0.05 1.2]);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 12)
legend('(f_1,g_1)','(f_2,g_2)','(f_3,g_3)','(f_4,g_4)')

