% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function GUP_drawspectrum(G,signal,plotspect)

% Calculates the spectral distribution of a signal on the graph $G$ 

% INPUT:    
% G            : The graph G
% signal       : The signal
% plotspect    : The following parameters are relevant:
%                part        : upper boundary of spectrum to be plotted
%                yhight      : upper bound of maximal value
%                fontsize    : fontsize


y = abs(G.U'*signal);

 if ~exist('plotspect','var')
      plotspect.part = 1/3;
      plotspect.yhight = max(y)+0.01;
      plotspect.fontsize = 18;    
 end
 
 if ~isfield(plotspect,'part')
      plotspect.part = 1/3;
 end
 
 if ~isfield(plotspect,'yhight')
      plotspect.yhight = max(y)+0.01;
 end
 
 if ~isfield(plotspect,'fontsize')
      plotspect.fontsize = 18;
 end

plot(1:G.N,y, 'color',[153, 51, 255]/255, 'LineWidth', 2)
axis([1 G.N*plotspect.part 0 plotspect.yhight]);
axis square;
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', plotspect.fontsize)
end



