% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function GUP_drawsignal(nodes,edges,signal,plotpar)

% This function draws a signal on the nodes of the graph G

% INPUT:    
% nodes        : Nodes of the graph G
% edges        : Edges of the graph G
% signal       : The signal on the Nodes of G
% plotpar      : The following parameters are relevant:
%                MM          : size of nodes
%                lb          : additional space for left and right boundary
%                ub          : additional space for upper and lower boundary
%                uaxis       : upper bound for representation of values
%                laxis       : lower bound for representation of values
%                fontsize    : fontsize
%                colorbar    : 'y' or 'n'

 if ~exist('plotpar','var')
      plotpar.MM = 4;
      plotpar.lb = 0.01;
      plotpar.ub = 0.01;
      plotpar.uaxis = 1/3;
      plotpar.laxis = 0;    
 end
 
 if ~isfield(plotpar,'MM')
      plotpar.MM = 4;
 end
 
 if ~isfield(plotpar,'lb')
     plotpar.lb = 0.01;
 end
 
 if ~isfield(plotpar,'ub')
     plotpar.ub = 0.01;
 end
 
 if ~isfield(plotpar,'uaxis')
     plotpar.uaxis = 1/3;
 end
 
  if ~isfield(plotpar,'laxis')
     plotpar.laxis = 0;
  end
  
  if ~isfield(plotpar,'fontsize')
     plotpar.fontsize = 18;
  end
  
  if ~isfield(plotpar,'colorbar')
     plotpar.colorbar = 'y';
  end
 
% plot the graph with nodes and edges

% draw the edges
for k = 1 : size(edges,1)
    i = edges(k,1);
    j = edges(k,2);    
    plot( [nodes(i,1) nodes(j,1)], [nodes(i,2) nodes(j,2)], 'color',[0.6,0.6,0.6], 'LineWidth', 1 );
    hold on
end
% draw the signal on the nodes
scatter(nodes(:,1),nodes(:,2), plotpar.MM*100, signal, '.');
hold off;
caxis([plotpar.laxis plotpar.uaxis] ) ;
axis square;
axis([min(nodes(:,1))-plotpar.lb max(nodes(:,1))+plotpar.lb min(nodes(:,2))-plotpar.ub max(nodes(:,2))+plotpar.ub]);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', plotpar.fontsize)
colormap(flipud(copper));
if plotpar.colorbar == 'y'
   colorbar;
end

end


