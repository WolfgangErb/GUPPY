% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function GUP_drawgraph(nodes,edges,plotpar)

% This function draws the graph based on the nodes and the edges of G

% INPUT:    
% nodes        : Nodes of the graph G
% edges        : Edges of the graph G
% plotpar      : The following parameters are relevant:
%                MM          : markersize of nodes
%                lb          : additional space for left and right boundary
%                ub          : additional space for upper and lower boundary


 if ~exist('plotpar','var')
      plotpar.MM = 5;
      plotpar.lb = 0.01;
      plotpar.ub = 0.01;
 end
 
  if ~isfield(plotpar,'MM')
      plotpar.MM = 5;
 end
 
 if ~isfield(plotpar,'lb')
     plotpar.lb = 0.01;
 end
 
 if ~isfield(plotpar,'ub')
     plotpar.ub = 0.01;
 end
 
% plot the graph with edges and nodes

% draw the edges
for k = 1 : size(edges,1)
    i = edges(k,1);
    j = edges(k,2);    
    plot( [nodes(i,1) nodes(j,1)], [nodes(i,2) nodes(j,2)], 'color',[0.6,0.6,0.6], 'LineWidth', 1 );
    hold on
end
% draw the nodes
plot( nodes(:,1),nodes(:,2),'o','color',[153, 51, 255]/255,'LineWidth',plotpar.MM,'MarkerSize',plotpar.MM)
hold off;
axis([min(nodes(:,1))-plotpar.lb max(nodes(:,1))+plotpar.lb min(nodes(:,2))-plotpar.ub max(nodes(:,2))+plotpar.ub]) ;
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)  
end

