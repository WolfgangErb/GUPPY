% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function A = GUP_adjmat(edges,N)

% This function calculates the adjacency matrix of an unweighted 
% undirected graph from the information of the edges

% INPUT:    
% edges        : The edges of the graph G
% N            : Number of vertices of the graph G
%
% OUTPUT:  
% A            : Adjacency matrix of the graph G

A = zeros(N,N);

for i = 1:size(edges,1)
    A(edges(i,1),edges(i,2)) = 1;
    A(edges(i,2),edges(i,1)) = 1;
end
