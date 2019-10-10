% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function [nodes,edges,A] = GUP_gengraph(type)

% Calculates different examples of graphs

% INPUT:    
% type         : 'sensor1'     : A simple sensor network
%              : 'sensor2'     : A larger sensor network
%              : 'bunny'       : The stanford bunny
%              : 'guppy'       : The guppy graph
%              : 'rand'        : A random sensor network
%              : 'star'        : A simple star graph
%
% OUTPUT:    
% nodes        : Nodes of the graph G
% edges        : Edges of the graph G
% A            : Adjacency of the matrix G


switch type
    
case 'sensor1'
    
load data_sensor1.mat

r = 1/6;

[edges,A] = GUP_NN(nodes,r);

case 'sensor2'
    
load data_sensor2.mat

r = 1/6;

[edges,A] = GUP_NN(nodes,r);

case 'bunny'
    
load data_bunny.mat

nodes = [bunny(:,1),bunny(:,2)];
thresh = 0.0025;
N = size(nodes,1);

stp = 0;
for i = 1 : N
    for j = i+1 : N
        if norm( nodes(i,:) - nodes(j,:), 2 ) <= thresh
            stp = stp + 1;
            idx(stp) = j;
        end
    end
end

nodes(idx,:)=[];

r = 0.01;

[edges,A] = GUP_NN(nodes,r);

case 'guppy'
    
load data_guppy.mat

thresh = 5;
N = size(nodes,1);

stp = 0;
for i = 1 : N
    for j = i+1 : N
        if norm( nodes(i,:) - nodes(j,:), 2 ) <= thresh
            stp = stp + 1;
            idx(stp) = j;
        end
    end
end

nodes(idx,:)=[];

r = 25;

[edges,A] = GUP_NN(nodes,r);

  
case 'rand'
    
N = 300;
thresh = 0.02;
shift = [0,0];
D = [1,1];
nodes = ones(N,1)*D.*rand(N,2)+ones(N,1)*shift;
stp = 0;
for i = 1 : N
    for j = i+1 : N
        if norm( nodes(i,:) - nodes(j,:), 2 ) <= thresh
            stp = stp + 1;
            idx(stp) = j;
        end
    end
end

nodes(idx,:)=[];

r = 1/6;

[edges,A] = GUP_NN(nodes,r);

case 'star'
    
N = 40;

nodes = zeros(N,2);
edges = ones(N-1,2);

nodes(2:N,:) = [cos((0:N-2)'*2*pi/(N-1)),sin((0:N-2)'*2*pi/(N-1))];
edges(:,2) = (2:N)';


A = GUP_adjmat(edges,N);

end
