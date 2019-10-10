% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function [f,g] = GUP_genfilter(G,type,par)

% Calculates different filter pairs for a graph G

% INPUT:  
% G            : The graph
% type         : 'P-P'         : A projection-projection pair
%              : 'D-L'         : A distance-Laplace pair
%              : 'D-P'         : A distance-projection pair
%              : 'M-D-P'       : A modified distance-projection pair
%
% par          : the relevant filter parameters are
%              : idxcen        : Index of center node (spatial filter)
%              : fcen          : Center node (spatial filter)
%              : frad          : Radius (spatial filter)
%              : gK            : Range of spectrum (spectral filter)
%              : alpha         : Parameter for modified spatial filter
%              : beta          : Parameter for spectral filter
%
% OUTPUT:    
% f            : Spatial filter
% g            : Spectral filter

N = G.N;
f = zeros(N,1);
g = zeros(N,1);

switch type
    
    case {'slepian','projection-projection','P-P'}

         for i = 1 : N
            if norm( G.nodes(i,:) - par.fcen, 2 ) <= par.frad
               f(i) = 1;
            end
         end

         g(1:par.gK) = ones(par.gK,1);
         
     case {'graphspread','distance-laplace','D-L',}
         dist = zeros(N,1);      
         for j = 1:N
            s = zeros(N,1);
            s(par.idxcen) = 1;
            
            while s(j) < 1e-12 && dist(j)<=N+1
                dist(j) = dist(j)+1;
                s = G.A*s;
            end
            
            if dist(j) == N+1
               dist(j) = inf;
            end
         end
         
         f = ones(N,1)-dist.^2/max(dist).^2;
         g = ones(N,1)-G.lambda/2;

     case {'distance-projection','D-P'}
         dist = zeros(N,1);      
         for j = 1:N
            s = zeros(N,1);
            s(par.idxcen) = 1;
            
            while s(j) < 1e-12 && dist(j)<=N+1
                dist(j) = dist(j)+1;
                s = G.A*s;
            end
            
            if dist(j) == N+1
               dist(j) = inf;
            end
         end
         
         f = ones(N,1)-dist/max(dist);
         g(1:par.gK) = ones(par.gK,1);
         
     case {'modified-distance-projection','M-D-P'}
         dist = zeros(N,1);      
         for j = 1:N
            s = zeros(N,1);
            s(par.idxcen) = 1;
            
            while s(j) < 1e-12 && dist(j)<=N+1
                dist(j) = dist(j)+1;
                s = G.A*s;
            end
            
            if dist(j) == N+1
               dist(j) = inf;
            end
         end
         
         f = ones(N,1)-dist.^(par.alpha)/max(dist).^(par.alpha);
         g(1:par.gK) = ones(par.gK,1);
         g = g.*(ones(N,1)-(G.lambda).^(par.beta)/2^(par.beta));
end


end