% GUPPY: a very simple toolbox for
% space-frequency decompositions and uncertainty principles on graphs
% (C) W. Erb 01.08.2019

function [xin, yin, xout, yout] = GUP_numrange(A,B,MM)

% This method approximates the numerical range of the 
% matrix pair (A,B) from inside and from outside with a
% polygon having MM vertices.

% Reference: this method is described in the following paper:
% Johnson, C.R. Numerical Determination of the Field of Values of a
% General Complex Matrix. SIAM J. Num. Anal. 15 , 3 (1978), 595–602.

% INPUT:    
% A,B          : pair of real symmetric matrices
% MM           : number of vertices of the polygonal approximants
%
% OUTPUT:  
% xin,yin      : x and y-coordinates of the nodes of the inner polygon
% xout,yout    : x and y-coordinates of the nodes of the outer polygon

if nargin == 2
    MM = 400;
end

CosTH = cos((0:MM-1)'*2*pi/MM);
SinTH = sin((0:MM-1)'*2*pi/MM);

xin = zeros(MM,1);
yin = zeros(MM,1);
sigma = zeros(MM,1);

% Calculation of inside polygon
for th=0:MM-1
    H=CosTH(th+1)*A + SinTH(th+1)*B;
    SymH=.5*(H+H'); 
    [U,Lambda]=eig(SymH);
    [eigth,idx] = sort(diag(Lambda),'descend');
    u = U(:,idx(1)); 
    sigma(th+1) = eigth(1);
    xin(th+1)=u'*A*u;
    yin(th+1)=u'*B*u;
end

% Calculation of outside polygon
sigmaup = (sigma*CosTH(2)-circshift(sigma,1))/SinTH(2);
xout = CosTH.*sigma-SinTH.*sigmaup;
yout = SinTH.*sigma+CosTH.*sigmaup;

end