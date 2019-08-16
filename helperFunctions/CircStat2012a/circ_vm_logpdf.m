function [p alpha] = circ_vm_logpdf(alpha, thetahat, kappa)

%This code was stolen from here: 

%{
http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
written by :  Francesco Montorsi
%}


% if no angles are supplied, 100 evenly spaced points around the circle are
% chosen
if nargin < 1 || isempty(alpha)
alpha = linspace(0, 2*pi, 101)';
alpha = alpha(1:end-1);
end
if nargin < 3
kappa = 1;
end
if nargin < 2
thetahat = 0;
end

alpha = alpha(:);

% evaluate pdf
C = -log( 2*pi*besseli(0,kappa) );
p = C + kappa*cos(alpha-thetahat);