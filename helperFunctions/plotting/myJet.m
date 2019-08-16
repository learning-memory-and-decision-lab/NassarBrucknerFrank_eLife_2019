function cMap = myJet(m)
%
% FUNCTION:
%   myJet.m
% 
% DESCRIPTION:
%   like jet but from pure blue [0 0 1] to pure red [1 0 0]
%
% INPUTS:
%   m = number of colors
%
% OUTPUTS:
%   cMap = the color map
%
% NOTES:
%   (1) written by jfburke on 04/2012 (john.fred.burke@gmail.com)
%
%

if nargin<1||(~exist('m')||isempty(m))
  m=64;
end

if m==1
  cMap = [0.5 1 0.5];
  return
end

if m==2
  cMap = [1 0 0; 0 0 1];
  return
end

if m==3
  cMap = [1 0 0;  0.5 1 0.5; 0 0 1];
  return
end

if m==5
  cMap = [1 0 0;  0 1 1; 0.5 1 0.5; 1 1 0; 0 0 1];
  return
end


K  = floor(m./3);
N1 = K;
N2 = m-K-N1-1;
N3 = m-2*K-2;

RAMP_up   = [(0:1:K)/K]';
RAMP_down = [(K:-1:0)/K]';
RED   = [zeros(N1,1); RAMP_up; ones(N2,1);];
GREEN = [RAMP_up; ones(N3,1); RAMP_down];
BLUE  = [ones(N2,1); RAMP_down; zeros(N1,1)];

cMap = [RED GREEN BLUE];