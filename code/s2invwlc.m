function s2inv=s2invwlc(N,NM,LAM,FA,k,d,ORDmax,ORD,ResLayer)
% function [val]=s2invwlc(N,NM,FA,LAM,k,d,ORDmax,ORD,ResLayer)

% Calculate the Fourier transform of the Green function
% for the wormlike chain in d-dimensions
%
% Andrew Spakowitz (4/14/15)

% Fill in unset optional values

switch nargin
    case 2
        d=3;
        ORDmax=20;
        ORD=20;
        ResLayer=500;
    case 3
        ORDmax=20;
        ORD=20;
        ResLayer=500;        
    case 4
        ORD=20;
        ResLayer=500;        
    case 5
        ResLayer=500;        
end

% If dimensions is 2, reset value to small perturbation above 2

if d==2
    d=2+1e-10;
end

% Calculate the s matrix
% 
% [SAA,SAB,SBA,SBB]=s2wlc(N,NM,FA,LAM,k,d,ORDmax,ORD,ResLayer);
% 
% DET=SAA.*SBB-SAB.*SBA;
% 
% % val=real(N*NM*(SAA+SBB+2*SAB)./DET);
% s2inv(1,1) = SBB./DET;
% s2inv(1,2) = -SAB./DET;
% s2inv(2,1) = -SAB./DET;
% s2inv(2,2) = SAA./DET;

s2=s2wlc(N,NM,LAM,FA,k,d,ORDmax,ORD,ResLayer);
DET=s2(1,1)*s2(2,2)-s2(1,2)*s2(2,1);

s2inv(1,1) = s2(2,2)./DET;
s2inv(1,2) = -s2(1,2)./DET;
s2inv(2,1) = -s2(2,1)./DET;
s2inv(2,2) = s2(1,1)./DET;
