function s2inv=s2invgaussian(N,NM,LAM,FA,k,d)

s2inv=zeros(2,2);

% Calculate the s matrix
% [SAA,SAB,SBA,SBB]=s2gaussian(N,NM,FA,LAM,k,d);
% DET=SAA*SBB-SAB*SBA;

% s2inv(1,1) = SBB./DET;
% s2inv(1,2) = -SAB./DET;
% s2inv(2,1) = -SAB./DET;
% s2inv(2,2) = SAA./DET;

s2=s2gaussian(N,NM,LAM,FA,k,d);
DET=s2(1,1)*s2(2,2)-s2(1,2)*s2(2,1);

s2inv(1,1) = s2(2,2)./DET;
s2inv(1,2) = -s2(1,2)./DET;
s2inv(2,1) = -s2(2,1)./DET;
s2inv(2,2) = s2(1,1)./DET;