function testYlm()

npts=100;
Lmax=20;
amp=zeros(Lmax+1,1);
Data=zeros(npts,1);
theta=linspace(0,pi(),npts);
for L=0:Lmax
    amp(L+1)=ylm(L,0,0,0);
    for j=1:npts
        Data(j)=Data(j)+amp(L+1)*ylm(L,0,theta(j),0);
    end
    
end
plot(theta,Data,'o-')
xlabel('theta')
ylabel('probability')
title(sprintf('Trying to fit delta function at theta=0 with L=0,1...%g',Lmax)) 

figure
plot(1:Lmax+1,amp,'.')
xlabel('L')
ylabel('amplitude')
title('For rigid rod the coefficent increases with L!')

if 0
L=1;
M=0;

% first test ylm at one point
ylm(L,M,0.56,1.23)

% now plot ylm
npts=20;
Data=zeros(npts,2*npts);
theta=linspace(0,pi(),npts);
phi=linspace(0,2*pi(),2*npts);

for i1=1:npts
    for i2=1:2*npts
        Data(i1,i2)=real(ylm(L,M,theta(i1),phi(i2)));
    end
end

contour(phi,theta,Data)
colorbar
xlabel('phi')
ylabel('theta')
title(sprintf('Ylm L=%g, M=%g',L,M))
format long


% Calculate D matrix
alpha=1;
beta=1;
format short
L_max=10;
D=WignerD_lm0(L_max,alpha,beta)
% Rotate plot
Data=zeros(npts,2*npts);

for i1=1:npts
    for i2=1:2*npts
        YLM=ylmVec(L,theta(i1),phi(i2));
        %Data(i1,i2)=YLM(M+1);
        Data(i1,i2)=YLM(1)*D(L+1,1) + 2*real( dot(D(L+1,2:L+1),YLM(2:end)) ); 
    end
end

figure
contour(phi,theta,Data)
xlabel('phi')
ylabel('theta')
title(sprintf('Ylm L=%g, M=%g, rotated by alpha=%g, beta=%g',L,M,alpha,beta))
colorbar
format long

% check meaning of plot
%{
for i1=1:npts
    for i2=1:2*npts
        Data(i1,i2)=theta(i1)+phi(i2);
    end
end

figure
contour(phi,theta,Data)
colorbar
xlabel('phi')
ylabel('theta')
title(sprintf('theta + phi '))
%}


Q1=[1;2;3];
Q2=[3,2,1];
Q3=[1,1,1];

Q1_n=Q1/norm(Q1,2);
Q2_n=Q2/norm(Q2,2);
Q3_n=Q3/norm(Q3,2);

aligned=10^-13; % angle in radians betwene to vectors before I assume they are the same
if (Q2_n-Q1_n) < aligned
    cosB1=0;
    alpha1=0;
    cosB2=dot(Q2_n,Q3_n);
elseif (Q2_n-Q1_n) < aligned
    cosB1=dot(Q2_n,Q1_n);
    alpha1=0;
    cosB2=0;
else
    cosB1=dot(Q2_n,Q1_n);
    cosB2=dot(Q2_n,Q3_n);
    v1=cross(Q2_n,Q1_n)/norm(cross(Q2_n,Q1_n));
    v2=cross(Q2_n,Q3_n)/norm(cross(Q2_n,Q3_n));
    alpha=acos(dot(v1,v2));
end 
firstD=WignerD_lm0(L_max,alpha1,acos(cosB1));
secondD=WingerD_lm0(L_max,alpha2,acos(cosB2));

for L1 =0:L_max
    for L2 = 0:L_max
        for M=1
        end
    end
end


end
end
function out = ylm(l,m,theta,phi)

P=legendre(l,cos(theta));

out=sqrt((2*l+1)*factorial(l-m)/(4*pi()*factorial(l+m)))*P(m+1)*exp(1i*m*phi);

end
function out = ylmVec(l,theta,phi)
    m=(0:l)';
    P=legendre(l,cos(theta));
    out=sqrt((2*l+1)*factorial(l-m)./(4*pi()*factorial(l+m)))...
        .*P.*exp(1i*m*phi);
end
function D = WignerD_lm0(L_max,alpha,beta)
%  Returns the Winger D matrix, D_l,m,0 of Euler rotations alpha, beta
%  Rows: l=0,1...L_max
%  Cols: m=0,1...l
%  Euler angle convention: z-y'-z''
    alpha=-alpha; beta=-beta;  %chage passive to active
    D=zeros(L_max+1,L_max+1)*NaN;
    for l=0:L_max
        D(l+1,1:l+1)=conj( sqrt(4*pi()/(2*l+1) )* ylmVec(l,beta,alpha));
    end
end