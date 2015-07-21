%{
Function: gamma3
This function calculates the cubic coefficient
in free energy expansion of Gaussian chain, Wormlike Chain
and Rigid rod

Usage: G=gamma3(Q,CHAIN,N,NM,LAM,FA)
Return:
    G, Cubic order expansion coeffiecient of free energy
of random copolymer melt
Parameters:
    Q, magnitude of wavevector in unit of 1/contour length
       Q input can be a vector
    CHAIN, which type of chain (1=Gaussian, 2=WLC, 3=Rigid Rod)
    N, number of monomers
    NM, number of Kuhn steps per monomer
    FA, fraction of A type monomers
    LAM, chemical correlation of random copolymers
Example: cubic order coefficient of ideal random copolymer with FA=0.75
    k=logspace(-1,3,50);
    CHAIN=3;
    N=5;
    NM=3;
    G=gamma3(k,CHAIN,N,NM,0,0.75);
    figure;semilogx(k,G)

Note cubic coefficient vanishes unless all followings
are satisfied
-- Q1, Q2 and Q3 have the same magnitude and
    form an equilateral triangle
-- FA, fraction of A type monomers is different from
    0.5. i.e. not approaching critical point FA

Shifan Mao 06/10/13
%}

function [G]=gamma3(Q,CHAIN,N,NM,LAM,FA)
%result to return
NQ=length(Q);  %Number of magnitudes of Qs
G=zeros(NQ,1);

%parameters for WLC calculations
ORDEig=5;
ORDL=5;
NumLayer=500;

%parameters :: wave vector
Q=Q./N./NM; %normalize Q magnitudes to contour length
d=3;        %three dimension solution
v=0.1;     % monomer volume
V=power(20,3);      % system volume
pref=power(2*pi,-9)*(V/v)*power(N*NM,2);

%combination matrix
M=combinator(2,6,'p','r');

for J=1:NQ
    Qmag=Q(J);
    G(J)=gam3(Qmag,N,NM,LAM,FA,d,M,CHAIN,ORDEig,ORDL,NumLayer);
    G(J)=pref*G(J);
end
end

function G=gam3(Qmag,N,NM,LAM,FA,d,M,CHAIN,ORDEig,ORDL,NumLayer)
    %sign indicator
    D=[1,-1];

    %wave vectors
    Q1=Qmag*[1,0,0];
    Q2=transpose(rotz(pi*2/3)*Q1(1:3)');
    Q3=-Q1-Q2;
    
    if CHAIN==1
        s3 = s3gaussian(N,NM,LAM,FA,Q1,Q2,Q3,d);
        s2inv = s2invgaussian(N,NM,LAM,FA,Qmag,d);
    elseif CHAIN==2
        s3 = s3wlc(N,NM,LAM,FA,Q1,Q2,Q3,ORDEig,ORDL,NumLayer);
        s2inv = s2invwlc(N,NM,LAM,FA,Qmag,d,ORDEig,ORDL,NumLayer);
    elseif CHAIN==3
        s3 = s3rigid(N,NM,LAM,FA,Q1,Q2,Q3,1);
        s2inv = s2invrigid(N,NM,LAM,FA,Qmag);
    end
    
    G=0;
    for I = 1:64
        G = G + real(...
                s3(M(I,1),M(I,2),M(I,3))*...
                s2inv(M(I,1),M(I,4))*...
                s2inv(M(I,2),M(I,5))*...
                s2inv(M(I,3),M(I,6))*...
                D(M(I,4))*D(M(I,5))*D(M(I,6)));
    end
end