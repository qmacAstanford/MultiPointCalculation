%{
Function : gamma2
Usage: G=gamma2(Q,CHAIN,N,NM,LAM,FA,CHI)
Return:
    G, Quadratic order expansion coeffiecient of free energy
of random copolymer melt
Parameters:
    Q, magnitude of wavevector in unit of 1/contour length
       Q input can be a vector
    CHAIN, which type of chain (1=Gaussian, 2=WLC, 3=Rigid Rod)
    N, number of monomers
    NM, number of Kuhn steps per monomer
    FA, fraction of A type monomers
    LAM, chemical correlation of random copolymers
    CHI, chemical incompatibility between A and B monomers, non-
        dimensionalized by monomer volume v (CHI*v)
Example:
    k=logspace(-1,4,50);
    CHAIN=3;
    N=10;
    NM=1;
    CHI=0./NM;
    G=gamma2(k,CHAIN,N,NM,-0.75,0.5,CHI);
    figure;semilogx(k,1./G)

This function calculates the quadratic coefficient
in free energy expansion of Gaussian chain, Wormlike Chain
and Rigid rod

Shifan Mao 06/10/13
%}

function [G]=gamma2(Q,CHAIN,N,NM,LAM,FA,CHI)
ORDEig=20;
ORDL=20;
NumLayer=500;

%parameters :: wave vector
Q=Q./N./NM;
d=3;        % three dimension solution
v=0.1;     % monomer volume
pref=power(2*pi,-3)*(1/v)*N*NM;
CHI=CHI./(N*NM); %from factoring prefactor

NQ=length(Q);  %Number of magnitudes of Qs
Qmag = Q; % Qmag=logspace(-1,3,NQ)/N/NM;
G=zeros(NQ,1);

%combination matrix
M=combinator(2,2,'p','r');

for J=1:NQ
    Q=Qmag(J);

    G(J)=gam2(Q,N,NM,LAM,FA,d,M,CHAIN,ORDEig,ORDL,NumLayer)-2*CHI;
    G(J)=pref*G(J);
end

end

function G=gam2(Q,N,NM,LAM,FA,d,M,CHAIN,ORDEig,ORDL,NumLayer)
    %sign indicator
    D=[1,-1];
    
    if CHAIN==1
        %Gaussian chain
        s2inv = s2invgaussian(N,NM,LAM,FA,Q,d);
    elseif CHAIN==2
        %Worm-like chain
        s2inv = s2invwlc(N,NM,LAM,FA,Q,d,ORDEig,ORDL,NumLayer);
    elseif CHAIN==3
        %Rigid rod
        s2inv = s2invrigid(N,NM,LAM,FA,Q);
    end
    
    G=0;
    for I = 1:4
        G = G + real(...
                s2inv(M(I,1),M(I,2))*...
                D(M(I,1))*D(M(I,2)));
    end
end