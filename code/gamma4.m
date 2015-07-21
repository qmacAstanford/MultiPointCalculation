%{
Function: gamma4
This function calculates the quartic coefficient
in free energy expansion of Gaussian chain, Wormlike Chain
and Rigid rod

Usage: G=gamma4(Q,Q1,Q2,Q3,Q4,R1,R2,R3,R4,CHAIN,N,NM,LAM,FA)
Return:
    G, Quartic order expansion coeffiecient of free energy
of random copolymer melt

Parameters:
    Q, magnitude of wavevector in unit of 1/contour length
       Q input can be a vector
    Q1, Q2, Q3, Q4, wavevectors. Only relative orientations
       matter, the magnitudes of Qs are all of Q
    R1, R2, R3, R4, replica indices. They can only take values
       of 1 to 4. The choice of Rs is only used to find
       equality relations of replicas
    CHAIN, which type of chain (1=Gaussian, 2=WLC, 3=Rigid Rod)
    N, number of monomers
    NM, number of Kuhn steps per monomer
    FA, fraction of A type monomers
    LAM, chemical correlation of random copolymers
Example: quartic order coefficient of ideal random copolymers FA=0.5
    with all replicas have same indices, and Q's form a square loop
    k=logspace(-1,3,50);
    Q1=[1,0,0];
    Q2=[0,-1,0];
    Q3=[-1,0,0];
    Q4=[0,1,0];
    CHAIN=1;
    N=10;
    NM=10;
    G=gamma4(k,Q1,Q2,Q3,Q4,...
                    1,1,1,1,...
                    CHAIN,N,NM,-0.75,0.5)
    figure;semilogx(k,G)


Note cubic coefficient vanishes unless all followings
are satisfied
-- Q1, Q2, Q3 and Q4 have the same magnitude and sum to zero
    (note Qs dont need to all lie on the same plane)
-- Replica indices are all equal, or two pairs of replica indices
    are equal. e.g. R1==R2 and R3==R4 without any other equality

Shifan Mao 06/10/13
%}

function [G]=gamma4(Q,Q1,Q2,Q3,Q4,R1,R2,R3,R4,CHAIN,N,NM,LAM,FA)
MIN=1e-10;
if (Q1+Q2+Q3+Q4)>=MIN
    disp('ERROR :: Qs must add up to zero')
    stop;
else
    %result to return
    NQ=length(Q);
    G=zeros(NQ,1);
    
    %parameters for WLC calculations
    ORDEig=5;
    ORDL=5;
    NumLayer=500;

    %parameters :: polymer

    %parameters :: wave vector magnitudes
    Q=Q./N./NM; %normalize Q magnitudes to contour length
    d=3;        %calculation in 3-dimension

    %combination matrix
    M1=combinator(2,8,'p','r');
    M2=combinator(2,9,'p','r');

    for J=1:NQ
        Qmag=Q(J);
        G(J)=gam4(Qmag,Q1,Q2,Q3,Q4,R1,R2,R3,R4,N,NM,LAM,FA,d,M1,M2,CHAIN,ORDEig,ORDL,NumLayer);
    end

end
end

function G=gam4(Q,Q1,Q2,Q3,Q4,R1,R2,R3,R4,N,NM,LAM,FA,d,M1,M2,CHAIN,ORDEig,ORDL,NumLayer)
    D=[1,-1];

    %wave vectors
    Q1=Q1./norm(Q1)*Q;
    Q2=Q2./norm(Q2)*Q;
    Q3=Q3./norm(Q3)*Q;
    Q4=Q4./norm(Q4)*Q;
    v=0.1;     % monomer volume
    V=power(20,3);      % system volume
    pref1=power(2*pi,-12)*(V/v)*power(N*NM,3);
    pref2=power(2*pi,-13)*(V/v)*power(N*NM,3)*v^2;
    
    %calculate single chain correlation functions
    if CHAIN==1
        s4 = s4gaussian(N,NM,LAM,FA,Q1,Q2,Q3,Q4,d);
        s3 = s3gaussian(N,NM,LAM,FA,Q1,Q2,Q3,d);
        s2 = s2gaussian(N,NM,LAM,FA,Q,d);
        s2inv = s2invgaussian(N,NM,LAM,FA,Q,d);
    elseif CHAIN==2
        s4 = s4wlc(N,NM,LAM,FA,Q1,Q2,Q3,Q4,ORDEig,ORDL,NumLayer);
        s3 = s3wlc(N,NM,LAM,FA,Q1,Q2,Q3,ORDEig,ORDL,NumLayer);
        s2inv = s2invwlc(N,NM,LAM,FA,Q,d,ORDEig,ORDL,NumLayer);
    elseif CHAIN==3
        s4=s4rigid(N,NM,LAM,FA,Q1,Q2,Q3,Q4,1);
        s3 = s3rigid(N,NM,LAM,FA,Q1,Q2,Q3,1);
        s2 = s2rigid(N,NM,LAM,FA,Q,1);
        s2inv = s2invrigid(N,NM,LAM,FA,Q);
    end
    
    % calculate cumulant average of chain correlation
    if (R1==R2 && R1==R3 && R1==R4)
        % all replicas the same
        for I1=1:2
            for I2=1:2
                for I3=1:2:q
                    for I4=1:2
                        s4(I1,I2,I3,I4) = s4(I1,I2,I3,I4)-3*power(s2(I1,I2),2);
                    end
                end
            end
        end
    elseif (R1==R2 && R3==R4 && R1~=R3)

    elseif (R1==R3 && R2==R4 && R1~=R2)

    elseif (R1==R4 && R2==R3 && R1~=R2)

    end
    
    G=0;
    for I = 1:256
        G = G + pref1*real(...
            s4(M1(I,1),M1(I,2),M1(I,3),M1(I,4))*...
            s2inv(M1(I,1),M1(I,5))*...
            s2inv(M1(I,2),M1(I,6))*...
            s2inv(M1(I,3),M1(I,7))*...
            s2inv(M1(I,4),M1(I,8))*...
            D(M1(I,5))*D(M1(I,6))*D(M1(I,7))*D(M1(I,8)));
    end

    for I = 1:512
        G = G + pref2*real(...
            -3*s3(M2(I,1),M2(I,2),M2(I,5))*...
               s3(M2(I,3),M2(I,4),M2(I,5))*...
               s2inv(M2(I,1),M2(I,6))*...
               s2inv(M2(I,2),M2(I,7))*...
               s2inv(M2(I,3),M2(I,8))*...
               s2inv(M2(I,4),M2(I,9))*...
               D(M2(I,6))*D(M2(I,7))*D(M2(I,8))*D(M2(I,9)));
    end
    
end