function GLMK=glmV2(K,EigK,ORDEig,ORDL,NumLayer)
% output GLM(lamda1+1,lamda2+1,mu+1,l+1)
% this function does not contain the small K limit
% K is the magnitude
% EigK: eigenvalues associated with K.  This is a ORDEig x ORDL matrix
% where the column denotes mu+1 value

if NumLayer-3 < ORDL
    error('make NumLayer bigger')
end

GLMK=zeros(ORDL,ORDL,ORDL,ORDEig)*NaN; % lam1+1, lam2+1, mu+1, L+1 
Denoms=zeros(ORDL,ORDL,ORDEig)*NaN;

for M=0:(ORDL-1) % M is mu
    % pre calculate a, rows refer to lambda+1
    lamVec=0:NumLayer+1;
    a=sqrt((lamVec-M).*(lamVec+M)./(4*lamVec.^2-1));
    
    for L=0:(ORDEig-1) % may be able to vectorize this loop for speed
        % pre calculate P, rows refer to lambda+1
        P=EigK(L+1,M+1)+lamVec.*(lamVec+1);
        
        % pre calculate \tilde j_\lamda^{\mu (+)}
        % row will refer to \lamda+1
        jp=zeros(NumLayer+1,1)*NaN;
        jp(NumLayer+1)=P(NumLayer+1)/K; % initialize recursion
        for lamPrime=(NumLayer-1):-1:M
            jp(lamPrime+1)=(P(lamPrime+1)/K) + ((a(lamPrime+1+1))^2 / jp(lamPrime+1+1));
%             if M==3 && lamPrime==5
%                 sprintf('and here we have %g',jp(lamPrime+1))
%                 sprintf('second term %g',((a(lamPrime+1+1))^2 / jp(lamPrime+1+1)))
%                 sprintf('first term %g',(P(lamPrime+1)/K) )
%             end
        end
        
        % pre calculate \tilde j_\lamda^{\mu (-)} 
        % row will refer to \lamda+1
        jm=zeros(ORDL+2,1)*NaN; % only need jm up to lam0+1
        jm(M+1)=P(M+1)/K; % initialize recursion
        for lamPrime=(M+1):(ORDL+1)
            jm(lamPrime+1)=P(lamPrime+1)/K + (a(lamPrime+1))^2 / jm(lamPrime+1-1);
        end
        
        % pre calculate \partial_p \tilde j_\lamda^{\mu (+)}
        % row will refer to \lamda+1
        djp=zeros(NumLayer+1,1)*NaN;
        djp(NumLayer+1)=1/K;
        for lamPrime=(NumLayer-1):-1:M
            djp(lamPrime+1)=1/K - (a(lamPrime+1+1))^2 * djp(lamPrime+1+1) ...
                            /(jp(lamPrime+1+1)^2);
        end
        
        % pre calculate \tilde j_\lamda^{\mu (-)} 
        % row will refer to \lamda+1
        djm=zeros(ORDL+2,1)*NaN; % only need jm up to lam0+1
        djm(M+1)=1/K; % initialize recursion
        for lamPrime=(M+1):(ORDL+1)
            djm(lamPrime+1)=1/K - (a(lamPrime+1))^2 * djm(lamPrime+1-1) ...
                           / (jm(lamPrime+1-1)^2);
        end
        
        for lam0=M:(ORDL-1)
            if lam0==M % then the first term is zero
                X1prime=( 0 + 1  ...
                          -a(lam0+2)^2*K*djp(lam0+1+1)/((jp(lam0+1+1))^2) );
                if a(lam0+1) ~= 0
                    error('should be zero')
                end
            else
                X1prime=( -a(lam0+1)^2*K*djm(lam0+1-1)/((jm(lam0+1-1))^2) +1  ...
                          -a(lam0+2)^2*K*djp(lam0+1+1)/((jp(lam0+1+1))^2) );
            end
            for lam=M:(ORDL-1)
                if lam==lam0
                    G=1/X1prime;
                elseif lam>lam0
                    % ranges from lam0+1 up to lam
                    X3=prod(jm((lam0+1+1):(lam+1)));
                    G=1i^(lam-lam0)*prod(a((lam0+1+1):(lam+1)))/(X1prime*X3);
                else  % lam0>lam
                    % ranges from lam0 -1 down to and including lam
                    X2=prod(jm((lam+1):(lam0+1-1)));
                    G=1i^(lam0-lam)*prod(a((lam+1):(lam0+1-1)))/(X1prime*X2);
                end
                GLMK(lam0+1,lam+1,M+1,L+1)=G;
                if M==2 && lam==3 && lam0==3 && L==1 && 0
%                      disp('djm quinns')
%                      disp(djm)
%                      disp('djp quinns')
%                      disp(djp(1:(ORDL+2)))
%                      disp('a quinns')
%                      disp(a(1:ORDL+2))
                     sprintf('Quinns version M=%g, lam0=%g, lam=%g, K=%g, L=%g, G=%g',M,lam0,lam,K, L ,G)
%                      format longEng
%                      disp('P,jp,djp,jm,djm')
%                      disp([P(1:7)',jp(1:7),djp(1:7),jm(1:7),djm(1:7)])
%                      disp('first term')
%                      disp(-a(lam0+1)^2*K*djm(lam0+1-1)/((jm(lam0+1-1))^2))
%                      disp('second term')
%                      disp(-a(lam0+2)^2*K*djp(lam0+1+1)/((jp(lam0+1+1))^2) )
                      Denom=a(lam0+1)^2*K*jm(lam0+1-1)^(-1) + P(lam0+1) + a(lam0+1+1)^2*K*jp(lam0+1+1)^(-1);
                      sprintf('The denominator at the pole is %g',Denom)
                      disp('---------------------------------------------')
                end
                if lam0 == M
                    Denom=P(lam0+1) + a(lam0+1+1)^2*K*jp(lam0+1+1)^(-1);
                else    
                    Denom=a(lam0+1)^2*K*jm(lam0+1-1)^(-1) + P(lam0+1) + a(lam0+1+1)^2*K*jp(lam0+1+1)^(-1);
                end

                  Denoms(lam0+1,M+1,L+1)=Denom;
                  
%                 if abs(Denom) > 10^-5
%                     sprintf('M=%g, lam0=%g, lam=%g, K=%g, L=%g, G=%g',M,lam0,lam,K, L ,G)
%                     sprintf('Denominator =%g',Denom)
%                     error('bad pole')
%                 end
            end
        end
    end
end
end
