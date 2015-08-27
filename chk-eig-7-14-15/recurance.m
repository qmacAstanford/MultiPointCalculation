function recurance(K,ORDEig,ORDL,NumLayer)


EigK=IntermediateKEigenValues(K,ORDEig,ORDL);
disp('Eigenvalues calculated by matrix methode (L+1,mu+1)')
disp(EigK)
% EigK:  (L+1,mu+1)

Denom=zeros(ORDL,ORDL,ORDEig)*NaN;
for M=0:(ORDL-1) % M is mu
    % pre calculate a, rows refer to lambda+1
    lamVec=0:NumLayer+1;
    a=sqrt((lamVec-M).*(lamVec+M)./(4*lamVec.^2-1));
    
    for L=0:(ORDEig-1) % may be able to vectorize this loop for speed
        % pre calculate P, rows refer to lambda+1
        P=EigK(L+1,M+1)+lamVec.*(lamVec+1);  %  --------- Change this line -----------
        %P=EigK(L+1,1)+lamVec.*(lamVec+1);
        
        % pre calculate \tilde j_\lamda^{\mu (+)}
        % row will refer to \lamda+1
        jp=zeros(NumLayer+1,1)*NaN;
        jp(NumLayer+1)=P(NumLayer+1)/K; % initialize recursion
        for lamPrime=(NumLayer-1):-1:M
            jp(lamPrime+1)=(P(lamPrime+1)/K) + ((a(lamPrime+1+1))^2 / jp(lamPrime+1+1));
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
        
        for lam0=0:(ORDL-1)
            if lam0<M
                Denom(lam0+1,M+1,L+1)=NaN;
            elseif lam0 == M
                Denom(lam0+1,M+1,L+1)=P(lam0+1) + a(lam0+1+1)^2*K*jp(lam0+1+1)^(-1);
            else    
                Denom(lam0+1,M+1,L+1)=a(lam0+1)^2*K*jm(lam0+1-1)^(-1) + P(lam0+1) + a(lam0+1+1)^2*K*jp(lam0+1+1)^(-1);
            end            
        end
    end
end
disp('Denominator (lam+1,mu+1,L+1)')
format long
disp(Denom)