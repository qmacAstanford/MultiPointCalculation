function GLMK=glm(K,EigK,ORDEig,ORDL,NumLayer)
% output GLM(lamda1+1,lamda2+1,mu+1,l+1)
% this function does not contain the small K limit
% K is the magnitude
% EigK: eigenvalues associated with K.  This is a ORDEig x ORDL matrix.
GLMK=zeros(ORDL,ORDL,ORDL,ORDEig);
MIN=1e-10;

if abs(K)<MIN
    for M=0:(ORDL-1)
    for L1=M:(ORDL-1)
        GLMK(L1+1,L1+1,M+1,:)=1;
    end
    end

else

for M=0:(ORDL-1)
    AL=zeros(NumLayer,1);
    for iEig=1:ORDEig
        
        WP=zeros(NumLayer,1);
        WM=zeros(NumLayer,1);
        dJp=zeros(NumLayer,1);
        dJm=zeros(NumLayer,1);
        WPPROD=zeros(NumLayer,1);
        
        LP=NumLayer-1;
        AL(NumLayer)=sqrt((LP-M)*(LP+M))/sqrt(4*LP^2-1);
        WP(NumLayer)=1/(EigK(iEig,M+1)+LP*(LP+1)); % Q.J.M. changed 8/3/15
        
        LM=M;
        if EigK(iEig,M+1)+LM*(LM+1)==0 % Q.J.M. changed 8/3/15
            WM(M+1)=0;
        else
            WM(M+1)=1/(EigK(iEig,M+1)+LM*(LM+1)); % Q.J.M. changed 8/3/15
        end
        AL(M+1)=sqrt((LM-M)*(LM+M))/sqrt(4*LM^2-1);
        
        dJp(NumLayer)=1;
        dJm(M+1)=1;
        WPPROD(NumLayer)=1i*K*AL(NumLayer)*WP(NumLayer);
        
        for n=(NumLayer-2):-1:M
            IP=n+1;
            LP=IP-1;
            IM=NumLayer-n+M;
            LM=IM-1;
            
            PL=EigK(iEig,M+1)+LP*(LP+1); % Q.J.M. changed 8/3/15
            AL(IP)=sqrt((LP-M)*(LP+M))/sqrt(4*LP^2-1);
            WP(IP)=1/(PL+WP(IP+1)*(AL(IP+1)*K)^2);
            dJp(IP)=1-dJp(IP+1)*(AL(IP+1)*K*WP(IP+1))^2;
            
            WPPROD(IP)=1i*K*AL(IP)*WP(IP);
            
            PL=EigK(iEig,M+1)+LM*(LM+1); % Q.J.M. changed 8/3/15
            AL(IM)=sqrt((LM-M)*(LM+M))/sqrt(4*LM^2-1);
            if WM(IM-1)==0
                WM(IM)=0;
            else
                WM(IM)=1/(PL+WM(IM-1)*(AL(IM)*K)^2);
            end
            dJm(IM)=1-dJm(IM-1)*(AL(IM)*K*WM(IM-1))^2;
            
        end
        
        for L1=M:(ORDL-1)
            for L2=L1:(ORDL-1)
                L=sort([L1 L2]);
                IM=L(1)+1;
                IP=L(2)+1;
                if L1==M
                    dWL0M=1/dJp(IM);
                else
                    dWL0M=1/(-(AL(IM)*K*WM(IM-1))^2*dJm(IM-1)+dJp(IM));
                end

                if IM==IP
                    GLMK(L1+1,L2+1,M+1,iEig)=dWL0M;
                else
                    GLMK(L1+1,L2+1,M+1,iEig)=dWL0M*prod(WPPROD((IM+1):IP));
                end
                GLMK(L2+1,L1+1,M+1,iEig)=GLMK(L1+1,L2+1,M+1,iEig);
%                  if M==2 && L1==3 && L2 ==3 && iEig==1 
%                      disp('dJm old')
%                      disp(dJm(1:(ORDL+2))/K)
%                      disp('dJp old')
%                      disp(dJp(1:(ORDL+2))/K)  
%                      disp('a old')
%                      disp(AL(1:ORDL+2)')    
%                      sprintf('old M=%g, K=%g, L=%g, GLMK=%g',M,K, iEig-1, GLMK(L1+1,L2+1,M+1,iEig))
%                  end                
            end
        end
        
    end
end
end
end