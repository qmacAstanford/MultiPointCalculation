function testing()
clc
    LAM=0.8;
    N=5;
    NM=1.4;
    R1=1.283;
    R12=1.43;
    R4=1.572;
    F=[0.7,0.3];
    
    Z1=exp(R1*NM);
    Z12=exp(R12*NM);
    Z3=exp(R4*NM);
    Z1L=Z1*LAM;
    Z12L=Z12*LAM;
    Z3L=Z3*LAM;

    sprintf('old')
    S4=zeros(2,2,2,2);
    DF=[1,-1];
    SDEL=ones(2,2,2,2);
    case8(S4,N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,DF,SDEL,1,2,3,4,0)
    
    sprintf('New')
    S4=zeros(2,2,2,2);
    case8v2(S4,N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,[1,2,3,4])
end


function S4=case8(S4,N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,DF,SDEL,A1,A2,A3,A4,~)
    if N<4
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];
    MIN=(10^-4)/NM;
    if max(abs([E1,E12,E3]))<MIN
       valeq=NM^4;
    elseif max(abs([E1,E12]))<MIN
        valeq=(NM^2/(2*E3^3))*(2*(E12+E3)*(expl(2,NM*E3)+expl(2,-NM*E3))...
            -2*NM*E12*E3*sinh(NM*E3));
    elseif max(abs([E12,E3]))<MIN
        valeq=(NM^2/(2*E1^3))*(2*(E12+E1)*(expl(2,NM*E1)+expl(2,-NM*E1))...
            -2*NM*E12*E1*sinh(NM*E1));
    elseif max(abs([E1,E3]))<MIN
        valeq=(NM^2/(E12^3))*(4*(E1+E12+E3)*sinh(0.5*NM*E12)^2 ...
            -NM*(E1+E3)*E12*sinh(NM*E12));
    elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=2.*E1.^(-2).*NM.^2.*((-1)+cosh(E1.*NM));
    elseif (abs(E1-E12)<MIN && abs(E12-E3)>=MIN)
       valeq=exp(1).^((-1).*(E1+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
          E1.*NM)+(-1).*exp(1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1) ...
          .*(E1+(-1).*E3).^(-1).*E3.^(-1).*NM;
    elseif (abs(E1-E3)<MIN && abs(E1-E12)>=MIN)
       valeq=16.*E1.^(-2).*(E1+(-1).*E12).^(-2).*sinh((1/2).*E1.*NM).^2.*sinh(( ...
        1/2).*(E1+(-1).*E12).*NM).^2;
    elseif (abs(E12-E3)<MIN && abs(E1-E3)>=MIN)
       valeq=exp(1).^((-1).*(E1+E12).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
          E1.*NM)+(-1).*exp(1).^(E12.*NM)).*((-1)+exp(1).^(E12.*NM)).*E1.^( ...
          -1).*(E1+(-1).*E12).^(-1).*E12.^(-1).*NM;
    else
%        valeq=exp(1).^((-1).*(E1+E12+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1) ...
%           .^(E1.*NM)+(-1).*exp(1).^(E12.*NM)).*(exp(1).^(E12.*NM)+(-1).*exp( ...
%           1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1).*(E1+(-1).*E12) ...
%           .^(-1).*(E12+(-1).*E3).^(-1).*E3.^(-1);    
        valeq=2*( -coshl(4,NM*E1)+coshl(4,NM*E12)-coshl(4,NM*E3)-coshl(4,NM*(E1-E12))+...
           coshl(4,NM*(E1-E3))-coshl(4,NM*(E12-E3))+coshl(4,NM*(-E1+E12-E3)) )...
           *( (E12-E3)*E3*E1*(E12-E1) )^(-1);    
    end
    
    % on different monomers
    valne=zeros(2,2,2);
    for I1=1:2
        for I2=1:2
            for I3=1:2
                ZT1=ZE1(I1);
                ZT12=ZE12(I2);
                ZT3=ZE3(I3);
                valne(I1,I2,I3)=case8JPart(ZT1,ZT12,ZT3,N);
                
            end
        end
    end

    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2
                    FV=[F(I1),F(I2),F(I3),F(I4)];
                    DFV=[DF(I1),DF(I2),DF(I3),DF(I4)];
                    
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      valeq*(FV(A1)*FV(A2)*FV(A3)*FV(A4)*valne(1,1,1)+...
                           FV(A1)*FV(A2)*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(1,1,2)+...
                           FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*FV(A4)*valne(1,2,1)+...
                           FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(1,2,2)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*FV(A4)*valne(2,1,1)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(2,1,2)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*FV(A4)*valne(2,2,1)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(2,2,2));
                    
                end
            end
        end
    end
end

function S4=case8v2(S4,N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,order)
    if N<4
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];
    MIN=(10^-4)/NM;
    if max(abs([E1,E12,E3]))<MIN
       valeq=NM^4;
    elseif max(abs([E1,E12]))<MIN
        valeq=(NM^2/(2*E3^3))*(2*(E12+E3)*(expl(2,NM*E3)+expl(2,-NM*E3))...
            -2*NM*E12*E3*sinh(NM*E3));
    elseif max(abs([E12,E3]))<MIN
        valeq=(NM^2/(2*E1^3))*(2*(E12+E1)*(expl(2,NM*E1)+expl(2,-NM*E1))...
            -2*NM*E12*E1*sinh(NM*E1));
    elseif max(abs([E1,E3]))<MIN
        valeq=(NM^2/(E12^3))*(4*(E1+E12+E3)*sinh(0.5*NM*E12)^2 ...
            -NM*(E1+E3)*E12*sinh(NM*E12));
    elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=2.*E1.^(-2).*NM.^2.*((-1)+cosh(E1.*NM));
    elseif (abs(E1-E12)<MIN && abs(E12-E3)>=MIN)
       valeq=exp(1).^((-1).*(E1+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
          E1.*NM)+(-1).*exp(1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1) ...
          .*(E1+(-1).*E3).^(-1).*E3.^(-1).*NM;
    elseif (abs(E1-E3)<MIN && abs(E1-E12)>=MIN)
       valeq=16.*E1.^(-2).*(E1+(-1).*E12).^(-2).*sinh((1/2).*E1.*NM).^2.*sinh(( ...
        1/2).*(E1+(-1).*E12).*NM).^2;
    elseif (abs(E12-E3)<MIN && abs(E1-E3)>=MIN)
       valeq=exp(1).^((-1).*(E1+E12).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
          E1.*NM)+(-1).*exp(1).^(E12.*NM)).*((-1)+exp(1).^(E12.*NM)).*E1.^( ...
          -1).*(E1+(-1).*E12).^(-1).*E12.^(-1).*NM;
    else
%        valeq=exp(1).^((-1).*(E1+E12+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1) ...
%           .^(E1.*NM)+(-1).*exp(1).^(E12.*NM)).*(exp(1).^(E12.*NM)+(-1).*exp( ...
%           1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1).*(E1+(-1).*E12) ...
%           .^(-1).*(E12+(-1).*E3).^(-1).*E3.^(-1);    
        valeq=2*( -coshl(4,NM*E1)+coshl(4,NM*E12)-coshl(4,NM*E3)-coshl(4,NM*(E1-E12))+...
           coshl(4,NM*(E1-E3))-coshl(4,NM*(E12-E3))+coshl(4,NM*(-E1+E12-E3)) )...
           *( (E12-E3)*E3*E1*(E12-E1) )^(-1);    
    end
    
    for ISet=0:7
        Ivec=de2bi(ISet,3);
        ZT1=ZE1(Ivec(1)+1);
        ZT12=ZE12(Ivec(2)+1);
        ZT3=ZE3(Ivec(3)+1); 
        
        valne=case8JPart(ZT1,ZT12,ZT3,N);
        
        for alphaSet=0:15
            aVec=de2bi(alphaSet,4);
            S4(aVec(1)+1,aVec(2)+1,aVec(3)+1,aVec(4)+1)=...
                  S4(aVec(1)+1,aVec(2)+1,aVec(3)+1,aVec(4)+1)+...
                  gamma4Prefac(Ivec,aVec,order,F(1),4)*...
                  valne*valeq;           
        end
    end
    
    
 %{
    % on different monomers
    valne=zeros(2,2,2);
    for I1=1:2
        for I2=1:2
            for I3=1:2
                ZT1=ZE1(I1);
                ZT12=ZE12(I2);
                ZT3=ZE3(I3);
                valne(I1,I2,I3)=case8JPart(ZT1,ZT12,ZT3,N);
                
            end
        end
    end

    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2

                    FV=[F(I1),F(I2),F(I3),F(I4)];
                    DFV=[DF(I1),DF(I2),DF(I3),DF(I4)];        
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      valeq*(FV(A1)*FV(A2)*FV(A3)*FV(A4)*valne(1,1,1)+...
                           FV(A1)*FV(A2)*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(1,1,2)+...
                           FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*FV(A4)*valne(1,2,1)+...
                           FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(1,2,2)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*FV(A4)*valne(2,1,1)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(2,1,2)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*FV(A4)*valne(2,2,1)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(2,2,2));                    
                end
            end
        end
    end
%}    
end