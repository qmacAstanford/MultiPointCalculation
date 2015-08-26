function S4=case8(N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,order)
    if N<4
        S4=zeros(2,2,2,2);
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
    
    S4=BinomialSum(valeq,valne,order,F);
    
end