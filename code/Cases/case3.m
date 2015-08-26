function S4=case3(N,NM,R1,R12,R3,Z12,Z12L,F,order)
% Case 3: J1==J2 <J3==J4
    if N<2
        S4=zeros(2,2,2,2);
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE12=Z12;
    ZE12L=Z12L;

    % on same monomer
    offset=10^-12;
    if abs(E1)<offset
        E1=offset;
    end
    if abs(E12)<offset
        E12=offset;
    end
    if abs(E3)<offset
        E3=offset;
    end
    MIN=(10^-4)/NM;
    flag=0;
    if max(abs([E1,E12,E3]))<MIN
        valeq=NM^4*(NM*(E1-E12+E3)+3)/12;
    elseif max(abs([E1,E12]))<MIN
        valeq=(-NM^2/(12*E3^3))*((NM^2*(2*E1-E12)+6*NM)*E3^2+...
                    2*NM*E3*((1+2*exp(NM*E3))*E12+(-expl(1,NM*E3))*E1)+...
                    6*(-expl(1,NM*E3))*(E3+E12));
    elseif max(abs([E12,E3]))<MIN
        valeq=(-NM^2/(12*E1^3))*((NM^2*(2*E3-E12)+6*NM)*E1^2+...
                    2*NM*E1*((1+2*exp(NM*E1))*E12+(-expl(1,NM*E1))*E3)+...
                    6*(-expl(1,NM*E1))*(E1+E12));
    elseif max(abs([E1,E3]))<MIN
        valeq=((2*NM^2+NM^3*(E1+E3))*exp(-NM*E12)*E12^3+...
                2*NM^2*(E1+E3)*E12^2+...
                (expl(1,-NM*E12))*(...
                      (3*NM^2*(E1+E3)+4*NM)*E12^2+...
                      4*NM*(E1+E3)*E12 ...
                      -2*(E1+E12+E3)*expl(1,NM*E12)))/(2*E12^5);
    elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(exp(1).^E1).^((-1).*NM).*log(exp(1).^E1).^(-4).*(1+(exp(1).^E1) ...
       .^NM.*((-1)+NM.*log(exp(1).^E1))).^2;
    elseif abs(E1-E12)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*log(exp(1).^E1).^(-2).*(1+(exp(1).^E1) ...
       .^NM.*((-1)+NM.*log(exp(1).^E1))).*(((-1)+(exp(1).^E1).^NM).*log( ...
       exp(1).^E1).^(-1)+(1+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E3).^( ...
       -1)).*(log(exp(1).^E1)+(-1).*log(exp(1).^E3)).^(-1);
    elseif abs(E1-E3)<MIN
       flag=1;
       valeq=(exp(1).^E12).^((-1).*NM).*log(exp(1).^E1).^(-2).*(log(exp(1).^E1) ...
       +(-1).*log(exp(1).^E12)).^(-2).*log(exp(1).^E12).^(-2).*(((-1)+( ...
       exp(1).^E12).^NM).*log(exp(1).^E1)+(-1).*((-1)+(exp(1).^E1).^NM).* ...
       log(exp(1).^E12)).^2;
    elseif abs(E12-E3)<MIN
       valeq=log(exp(1).^E1).^(-1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^( ...
       -1).*log(exp(1).^E12).^(-3).*((-1).*((-1)+(exp(1).^E12).^NM).*log( ...
       exp(1).^E1)+((-1)+(exp(1).^E1).^NM).*log(exp(1).^E12)).*((-1)+( ...
       exp(1).^E12).^((-1).*NM)+NM.*log(exp(1).^E12));
    else
       valeq=(-1).*(exp(1).^E12).^((-1).*NM).*log(exp(1).^E1).^(-1).*(log(exp( ...
       1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1).^E12).^(-1).*((( ...
       -1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+(-1).*((-1)+(exp(1).^E1) ...
       .^NM).*log(exp(1).^E12)).*(((-1)+(exp(1).^E12).^NM).*log(exp(1) ...
       .^E12).^(-1)+(1+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E3).^(-1)).*( ...
       log(exp(1).^E12)+(-1).*log(exp(1).^E3)).^(-1);
    end
    if valeq<0 && isreal(E1) && isreal(E12) && isreal(E3)
        sprintf('E1=%g, E12=%g,E3=%g, valeq=%g, flag=%d, NM=%d',E1,E12,E3,valeq,flag,NM)
        disp(E1)
        disp(E12)
        error('value cannot be negitive')
    end
    
    % on different monomers
    valne1=twoSum(ZE12,N);
    valne2=twoSum(ZE12L,N);
    
    valne=zeros(2,2,2);
    valne(:,1,:)=ones(2,2)*valne1;
    valne(:,2,:)=ones(2,2)*valne2;
    
    S4=BinomialSum(valeq,valne,order,F);
end