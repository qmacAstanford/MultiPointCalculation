clc


npts=50;
x=logspace(-9,0.1,npts);

old=zeros(npts,1);
new=old;
zabc=old;
zac=old;
zab=old;
zbc=old;
abc=old;
ab=old;
ac=old;
bc=old;
chosen=old;

for pt=1:npts

NM=10;
E1=-x(pt);
E12=0;
E3=-x(pt);
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

    %if max(abs([E1,E12,E3]))<MIN
        zabc(pt)=NM^4*(NM*(E1+E12+E3)+5)/120;
        
    %elseif max(abs([E1,E12]))<MIN
        zab(pt)=chicken(E1,E12,E3,NM);
        
    %elseif max(abs([E12,E3]))<MIN
        zbc(pt)=chicken(E12,E3,E1,NM);
        
    %elseif max(abs([E1,E3]))<MIN
        zac(pt)=chicken(E1,E3,E12,NM);      
    %elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
        abc(pt)=(1/2).*E1.^(-4).*((-2).*(3+E1.*NM)+exp(1).^(E1.*NM).*(6+E1.*NM.*(( ...
         -4)+E1.*NM)));
    %elseif (abs(E1-E12)<MIN )%&& abs(E12-E3)>=MIN) 
        ab(pt)=E1.^(-3).*(E1+(-1).*E3).^(-2).*E3.^(-2).*(2.*((-1)+exp(1).^(E1.* ...
          NM)).*E3.^3+(2+exp(1).^(E1.*NM)).*E1.^2.*E3.^2.*NM+E1.^3.*((-1)+ ...
          exp(1).^(E3.*NM)+(-1).*E3.*NM)+(-1).*E1.*E3.^2.*((-3)+E3.*NM+exp( ...
          1).^(E1.*NM).*(3+E3.*NM)));
    %elseif (abs(E1-E3)<MIN )%&& abs(E1-E12)>=MIN)
%         ac(pt)=E1.^(-3).*(E1+(-1).*E12).^(-2).*E12.^(-2).*(2.*((-1)+exp(1).^(E1.* ...
%           NM)).*E12.^3+(2+exp(1).^(E1.*NM)).*E1.^2.*E12.^2.*NM+E1.^3.*((-1)+ ...
%           exp(1).^(E12.*NM)+(-1).*E12.*NM)+(-1).*E1.*E12.^2.*((-3)+E12.*NM+ ...
%           exp(1).^(E1.*NM).*(3+E12.*NM)));
        ac(pt)=((-3*E1*expl(4,NM*E1)+NM*E1^2*expl(3,NM*E1))*E12^2+...
          (2*expl(4,NM*E1)-NM*E1*expl(3,NM*E1))*E12^3+...
          E1^3*expl(4,NM*E12))/(E1^3*E12^2*(E1-E12)^2);
    %elseif (abs(E12-E3)<MIN )%&& abs(E1-E3)>=MIN)
       bc(pt)=E1.^(-2).*(E1+(-1).*E12).^(-2).*E12.^(-3).*(((-1)+exp(1).^(E1.*NM) ...
          ).*E12.^3+(-1).*E1.*E12.^3.*NM+E1.^2.*E12.*(3+2.*E12.*NM+exp(1).^( ...
          E12.*NM).*((-3)+E12.*NM))+(-1).*E1.^3.*(2+E12.*NM+exp(1).^(E12.* ...
          NM).*((-2)+E12.*NM)));
    %else
       old(pt)=E1.^(-2).*(E1+(-1).*E12).^(-1).*E12.^(-2).*(E1+(-1).*E3).^(-1).*( ...
          E12+(-1).*E3).^(-1).*E3.^(-2).*(((-1)+exp(1).^(E1.*NM)).*E12.^2.*( ...
          E12+(-1).*E3).*E3.^2+E1.*E12.^2.*E3.^2.*((-1).*E12+E3).*NM+E1.^3.* ...
          ((-1).*((-1)+exp(1).^(E12.*NM)).*E3.^2+E12.*E3.^2.*NM+E12.^2.*(( ...
          -1)+exp(1).^(E3.*NM)+(-1).*E3.*NM))+E1.^2.*(((-1)+exp(1).^(E12.* ...
          NM)).*E3.^3+(-1).*E12.*E3.^3.*NM+E12.^3.*(1+(-1).*exp(1).^(E3.*NM) ...
          +E3.*NM)));
    %end
  %}  
    
    
    MIN=10^-5;
    if max(abs([E1,E12,E3]))<MIN
        valeq=NM^4*(2*NM*E1-NM*E12+6)/12;
    elseif max(abs([E1,E12]))<MIN
        valeq=NM^2*exp(NM*E3)*expl(1,-NM*E3)*(3*(E12+E3)*expl(2,-NM*E3)-3*NM*E3^2 ...
              +NM*(E1+E12)*E3*expl(1,-NM*E3))/(6*E3^3);
    elseif max(abs([E12,E3]))<MIN
        valeq=(NM^2/(E1^3))*((E1+E12)*expl(2,NM*E1)-0.5*NM*E1*E12*expl(1,NM*E1));
    elseif max(abs([E1,E3]))
        valeq=(NM/(2*E12^4))*(2*(E1+E12+E3)*...
               (expl(3,NM*E12)+expl(3,-NM*E12)+NM*E12*expl(2,-NM*E12))+...
               NM^2*(E1+E3)*E12^2*expl(1,-NM*E12)+...
               NM*E12*E3*(expl(2,-NM*E12)-expl(2,NM*E12)));
    elseif(abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*NM.*log(exp(1) ...
       .^E1).^(-3).*(1+(exp(1).^E1).^NM.*((-1)+NM.*log(exp(1).^E1)));
    elseif abs(E1-E12)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*(1+(-1).*(exp(1).^E3).^((-1).*NM)).*(( ...
       exp(1).^E1).^NM+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E1).^(-2).*( ...
       1+(exp(1).^E1).^NM.*((-1)+NM.*log(exp(1).^E1))).*(log(exp(1).^E1)+ ...
       (-1).*log(exp(1).^E3)).^(-1).*log(exp(1).^E3).^(-1);
    elseif abs(E1-E3)<MIN
       valeq=(exp(1).^E12).^((-1).*NM).*(1+(-1).*(exp(1).^E1).^((-1).*NM)).*(( ...
       exp(1).^E1).^NM+(-1).*(exp(1).^E12).^NM).*log(exp(1).^E1).^(-2).*( ...
       log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-2).*log(exp(1).^E12).^( ...
       -1).*((-1).*((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+((-1)+(exp( ...
       1).^E1).^NM).*log(exp(1).^E12));
    elseif abs(E12-E3)<MIN
       valeq=(1+(-1).*(exp(1).^E12).^((-1).*NM)).*NM.*log(exp(1).^E1).^(-1).*( ...
       log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1).^E12).^( ...
       -2).*((-1).*((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+((-1)+(exp( ...
       1).^E1).^NM).*log(exp(1).^E12));
    else
        valeq=expl(1,NM*E1)*(expl(1,-NM*E1)-expl(1,-NM*E12))*...
              ( expl(2,NM*E3)*E12 - expl(2,NM*E12)*E3 )/...
              (E1*E12*E3*(E12-E3)*(E1-E12));
%        valeq=(-1).*(exp(1).^E12).^((-1).*NM).*(1+(-1).*(exp(1).^E3).^((-1).*NM) ...
%        ).*((exp(1).^E12).^NM+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E1).^( ...
%        -1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1) ...
%        .^E12).^(-1).*(((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+(-1).*(( ...
%        -1)+(exp(1).^E1).^NM).*log(exp(1).^E12)).*(log(exp(1).^E12)+(-1).* ...
%        log(exp(1).^E3)).^(-1).*log(exp(1).^E3).^(-1);
    end  
   
    chosen(pt)=valeq;
    
end

p=semilogx(x,chosen,'-',x,old,'o',x,new,'+',x,zabc,'o',x,zac,'+',x,zab,'o',x,zbc,'.',x,abc,'o',x,ab,'.',x,ac,'.',x,bc,'.');
L=legend('chosen','old','new','zabc','zac','zab','zbc','abc','ab','ac','bc');
ylim([0,max(chosen)])
set(L,'location','eastoutside')
set(p(1),'LineWidth',3)







