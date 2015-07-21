function temp()

a=0.99999;
b=0.867362;
c=0.897869;
N=10;

out=case8JPart(a,b,c,N)

% a=0.83636;
% cs=logspace(-4,2,100);
% N=10;
% data=zeros(100,1);
% for j=1:100
%     c=cs(j);
%     b=c;
%     data(j)=case8JPart(a,b,c,N);
% end
% 
% loglog(cs,data,'+')
% xlabel('c')
% ylabel('out')
% 

end
function out=case8JPart(a,b,c,N)
global Flage
    x1=log(a);
    x2=log(b);
    x3=log(c);
    if( min([a,b,c])<0)
        error('inputs must be positive')
    end
    if a*b*c==0
        out=0;
        return 
    end
    
    
    if 0 %max([a,b,c])<0.5 %I included this for speedup but in rarly gets used 
        ZT1=a;
        ZT12=b;
        ZT3=c;
        out=((-1)+ZT1).^(-2).*(ZT1+(-1).*ZT12).^(-1).*((-1)+ZT12).^(-2).*(ZT1+ ...
        (-1).*ZT3).^(-1).*((-1)+ZT3).^(-2).*((-1).*ZT12+ZT3).^(-1).*(3.* ...
        ZT1.^3.*ZT12.^2.*ZT3+(-1).*N.*ZT1.^3.*ZT12.^2.*ZT3+(-2).*ZT1.^4.* ...
        ZT12.^2.*ZT3+N.*ZT1.^4.*ZT12.^2.*ZT3+(-1).*ZT1.^(1+N).*ZT12.^2.* ...
        ZT3+(-3).*ZT1.^2.*ZT12.^3.*ZT3+N.*ZT1.^2.*ZT12.^3.*ZT3+ZT1.^4.* ...
        ZT12.^3.*ZT3+(-1).*N.*ZT1.^4.*ZT12.^3.*ZT3+2.*ZT1.^(1+N).* ...
        ZT12.^3.*ZT3+2.*ZT1.^2.*ZT12.^4.*ZT3+(-1).*N.*ZT1.^2.*ZT12.^4.* ...
        ZT3+(-1).*ZT1.^3.*ZT12.^4.*ZT3+N.*ZT1.^3.*ZT12.^4.*ZT3+(-1).* ...
        ZT1.^(1+N).*ZT12.^4.*ZT3+ZT1.^2.*ZT12.^(1+N).*ZT3+(-2).*ZT1.^3.* ...
        ZT12.^(1+N).*ZT3+ZT1.^4.*ZT12.^(1+N).*ZT3+(-3).*ZT1.^3.*ZT12.* ...
        ZT3.^2+N.*ZT1.^3.*ZT12.*ZT3.^2+2.*ZT1.^4.*ZT12.*ZT3.^2+(-1).*N.* ...
        ZT1.^4.*ZT12.*ZT3.^2+ZT1.^(1+N).*ZT12.*ZT3.^2+3.*ZT1.*ZT12.^3.* ...
        ZT3.^2+(-1).*N.*ZT1.*ZT12.^3.*ZT3.^2+N.*ZT1.^4.*ZT12.^3.*ZT3.^2+( ...
        -3).*ZT1.^(1+N).*ZT12.^3.*ZT3.^2+(-2).*ZT1.*ZT12.^4.*ZT3.^2+N.* ...
        ZT1.*ZT12.^4.*ZT3.^2+(-1).*N.*ZT1.^3.*ZT12.^4.*ZT3.^2+2.*ZT1.^(1+ ...
        N).*ZT12.^4.*ZT3.^2+(-1).*ZT1.*ZT12.^(1+N).*ZT3.^2+3.*ZT1.^3.* ...
        ZT12.^(1+N).*ZT3.^2+(-2).*ZT1.^4.*ZT12.^(1+N).*ZT3.^2+3.*ZT1.^2.* ...
        ZT12.*ZT3.^3+(-1).*N.*ZT1.^2.*ZT12.*ZT3.^3+(-1).*ZT1.^4.*ZT12.* ...
        ZT3.^3+N.*ZT1.^4.*ZT12.*ZT3.^3+(-2).*ZT1.^(1+N).*ZT12.*ZT3.^3+(-3) ...
        .*ZT1.*ZT12.^2.*ZT3.^3+N.*ZT1.*ZT12.^2.*ZT3.^3+(-1).*N.*ZT1.^4.* ...
        ZT12.^2.*ZT3.^3+3.*ZT1.^(1+N).*ZT12.^2.*ZT3.^3+ZT1.*ZT12.^4.* ...
        ZT3.^3+(-1).*N.*ZT1.*ZT12.^4.*ZT3.^3+N.*ZT1.^2.*ZT12.^4.*ZT3.^3+( ...
        -1).*ZT1.^(1+N).*ZT12.^4.*ZT3.^3+2.*ZT1.*ZT12.^(1+N).*ZT3.^3+(-3) ...
        .*ZT1.^2.*ZT12.^(1+N).*ZT3.^3+ZT1.^4.*ZT12.^(1+N).*ZT3.^3+(-2).* ...
        ZT1.^2.*ZT12.*ZT3.^4+N.*ZT1.^2.*ZT12.*ZT3.^4+ZT1.^3.*ZT12.*ZT3.^4+ ...
        (-1).*N.*ZT1.^3.*ZT12.*ZT3.^4+ZT1.^(1+N).*ZT12.*ZT3.^4+2.*ZT1.* ...
        ZT12.^2.*ZT3.^4+(-1).*N.*ZT1.*ZT12.^2.*ZT3.^4+N.*ZT1.^3.*ZT12.^2.* ...
        ZT3.^4+(-2).*ZT1.^(1+N).*ZT12.^2.*ZT3.^4+(-1).*ZT1.*ZT12.^3.* ...
        ZT3.^4+N.*ZT1.*ZT12.^3.*ZT3.^4+(-1).*N.*ZT1.^2.*ZT12.^3.*ZT3.^4+ ...
        ZT1.^(1+N).*ZT12.^3.*ZT3.^4+(-1).*ZT1.*ZT12.^(1+N).*ZT3.^4+2.* ...
        ZT1.^2.*ZT12.^(1+N).*ZT3.^4+(-1).*ZT1.^3.*ZT12.^(1+N).*ZT3.^4+(-1) ...
        .*ZT1.^2.*ZT12.*ZT3.^(1+N)+2.*ZT1.^3.*ZT12.*ZT3.^(1+N)+(-1).* ...
        ZT1.^4.*ZT12.*ZT3.^(1+N)+ZT1.*ZT12.^2.*ZT3.^(1+N)+(-3).*ZT1.^3.* ...
        ZT12.^2.*ZT3.^(1+N)+2.*ZT1.^4.*ZT12.^2.*ZT3.^(1+N)+(-2).*ZT1.* ...
        ZT12.^3.*ZT3.^(1+N)+3.*ZT1.^2.*ZT12.^3.*ZT3.^(1+N)+(-1).*ZT1.^4.* ...
        ZT12.^3.*ZT3.^(1+N)+ZT1.*ZT12.^4.*ZT3.^(1+N)+(-2).*ZT1.^2.* ...
        ZT12.^4.*ZT3.^(1+N)+ZT1.^3.*ZT12.^4.*ZT3.^(1+N));
    elseif ((1-a)<0.1/N) && ((1-b)<0.1/N) && ((1-c)<0.1/N && Flage) 
        out=-N*(N-1)*(N-2)*(N-3)*((N+1)*(3-a-b-c)-5)/120;
    else
        [a,b,c]=offset(a,b,c,3.0e-6); %keeps a and b and c different
        
        outs=[0;0;0;0];
        v1=[x1,x1,x3,x3,x1-x3,3*x2-N*x2];
        v2=[expl(2,x1),expl(2,x1),expl(2,x3),expl(2,x3),...
            expl(2,x1)-expl(2,x3),expl(2,3*x2)-expl(2,N*x2)];
        outs(1)= hooperHumperdink(v1,v2);

        v1=[x2,x2,x3,x3,x2-x3,3*x1-N*x1];
        v2=[expl(2,x2),expl(2,x2),expl(2,x3),expl(2,x3),...
            expl(2,x2)-expl(2,x3),expl(2,3*x1)-expl(2,N*x1)];
        outs(2)=(-1)*hooperHumperdink(v1,v2);

        v1=[x1,x1,x2,x2,x1-x2,3*x3-N*x3];
        v2=[expl(2,x1),expl(2,x1),expl(2,x2),expl(2,x2),...
            expl(2,x1)-expl(2,x2),expl(2,3*x3)-expl(2,N*x3)];
        outs(3)=(-1)*hooperHumperdink(v1,v2);

        v1=[x1,x2,x3,x1-x2,x1-x3,x2-x3];
        v2=[expl(2,x1),expl(2,x2),expl(2,x3),...
            expl(2,x1)-expl(2,x2),...
            expl(2,x1)-expl(2,x3),...
            expl(2,x2)-expl(2,x3)];
        outs(4)=hooperHumperdink(v1,v2)*(3-N);
        out=sum(outs);
        
        if 0 % abs(max(outs)/out)>1e6
            sprintf('a=%g, b=%g, c=%g',a,b,c)
            disp('outs')
            disp(outs)
            disp('out')
            disp(out)
            error('disaster')
        end
        out=out*a*b*c*( (a-1)^2*(b-1)^2*(c-1)^2*(a-b)*(a-c)*(b-c) )^(-1);
    end 
end