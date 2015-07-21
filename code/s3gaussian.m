function S3=s3gaussian(N,NM,LAM,FA,Q1,Q2,Q3,d)
%S3 is a three point correlation function
%For example:
%   S3(1,1,1)=SAAA
%   S3(1,1,2)=SAAB

S3=zeros(2,2,2);
MIN=5e-6;
%digits(5);

if sum(power(Q1+Q2+Q3,2)) <= MIN
    
    % Evaluate the quantities for s3 calculation
    FB=1-FA;
    F=[FA,FB];
    DF=[1,-1];
    
    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    Q3MAG=sqrt(sum(power(Q3,2)));
    
    R1=-Q1MAG*Q1MAG/(2*d);
    R2=-Q2MAG*Q2MAG/(2*d);
    R3=-Q3MAG*Q3MAG/(2*d);
    
    Z1=exp(R1*NM);
    Z2=exp(R2*NM);
    Z3=exp(R3*NM);
    
    Z1L=exp(R1*NM)*LAM;
    Z2L=exp(R2*NM)*LAM;
    Z3L=exp(R3*NM)*LAM;

    % Case 1: J1=J2=J3
        S3=case1(S3,N,NM,R1,R2,Z1,Z2,F,MIN);
        S3=case1(S3,N,NM,R1,R3,Z1,Z3,F,MIN);
        S3=case1(S3,N,NM,R2,R3,Z2,Z3,F,MIN);

    % Case 2: J1<J2=J3
        % J1<J2=J3
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(2,1,1)=1;
        SDEL(1,2,2)=1;
        S3=case2(S3,N,NM,R1,R2,Z1,Z2,Z1L,F,DF,SDEL,1,2,MIN);
        S3=case2(S3,N,NM,R1,R3,Z1,Z3,Z1L,F,DF,SDEL,1,2,MIN);

        % J2<J1=J3
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(2,1,2)=1;
        SDEL(1,2,1)=1;
        S3=case2(S3,N,NM,R2,R3,Z2,Z3,Z2L,F,DF,SDEL,2,1,MIN);
        S3=case2(S3,N,NM,R2,R1,Z2,Z1,Z2L,F,DF,SDEL,2,1,MIN);

        % J3<J1=J2
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(1,1,2)=1;
        SDEL(2,2,1)=1;
        S3=case2(S3,N,NM,R3,R1,Z3,Z1,Z3L,F,DF,SDEL,3,1,MIN);
        S3=case2(S3,N,NM,R3,R2,Z3,Z2,Z3L,F,DF,SDEL,3,1,MIN);

    % Case 3: J1=J2<J3
        % J1=J2<J3
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(1,1,2)=1;
        SDEL(2,2,1)=1;
        S3=case3(S3,N,NM,R1,R3,Z1,Z3,Z3L,F,DF,SDEL,2,3,MIN);
        S3=case3(S3,N,NM,R2,R3,Z2,Z3,Z3L,F,DF,SDEL,2,3,MIN);

        % J1=J3<J2
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(1,2,1)=1;
        SDEL(2,1,2)=1;
        S3=case3(S3,N,NM,R1,R2,Z1,Z2,Z2L,F,DF,SDEL,3,2,MIN);
        S3=case3(S3,N,NM,R3,R2,Z3,Z2,Z2L,F,DF,SDEL,3,2,MIN);

        % J2=J3<J1
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(1,2,2)=1;
        SDEL(2,1,1)=1;
        S3=case3(S3,N,NM,R2,R1,Z2,Z1,Z1L,F,DF,SDEL,3,1,MIN);
        S3=case3(S3,N,NM,R3,R1,Z3,Z1,Z1L,F,DF,SDEL,3,1,MIN);

    % Case 4: J1<J2<J3
        SDEL=ones(2,2,2);
        % J1<J2<J3
        S3=case4(S3,N,NM,R1,R3,Z1,Z3,Z1L,Z3L,F,DF,SDEL,1,2,3,MIN);
        
        %J3<J2<J1
        S3=case4(S3,N,NM,R3,R1,Z3,Z1,Z3L,Z1L,F,DF,SDEL,3,2,1,MIN);
        
        %J1<J3<J2
        S3=case4(S3,N,NM,R1,R2,Z1,Z2,Z1L,Z2L,F,DF,SDEL,1,3,2,MIN);
        
        %J2<J3<J1
        S3=case4(S3,N,NM,R2,R1,Z2,Z1,Z2L,Z1L,F,DF,SDEL,2,3,1,MIN);
        
        %J2<J1<J3
        S3=case4(S3,N,NM,R2,R3,Z2,Z3,Z2L,Z3L,F,DF,SDEL,2,1,3,MIN);
        
        %J3<J1<J2
        S3=case4(S3,N,NM,R3,R2,Z3,Z2,Z3L,Z2L,F,DF,SDEL,3,1,2,MIN);

    S3=real(S3);
end
end

function S3=case1(S3,N,NM,R1,R2,Z1,Z2,F,MIN)
% Case 1: J1=J2=J3
    % on same monomer
    if abs(R1-R2)<MIN
        if abs(R1)<MIN
            valeq=power(NM,3)/6;
        else
            valeq=R1.^(-3).*(2+NM.*R1+((-2)+NM.*R1).*Z1);
        end
    else
        valeq=R1.^(-2).*(R1+(-1).*R2).^(-1).*R2.^(-2).*((-1).*NM.*R1.*R2.^2+ ...
            R2.^2.*((-1)+Z1)+R1.^2.*(1+NM.*R2+(-1).*Z2));
    end
%     valeq=finite(valeq);
    
    S3(1,1,1)=S3(1,1,1)+2*F(1)*N*valeq;
    S3(2,2,2)=S3(2,2,2)+2*F(2)*N*valeq;
end

function S3=case2(S3,N,NM,R1,R2,Z1,Z2,Z1L,F,DF,SDEL,A1,A2,MIN)
% Case 2: J1<J2=J3
    % on same monomer
    if abs(R1-R2)<MIN
        if abs(R1)<MIN
            valeq=power(NM,3)/2;
        else
            valeq=R1.^(-3).*((-1)+Z1).*Z1.^(-1).*(1+((-1)+NM.*R1).*Z1);
        end
    else
        valeq=R1.^(-2).*(R1+(-1).*R2).^(-1).*R2.^(-1).*((-1)+Z1).*Z1.^(-1).*(R1+ ...
            R2.*((-1)+Z1)+(-1).*R1.*Z2);
    end
    
    % on different monomers
    if abs(Z1-1)<MIN
        valne1=(1/2).*((-1).*N+N.^2);
    else
        valne1=(-1).*((-1)+Z1).^(-2).*Z1.*(1+N.*((-1)+Z1)+(-1).*Z1.^N);
    end
    
    if abs(Z1L-1)<MIN
        valne2=(1/2).*((-1).*N+N.^2);
    else
        valne2=(-1).*((-1)+Z1L).^(-2).*Z1L.*(1+N.*((-1)+Z1L)+(-1).*Z1L.^N);
    end
    valeq=finite(valeq);
    valne1=finite(valne1);
    valne2=finite(valne2);
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                FV=[F(I1),F(I2),F(I3)];
                DFV=[DF(I1),DF(I2),DF(I3)];
                S3(I1,I2,I3)=S3(I1,I2,I3)+SDEL(I1,I2,I3)*...
                           valeq*(FV(A1)*FV(A2)*valne1+...
                             FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*valne2);
            end
        end
    end
end

function S3=case3(S3,N,NM,R1,R2,Z1,Z2,Z2L,F,DF,SDEL,A2,A3,MIN)
% Case 3: J1=J2<J3
    % on same monomer
    if abs(R1-R2)<MIN
        if abs(R1)<MIN
            valeq=power(NM,3)/2;
        else
            valeq=R1.^(-3).*((-1)+Z1).*Z1.^(-1).*(1+((-1)+NM.*R1).*Z1);
        end
    else
        valeq=R1.^(-1).*R2.^(-2).*((-1).*R1+R2).^(-1).*(R2+(-1).*R2.*Z1+R1.*(( ...
            -1)+Z2)).*((-1)+Z2).*Z2.^(-1);
    end
    
    % on different monomers
    if abs(Z2-1)<MIN
        valne1=(1/2).*((-1).*N+N.^2);
    else
        valne1=(-1).*((-1)+Z2).^(-2).*Z2.*(1+N.*((-1)+Z2)+(-1).*Z2.^N);
    end
    if abs(Z2L-1)<MIN
        valne2=(1/2).*((-1).*N+N.^2);
    else
        valne2=(-1).*((-1)+Z2L).^(-2).*Z2L.*(1+N.*((-1)+Z2L)+(-1).*Z2L.^N);
    end
    valeq=finite(valeq);
    valne1=finite(valne1);
    valne2=finite(valne2);
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                FV=[F(I1),F(I2),F(I3)];
                DFV=[DF(I1),DF(I2),DF(I3)];
                S3(I1,I2,I3)=S3(I1,I2,I3)+SDEL(I1,I2,I3)*...
                           valeq*(FV(A2)*FV(A3)*valne1+...
                             FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*valne2);
            end
        end
    end
end

function S3=case4(S3,N,NM,R1,R2,Z1,Z2,Z1L,Z2L,F,DF,SDEL,A1,A2,A3,MIN)
% Case 4: J1<J2<J3
    ZE1=[Z1,Z1L];
	ZE2=[Z2,Z2L];

    % on same monomer
    if abs(R1-R2)<MIN
        if abs(R1)<MIN
            valeq=power(NM,3);
        else
            valeq=2.*NM.*R1.^(-2).*((-1)+(1/2).*(Z1.^(-1)+Z1));
            
        end
    else
%         valeq=R1.^(-1).*(R1+(-1).*R2).^(-1).*R2.^(-1).*((-1)+Z1).*Z1.^(-1).*(Z1+ ...
%           (-1).*Z2).*((-1)+Z2).*Z2.^(-1);
          valeq=(R1.*R2.*(R2-R1)).^(-1).*...
                (-1).*(expl(1,-R1*NM).*(expl(2,R2*NM)+expl(2,-R2.*NM))-...
                expl(1,-R2*NM).*(expl(2,R1*NM)+expl(2,-R1.*NM)));
    end
%    valeq(valeq>1e5)=0;
    
    % on different monomers
%     MIN=abs(exp(NM)*exp(R1)*(1-exp(MIN)));
    valne=zeros(2,2);
    %valne(1,1) = no lambda term
    %valne(1,2) = lambda^(J3-J2) term
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT2=ZE2(I2);
            if (abs(ZT1-1)<MIN && abs(ZT2-1)>=MIN)
                valne(I1,I2)=(-1/2).*((-1)+ZT2).^(-3).*ZT2.*(2+N.^2.*((-1)+ZT2).^2+(-2).* ...
                ZT2.^N+(-1).*N.*(3+(-4).*ZT2+ZT2.^2));
            elseif (abs(ZT1-1)>=MIN && abs(ZT2-1)<MIN)
                valne(I1,I2)=(1/2).*(2+(-3).*N+N.^2).*((-1)+ZT1).^(-2).*ZT1+...
                    (-1/2).*((-1)+ZT1).^(-3).*ZT1.*(2.*ZT1+N.*ZT1+(-1).*N.^2.*ZT1+(-1) ...
                    .*N.*ZT1.^2+N.^2.*ZT1.^2+(-2).*ZT1.^N);
            elseif (abs(ZT1-1)<MIN && abs(ZT2-1)<MIN)
                valne(I2,I2)=(1/6).*N.*(2+(-3).*N+N.^2);
            elseif abs(ZT1-ZT2)<MIN
                valne(I1,I2)=((-1)+ZT1).^(-3).*ZT1.*(2.*ZT1+(-1).*N.*ZT1+N.*ZT1.^2+(-1).*N.* ...
                ZT1.^N+(-2).*ZT1.^(1+N)+N.*ZT1.^(1+N));
            else
                if (ZT1>0.2 && ZT2>0.2)
                    valne(I1,I2)=((ZT1-1)^2*(ZT1-ZT2)*(ZT2-1)^2)^(-1)* ...
                                  ZT2*ZT1*...
                                  s3case4JnearOne(N,log(ZT1),log(ZT2));
                else
                    valne(I1,I2)=((-1)+ZT1).^(-2).*ZT1.*(ZT1+(-1).*ZT2).^(-1).*((-1)+ZT2).^(-2).* ...
                     ZT2.*(ZT1.^N.*((-1)+ZT2).^2+(-1).*((-2)+N).*ZT2+((-1)+N).*ZT2.^2+( ...
                      -1).*ZT2.^N+(-1).*ZT1.^2.*((-1)+N+(-1).*N.*ZT2+ZT2.^N)+ZT1.*((-2)+ ...
                      N+(-1).*N.*ZT2.^2+2.*ZT2.^N));
                end
            end
        
            if (valne<0)
                disp('Error :: Numerial issue with negative result from summation')
            end
        end
    end
    %valeq=finite(valeq);
    %valne=finite(valne);

    for I1=1:2
        for I2=1:2
            for I3=1:2
                FV=[F(I1),F(I2),F(I3)];
                DFV=[DF(I1),DF(I2),DF(I3)];
                S3(I1,I2,I3)=S3(I1,I2,I3)+SDEL(I1,I2,I3)*...
                   valeq*(FV(A1)*FV(A2)*FV(A3)*valne(1,1)+...
                      FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*valne(1,2)+...
                      FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*valne(2,1)+...
                      FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*valne(2,2));
            end
        end
    end
end
function out=s3case4JnearOne(N,x1,x2)
out=0;
out=out+expl(1,N*x1)*(expl(1,x2))^2 + expl(2,x2)*(expl(2,x2)+2*x2);
out=out-(N-2)*expl(3,x2) + (N-1)*expl(3,2*x2) - expl(3,N*x2);
out=out-(2*expl(1,x1)+expl(1,x1)^2)*...
        (-N*expl(3,x2)+expl(3,N*x2)+0.5*N*(N-1)*x2^2)+...
        N*expl(3,x2)-expl(3,N*x2);
out=out+expl(1,x1)*(-2*N*(expl(3,x2)+0.5*x2^2)-N*expl(1,x2)^2+2*expl(2,N*x2))...
       -2*N*expl(3,x2)-N*(2*x2*expl(2,x2)+expl(2,x2)^2)+2*expl(3,N*x2);
end
function xf=finite(x)
    MIN=-1e10;
    MAX=1e10;
    xf=x;
    xf(xf>MAX | xf<MIN)=0;
    xf(isnan(xf))=0;
end