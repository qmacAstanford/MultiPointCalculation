clear;
close all

% polymer parameters
CHAIN=1;
N=10;
NM=10;
CHI=0;
chid=-0.1; % difference from spinodal chi

FA=0.6;
LAM=0.;
NFA=50;
FAV=linspace(0.1,0.5,NFA);
NQ=6;

% results
B = zeros(NFA,1);
C = zeros(NFA,1);

% wavevectors
if NQ==2
    QS1=[1,0,0];
    QS2=[-1,0,0];
    QS=[QS1',QS2'];
elseif NQ==3
    QS1=[1,0,0];
    QS2=(1/2)*[-1,+sqrt(3),0];
    QS3=(1/2)*[-1,-sqrt(3),0];
    QS4=-QS1;
    QS5=-QS2;
    QS6=-QS3;
    QS=[QS1',QS2',QS3',QS4',QS5',QS6'];
    NQ=6;
elseif NQ==6
    QS1=power(2,-1/2)*[1,1,0];
    QS2=power(2,-1/2)*[0,1,1];
    QS3=power(2,-1/2)*[1,0,1];
    QS4=power(2,-1/2)*[-1,1,0];
    QS5=power(2,-1/2)*[0,1,-1];
    QS6=power(2,-1/2)*[1,0,-1];
    QS7=-QS1;
    QS8=-QS2;
    QS9=-QS3;
    QS10=-QS4;
    QS11=-QS5;
    QS12=-QS6;
    QS=[QS1',QS2',QS3',QS4',QS5',QS6',...
        QS7',QS8',QS9',QS10',QS11',QS12'];
    NQ=12;
end

for IFA=1:NFA
FA=FAV(IFA)
    
% amplitude
AMP=linspace(0,50,100)/sqrt(NQ);

% find critical k
G=@(k) gamma2(k,CHAIN,N,NM,LAM,FA,CHI);
ks = fminbnd(G,1e-3,1e3);
% figure;hold
% fplot(G,[1e-3,1e3])
% plot([ks,ks],[0,2],'r--')

% find spinodal chi
chis = -gamma2(ks,CHAIN,N,NM,LAM,FA,CHI);

% calculate free energy
% coefficients
gam2 = chid*1;

gam3 = 0;
MIN=1e-10;
for J1=1:NQ
    for J2=1:NQ
        for J3=1:NQ
            Q1=QS(:,J1);
            Q2=QS(:,J2);
            Q3=QS(:,J3);
            if sum(power(Q1+Q2+Q3,2)) <= MIN
                gam3 = gam3+gamma3(ks,CHAIN,N,NM,LAM,FA);
            end
        end
    end
end

gam4 = 0;
MIN=1e-10;
for J1=1:NQ
    for J2=1:NQ
        for J3=1:NQ
            for J4=1:NQ
                Q1=QS(:,J1);
                Q2=QS(:,J2);
                Q3=QS(:,J3);
                Q4=QS(:,J4);
                if sum(power(Q1+Q2+Q3+Q4,2)) <= MIN
                    gam4 = gam4+gamma4(ks,Q1,Q2,Q3,Q4,...
                                        1,1,1,1,...
                                        CHAIN,N,NM,LAM,FA);
                end
            end
        end
    end
end

B(IFA) = -(1/6)*gam3;
C(IFA) = -(1/24)*gam4;

end

% free energy
IFA=NFA;
A=chid;
bF=A*power(AMP,2)+B(IFA)*power(AMP,3)+C(IFA)*power(AMP,4);
figure;plot(AMP,bF)

% k=logspace(-1,4,50);
% 
% Q1=[1,0,0];
% Q2=[0,-1,0];
% Q3=[-1,0,0];
% Q4=[0,1,0];
% 
% figure(1);hold
% figure(2);hold
% figure(3);hold
% 
% for ii=1:NLAM
%     col = (ii-1)./(NLAM-1);
%     LAM=LAMV(ii);
%     
%     G=gamma2(k,CHAIN,N,NM,LAM,FA,CHI);
%     figure(1);plot(k,1./G,'color',[col 0 1-col])
% 
%     if FA~=0.5
%         G=gamma3(k,CHAIN,N,NM,LAM,FA);
%         figure(2);plot(k,G,'color',[col 0 1-col])
%     end
% 
%     G=gamma4(k,Q1,Q2,Q3,Q4,...
%                     1,1,1,1,...
%                     CHAIN,N,NM,LAM,FA);
%     figure(3);plot(k,G,'color',[col 0 1-col])
% end
% 
% figure(1);set(gca,'xscale','log')
% figure(2);set(gca,'xscale','log')
% figure(3);set(gca,'xscale','log')