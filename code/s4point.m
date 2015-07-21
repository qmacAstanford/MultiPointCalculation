% function [s4aaaa,Qmag,Thet]=s4point
%this function calculates the cubic coefficient
%in free energy expansion

%which chain? (1=Gaussian, 2=WLC, 3=Rigid Rod)
CHAIN=1;
ORDEig=5;
ORDL=5;
NumLayer=500;

%parameters :: polymer
N=5;
NM=1.;
FA=1.;
LAM=0.99;
d=3;

%parameters :: wave vector
NQ=100;NT=100;
Qmag=logspace(0,2,NQ)/N/NM;
pert=0.;
Thet=linspace(-2*pi+pert,2*pi-pert,NT);
s4aaaa=zeros(NQ,NT);
s4aaab=zeros(NQ,NT);
s4aaba=zeros(NQ,NT);
s4abaa=zeros(NQ,NT);

for I=1:NT
    I
    for J=1:NQ
        T=Thet(I);
        Q=Qmag(J);
        
        s4=s4calc(Q,T,N,NM,LAM,FA,d,CHAIN,ORDEig,ORDL,NumLayer);
        s4aaaa(I,J)=s4(1,1,1,1);
        s4aaab(I,J)=s4(1,1,1,2);
        s4aaba(I,J)=s4(1,1,2,1);
        s4abaa(I,J)=s4(2,1,1,1);
    end
end

%make a plot
figure;
if NT==1
    set(gca,'fontsize',15)
    plot(Qmag*N*NM,s4aaaa)
    xlabel('Q')
    zlabel('\Gamma')
else
    set(gca,'fontsize',15)
    surf(Qmag*N*NM,Thet,s4aaaa,'edgecolor','none','LineStyle','none','FaceLighting','phong');
    set(gca,'xscale','log')
    view([-90,90])
    xlabel('kL')
    ylabel('\theta')
    ylim([min(Thet),max(Thet)])
    zlabel('S_{123}')
    set(gca, 'CLim', [0,1.]);
    colorbar
end
% 
% function s4=s4calc(Q,T,N,NM,LAM,FA,d,CHAIN,ORDEig,ORDL,NumLayer)
%     %wave vectors
%     Q1=Q*[1,0,0];
%     Q2=transpose(rotz(T)*Q1(1:3)');
%     Q3=-Q1;
%     Q4=-Q2;
%     
%     if CHAIN==1
%         s4 = s4gaussian(N,NM,LAM,FA,Q1,Q2,Q3,Q4,d);
%     elseif CHAIN==2
%         s4 = s4wlc(N,NM,LAM,FA,Q1,Q2,Q3,Q4,ORDEig,ORDL,NumLayer);
%     elseif CHAIN==3
%         s4 = s4rigid(N,NM,LAM,FA,Q1,Q2,Q3,Q4,1);
%     end
%     s4=s4/power(N*NM,4);
%     
% end