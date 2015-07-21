%this code tests the calculation of 2-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains

clear;

%Which chains?
WLC_ON=0;
GS_ON=1;
RR_ON=0; %read from pre-calculated if RR_ON=2

%Chain structural information
N=100;
NM_WL=logspace(-1,2,6);
NM_GS =100;
NM_RR =1;

%Chain chemical information
FA=1.;
LAM=.99;

%Calculation parameters
d=3;
ORDmax=20;
ORD=20;
ResLayer=500;

%wavevector
k=logspace(-2,3,100)';

figure;hold;set(gca,'fontsize',15);leg=[];
%%%% WLC %%%%
if (WLC_ON==1)
    cnt=1;
    leg=[];
    s2aa=zeros(length(k),1);
    s2ab=zeros(length(k),1);
    s2bb=zeros(length(k),1);
    for NM_WLC=NM_WL
        col=(cnt-1)/(length(NM_WL)-1)
        
        for ii=1:length(k)
            s2 = s2wlc(N,NM_WLC,LAM,FA,k(ii)/N/NM_WLC,d,ORDmax,ORD,ResLayer);
            
            s2=s2/power(N*NM_WLC,2);
            s2aa(ii) = s2(1,1);
            s2ab(ii) = s2(1,2);
            s2bb(ii) = s2(2,2);
        end
        plot(k,s2aa,k,s2ab,'linewidth',2,'color',[0 0 1-col]);
    
%         plot(k,s2aa,'b-','linewidth',2,'color',[0 0 1-col]);
%         leg=[leg {['WLC NM=',sprintf('%.2f',NM_WLC)]}];
        cnt=cnt+1;
    end
end

%%%% Gaussian Chain %%%%
if (GS_ON==1)
    g2aa=zeros(length(k),1);
    g2ab=zeros(length(k),1);
    g2bb=zeros(length(k),1);
    for ii=1:length(k)
        g2 = s2gaussian(N,NM_GS,LAM,FA,k(ii)/N/NM_GS,d);
        g2 = g2./power(N*NM_GS,2);
        
        g2aa(ii) = g2(1,1);
        g2ab(ii) = g2(1,2);
        g2bb(ii) = g2(2,2);
    end
    plot(k,g2aa,'k-',k,g2ab,'k--','linewidth',2);
    leg=[leg {['Gaussian S_{AA}']}];
    leg=[leg {['Gaussian S_{AB}']}];

%     plot(k,g2aa/power(N*NM_GS,2),'k-','linewidth',2);
%     leg=[leg {['Gaussian']}];
end

%%%% Rigid Rod %%%%
if (RR_ON==1)
    r2aa=zeros(length(k),1);
    r2ab=zeros(length(k),1);
    r2bb=zeros(length(k),1);
    for ii=1:length(k)
        r2 = s2rigid(N,NM_RR,LAM,FA,k(ii)/N/NM_RR,1);
        r2 = r2./power(N*NM_RR,2);
        
        r2aa(ii) = r2(1,1);
        r2ab(ii) = r2(1,2);
        r2bb(ii) = r2(2,2);
    end
    plot(k,r2aa,'b-',...
        k,r2ab,'b--','linewidth',2);
    leg=[leg {'Rigid Rod S_{AA} '}];
    leg=[leg {'Rigid Rod S_{AB} '}];

%     plot(k,r2aa/power(N*NM_RR,2),'b-','linewidth',2);
%     leg=[leg {'Rigid Rod'}];
elseif (RR_ON==2)
    r2aa = load('data/R2_N100NM1');
    plot(r2aa(:,1),r2aa(:,2),'b-','linewidth',2);
    leg=[leg {'Rigid Rod'}];
end

legend(leg);
set(gca,'xscale','log');set(gca,'yscale','linear');ylim([0,1]);
xlabel('Normalized Wavevector kL');ylabel('Normalized Structure Factor S/L^2')