%this code tests the calculation of 3-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains

clear;

%Which chains?
WLC_ON=0;
GS_ON=1;
RR_ON=0; %read from pre-calculated if RR_ON=2

%Chain structural information
N=100;
NM_WL=logspace(-1,1,2);
NM_GS=10;
NM_RR=1;

%Chain chemical information
FA=.5;
LAM=-.75;
 
%Calculation parameters
d=3;
ORDEig=5;
ORDL=5;
NumLayer=500;

%wavevector and structure factor
QM=logspace(-2,3,500)';
Q1=zeros(length(QM),1);
Q2=zeros(length(QM),1);
Q3=zeros(length(QM),1);
ang=120;
for ii=1:length(QM)
    pert=[1,-1,1]*QM(ii)*1e-2;
    Q1(ii,1:3)=QM(ii)*[1,0,0]+pert;
    Q2(ii,1:3)=transpose(rotz(ang)*Q1(ii,1:3)')-pert;
    Q3(ii,1:3)=-Q1(ii,1:3)-Q2(ii,1:3);
end

%begin making plots
figure;hold;set(gca,'fontsize',15);leg=[];

%%%% WLC %%%%
if (WLC_ON==1)
    cnt=1;
    s3aaa=zeros(length(QM),1);
    s3aab=zeros(length(QM),1);
    s3aba=zeros(length(QM),1);
    s3baa=zeros(length(QM),1);
    for NM_WLC=NM_WL
        col=(cnt-1)/(length(NM_WL)-1)
        for ii=1:length(QM)
            s3=s3wlc(N,NM_WLC,LAM,FA,Q1(ii,1:3)/N/NM_WLC,...
                                     Q2(ii,1:3)/N/NM_WLC,...
                                     Q3(ii,1:3)/N/NM_WLC,ORDEig,ORDL,NumLayer);
            s3=s3/power(N*NM_WLC,3);
            s3aaa(ii)=s3(1,1,1);
            s3aab(ii)=s3(1,1,2);
            s3aba(ii)=s3(1,2,1);
            s3baa(ii)=s3(2,1,1);
        end
        plot(QM,s3aaa,'-','linewidth',2,'color',[0 0 1-col]);
%         plot(QM,s3aaa,'k-',...
%              QM,s3aab,'b-',...
%              QM,s3aba,'rx-',...
%              QM,s3baa,'c-','linewidth',2);
        leg=[leg {['WLC NM=',sprintf('%.2f',NM_WLC)]}];
        cnt=cnt+1;
    end
end

%%%% Gaussian Chain %%%%
if (GS_ON==1)
    g3aaa=zeros(length(QM),1);
    g3aab=zeros(length(QM),1);
    g3aba=zeros(length(QM),1);
    g3baa=zeros(length(QM),1);
    for ii=1:length(QM)
        g3=s3gaussian(N,NM_GS,LAM,FA,Q1(ii,1:3)/N/NM_GS,...
                                     Q2(ii,1:3)/N/NM_GS,...
                                     Q3(ii,1:3)/N/NM_GS,d);
        g3=g3/power(N*NM_GS,3);
        g3aaa(ii)=g3(1,1,1);
        g3aab(ii)=g3(1,1,2);
        g3aba(ii)=g3(1,2,1);
        g3baa(ii)=g3(2,1,1);
    end
%     plot(QM,g3aaa,'k-',...
%          QM,g3aba,'k--','linewidth',2);
%     leg=[leg {['Gaussian S_{AA}']}];
%     leg=[leg {['Gaussian S_{AB}']}];
    plot(QM,g3aaa,'k-','linewidth',2);
    leg=[leg {['Gaussian']}];
end
% 
% %%%% Gaussian Chain :: simple version (based on s2gaussian.m) %%%%
% if (GS_ON==1)
%     for ii=1:length(QM)
%         g3=s3gaussian2(N,NM_GS,LAM,FA,Q1(ii,1:3)/N/NM_GS,...
%                                      Q2(ii,1:3)/N/NM_GS,...
%                                      Q3(ii,1:3)/N/NM_GS,d);
%         g3aaa(ii)=g3(1,1,1);
%     end
%     plot(QM,g3aaa/power(N*NM_GS,4)/FA,'b-','linewidth',2);
%     leg=[leg {['Gaussian NM=',sprintf('%.2f',NM_GS)]}];
% end

%%%% Rigid Rod %%%%
if (RR_ON==1)
    r3aaa=zeros(length(QM),1);
    r3aab=zeros(length(QM),1);
    r3aba=zeros(length(QM),1);
    r3baa=zeros(length(QM),1);
    for ii=1:length(QM)
       ii
       r3=s3rigid(N,NM_RR,LAM,FA,Q1(ii,1:3)/N/NM_RR,...
                                 Q2(ii,1:3)/N/NM_RR,...
                                 Q3(ii,1:3)/N/NM_RR,1);
       r3aaa(ii)=r3(1,1,1);
       r3aab(ii)=r3(1,1,2);
       r3aba(ii)=r3(1,2,1);
       r3baa(ii)=r3(2,1,1);
    end
%     plot(QM,r3aaa/power(N*NM_RR,3),'b-',...
%          QM,r3aba/power(N*NM_RR,3),'b--','linewidth',2);
%     leg=[leg {'Rigid Rod S_{AA} '}];
%     leg=[leg {'Rigid Rod S_{AB} '}];
    plot(QM,r3aaa/power(N*NM_RR,3),'b-','linewidth',2);
    leg=[leg {'Rigid Rod'}];
elseif (RR_ON==2)
    r3aaa = load('data/R3_N100NM1');
    plot(r3aaa(:,1),r3aaa(:,2),'k--','linewidth',2);
    leg=[leg {'Rigid Rod'}];
end

legend(leg);
% title(['f_A = ', num2str(FA), '\lambda = ', num2str(LAM)])
set(gca,'xscale','log');set(gca,'yscale','linear');ylim([0,1.4]);
xlabel('Normalized Wavevector kL');
ylabel('Normalized Structure Factor S/L^3')