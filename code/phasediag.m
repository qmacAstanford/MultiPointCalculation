clear
close all

CHIV=linspace(0.01,0.1,10);

% TEST 1 :: lamellar phase
load('data/Gam_GS_N10NM10_NQ2.mat');
load('data/test2.mat');
B2=B;
C2=C;
C2=C2-min(C2);

% TEST 2 :: cylinder phase
load('data/Gam_GS_N10NM10_NQ3.mat');
load('data/test3.mat');
B3=B;
C3=C;
C3=C3-min(C3);

% TEST 3 :: cylinder phase
load('data/Gam_GS_N10NM10_NQ6.mat');
load('data/test6.mat');
B6=B;
C6=C;
C6=C6-min(C6);

options = optimset('Display','iter','TolX',1e-8,'MaxIter',1500); % show iterations
NF = length(FAV);

% free energy
for I=NF-1:NF-1
    for J=1:length(CHIV)
    FA=FAV(I);
    CHI=CHIV(J);

    % example at small CHI
    A=-CHI*power(2*pi,-3);

    F2=@(psi2) A*power(psi2,2)+B2(I)*power(psi2,3)+C2(I)*power(psi2,4);
    dF2=@(psi2) 2*A+3*B2(I)*power(psi2,1)+4*C2(I)*power(psi2,2);
    if C2(I)>0
        x2=fzero(dF2,[1e-10,100]);
    else
        x2=0;
    end

    F3=@(psi3) A*power(psi3,2)+B3(I)*power(psi3,3)+C3(I)*power(psi3,4);
    dF3=@(psi3) 2*A+3*B3(I)*power(psi3,1)+4*C3(I)*power(psi3,2);
    if B3(I)>0 && C3(I)>0
        x3=fzero(dF3,[1e-10,10]);
    else
        x3=0;
    end

    F6=@(x) A*power(x,2)+B6(I)*power(x,3)+C6(I)*power(x,4);
    dF6=@(psi6) 2*A+3*B6(I)*power(psi6,1)+4*C6(I)*power(psi6,2);

    if B6(I)>0 && C6(I)>0
        x6=fzero(dF6,[1e-10,10]);
    else
        x6=0;
    end

    % find OOT chi
    f6(J)=F6(x6);
    f2(J)=F2(x2);
    f3(J)=F3(x3);

    % % making plots
    % figure;hold
    % CHI=CHI23(3);
    % x3=CHI23(2);
    % A=-CHI*power(2*pi,-3);
    % F3=@(psi3) A*power(psi3,2)+B3(I)*power(psi3,3)+C3(I)*power(psi3,4);
    % fplot(F3,[0,x3+x3/3],'b-')
    % plot(x3,F3(x3),'r.')
    % 
    % figure;hold
    % CHI=CHI23(3);
    % x2=CHI23(1);
    % A=-CHI*power(2*pi,-3);
    % F2=@(psi2) A*power(psi2,2)+B2(I)*power(psi2,3)+C2(I)*power(psi2,4);
    % fplot(F2,[0,x2+x2/3],'b-')
    % plot(x2,F2(x2),'r.')


    end
end

figure;
plot(CHIV,f2,CHIV,f3,CHIV,f6)
legend('2','3','6')


% figure;
% plot(FAV(1:end-1),x2)
% 
% figure;
% plot(FAV(1:end-1),f2,FAV(1:end-1),f6)