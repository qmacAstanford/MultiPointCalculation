clear;
close all

k = load('data/N5NM1/k');
s = load('data/N5NM1/s');
t = load('data/N5NM1/t');

figure;surface(k,t,s);
set(gca,'xscale','log')
view([-90,90])
set(gca, 'CLim', [0,1.]);
colorbar
ylim([min(t),max(t)])
xlim([min(k),max(k)])

k = load('data/sfactor_lp_test2/LP_1.0/k2');
s = load('data/sfactor_lp_test2/LP_1.0/S2');
t = load('data/sfactor_lp_test2/LP_1.0/theta2');

figure;surface(k,t,s);
view([-90,90])