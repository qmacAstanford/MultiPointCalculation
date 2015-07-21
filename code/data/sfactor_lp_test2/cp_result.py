from numpy import *
import os

LPs=loadtxt('LPs')

for lp in LPs:
    os.system('scp -r brad@tower1:~/sfactor_lp_test_12_16_2014/LP_' + str(lp) +'/k2 LP_' + str(lp) + '/')
    os.system('scp -r brad@tower1:~/sfactor_lp_test_12_16_2014/LP_' + str(lp) +'/theta2 LP_' + str(lp) + '/')
    os.system('scp -r brad@tower1:~/sfactor_lp_test_12_16_2014/LP_' + str(lp) +'/S2 LP_' + str(lp) + '/')
