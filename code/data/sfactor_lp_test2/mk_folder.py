from numpy import *
import os
LPs=loadtxt('LPs')

for lp in LPs:
    os.system('mkdir LP_' + str(lp))
