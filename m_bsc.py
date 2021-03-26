import shelve

import numpy as np
from LDPC_lib import LDPC_Code
import matplotlib.pyplot as plt
import sys


K = 10000 # information block size
pckts = 20# number of packets

print("K = " + str(K),flush=True)
print("Packets = " + str(pckts),flush=True)


DegBits = {}
DegBits['deg'] = np.array([3,4], dtype=np.int32)
DegBits['dist'] = np.array([0.5,0.5])
DegBits['deg'] = np.array([3], dtype=np.int32)
DegBits['dist'] = np.array([1.0])

DegChk = {}
DegChk['deg'] = np.array([4,5], dtype=np.int32)
DegChk['dist'] = np.array([0.5,0.5])

DegChk['deg'] = np.array([6], dtype=np.int32)
DegChk['dist'] = np.array([1])

LDPC = LDPC_Code(DegBits, DegChk, K)

prb_entropy = np.linspace(0,0.5,50001)
for prb_idx in range(len(prb_entropy)-1,0,-1):
    pr = prb_entropy[prb_idx]
    capacity_bsc = 1 - (-pr * np.log2(pr) - (1 - pr) * np.log2(1 - pr))
    if LDPC.CodeRate<capacity_bsc:
        pflip_capacity = pr
        break
print("Rate=" + str(LDPC.CodeRate) + " N=" + str(LDPC.N) + " K=" + str(LDPC.K),flush=True)
print("Capacity -> Pflip=" + str(pflip_capacity),flush=True)

p_vector = np.linspace(0.05,0.12,20)
simuber = []
deber=[]
for pflip in p_vector:
    deber.append(LDPC.densityevolution(channel='BSC',pflip=pflip,pts=401,max_ite=100,max_pdf=100))

plt.figure(1)
plt.semilogy(p_vector,np.array(deber),'k--x')
plt.grid()
plt.title('LDPC - Regular (3,6) - K = '+str(K))

for pflip in p_vector:
    ers = []
    for p in range(pckts):
        pck_ber = LDPC.SimulateBSCPacket(pflip,maxite=200)
        ers.append(pck_ber)
        print("Packet BER=" + str(np.mean(pck_ber)),flush=True)
        print("Ensemble BER=" + str(np.mean(ers)),flush=True)
    print("AWGN BER=" + str(np.mean(ers)),flush=True)
    simuber.append(np.mean(ers))
plt.semilogy(p_vector,np.array(simuber),'r--^')
plt.xlabel('Crossover probability p')
plt.ylabel('BER')
# plt.show()
plt.savefig('bsc_ber.png')

##
filename ='./shelve_bsc.out'
my_shelf = shelve.open(filename,'n') # 'n' for new

for key in dir():
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()
