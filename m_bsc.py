import numpy as np
from LDPC_lib import LDPC_Code
import matplotlib.pyplot as plt

import sys

np.random.seed(0)

K = 10000 # information block size
pckts = 100# number of packets
pts_de = 501
max_ite = 150
max_pdf = 100

print("Points DE =" +str(pts_de),flush=True)
print("Max Ite DE ="+str(max_ite),flush=True)
print("Max PDF DE="+str(max_pdf),flush=True)
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
print("Pflip ="+str(p_vector))
simuber = []
deber=[]
for pflip in p_vector:
    deber.append(LDPC.densityevolution(channel='BSC',pflip=pflip,pts=pts_de,max_ite=max_ite,max_pdf=max_pdf))

plt.figure(1)
plt.semilogy(p_vector,np.array(deber),'k--x',label="Density Evolution")
plt.grid()
plt.title('LDPC - Regular (3,6) - K = '+str(K))

for pflip in p_vector:
    ers = []
    for p in range(pckts):
        pck_ber = LDPC.SimulateBSCPacket(pflip,maxite=200)
        ers.append(pck_ber)
        print("Packet ="+str(p)+" Pflip ="+str(pflip)+" Packet BER=" + str(np.mean(pck_ber)),flush=True)
        print("Ensemble BER=" + str(np.mean(ers)),flush=True)
    print("BSC BER=" + str(np.mean(ers)+" pflip="+str(pflip)),flush=True)
    simuber.append(np.mean(ers))

print(pflip_capacity)
print(p_vector)
print(np.array(deber))
print(np.array(simuber))

plt.semilogy(p_vector,np.array(simuber),'r--^',label="Simulation")
plt.xlabel('Crossover probability p')
plt.ylabel('BER')
# plt.show()
plt.legend(loc="upper left")

plt.savefig('bsc_ber.png')

