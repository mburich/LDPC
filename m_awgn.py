import numpy as np
from LDPC_lib import LDPC_Code
import matplotlib.pyplot as plt

import sys

np.random.seed(0)

K = 10000 # information block size
pckts = 100# number of packets

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

print("Rate=" + str(LDPC.CodeRate) + " N=" + str(LDPC.N) + " K=" + str(LDPC.K),flush=True)
Capacity_awgn = 10 * np.log10((2 ** (2 * LDPC.CodeRate) - 1) / (LDPC.CodeRate * 2))
print("Capacity_awgn=" + str(Capacity_awgn) + 'dB',flush=True)

EbN0gap_vector = np.linspace(0.8,1.4,20)

simuber = []
deber=[]
for EbN0gap in EbN0gap_vector:
    EbN0 = EbN0gap+Capacity_awgn
    N0 = (10 ** (-(EbN0) / 10)) / LDPC.CodeRate
    print("EbN0 = " + str(EbN0),flush=True)
    print("N0 = " + str(N0),flush=True)
    print("EbN0gap = " + str(EbN0gap),flush=True)
    deber.append(LDPC.densityevolution(channel='AWGN',N0=N0,pts=501,max_ite=150,max_pdf=100))

plt.figure(1)
plt.semilogy(EbN0gap_vector,np.array(deber),'k--x',label="Density Evolution")
plt.grid()
plt.title('LDPC - Regular (3,6) - K = '+str(K))

for EbN0gap in EbN0gap_vector:
    EbN0 = EbN0gap+Capacity_awgn
    N0 = (10 ** (-(EbN0) / 10)) / LDPC.CodeRate
    print("EbN0 = " + str(EbN0),flush=True)
    print("N0 = " + str(N0),flush=True)
    print("EbN0gap = " + str(EbN0gap),flush=True)

    ers = []
    for p in range(pckts):
        pck_ber = LDPC.SimulateAWGNPacket(N0,maxite=200)
        ers.append(pck_ber)
        print("Packet ="+str(p)+" EbN0gap = " + str(EbN0gap)+" Packet BER=" + str(np.mean(pck_ber)),flush=True)
        print("Ensemble BER=" + str(np.mean(ers)),flush=True)
    print("AWGN BER=" + str(np.mean(ers)),flush=True)
    simuber.append(np.mean(ers))

print(Capacity_awgn)
print(EbN0gap_vector)
print(np.array(deber))
print(np.array(simuber))


plt.semilogy(EbN0gap_vector,np.array(simuber),'r--^',label="Simulation")
plt.xlabel('EbN0 [dB]')
plt.ylabel('BER')
plt.legend(loc="upper left")
# plt.show()
plt.savefig('awgn_ber.png')

