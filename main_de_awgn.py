import numpy as np
from LDPC_lib import LDPC_Code
import matplotlib.pyplot as plt

import sys

np.random.seed(0)

pts_de = 2501
max_ite = 150
max_pdf = 40

print("Points DE =" +str(pts_de),flush=True)
print("Max Ite DE ="+str(max_ite),flush=True)
print("Max PDF DE ="+str(max_pdf),flush=True)
print("Res="+str(2*max_pdf/pts_de))
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

LDPC = LDPC_Code(DegBits, DegChk)

print("Rate=" + str(LDPC.CodeRate) + " N=" + str(LDPC.N) + " K=" + str(LDPC.K),flush=True)
Capacity_awgn = 10 * np.log10((2 ** (2 * LDPC.CodeRate) - 1) / (LDPC.CodeRate * 2))
print("Capacity_awgn=" + str(Capacity_awgn) + 'dB',flush=True)

EbN0gap_vector = np.linspace(0.8,1.4,20)
# p_vector = np.array([p_vector[-1]])

print("EbN0gap_vector ="+str(EbN0gap_vector))

deber=[]
for EbN0gap in EbN0gap_vector:
    EbN0 = EbN0gap+Capacity_awgn
    N0 = (10 ** (-(EbN0) / 10)) / LDPC.CodeRate
    print("EbN0 = " + str(EbN0),flush=True)
    print("N0 = " + str(N0),flush=True)
    print("EbN0gap = " + str(EbN0gap),flush=True)
    deber.append(LDPC.densityevolution(channel='AWGN',N0=N0,pts=pts_de,max_ite=max_ite,max_pdf=max_pdf))


print("Capacity_awgn=" + str(Capacity_awgn) + 'dB',flush=True)
print("EbN0gap ="+str(EbN0gap_vector))
print("EbN0 ="+str(EbN0gap_vector+Capacity_awgn))

print("AWGN Density Evolution BER="+str(np.array(deber)))
