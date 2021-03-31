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
# p_vector = np.array([p_vector[-1]])

print("Pflip ="+str(p_vector))

deber=[]
for pflip in p_vector:
    deber.append(LDPC.densityevolution(channel='BSC', pflip=pflip, pts=pts_de, max_ite=max_ite, max_pdf=max_pdf))
print("Capacity pflip="+str(pflip_capacity))
print("Pflip ="+str(p_vector))
print("Density Evolution BER="+str(np.array(deber)))
