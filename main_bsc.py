import numpy as np
from LDPC_lib import LDPC_Code

# import sys
# if  len(sys.argv)!=4:
#     print("Not enough arguments",flush=True)
#     sys.exit()
#
# pflipgap = float(sys.argv[1]) #gap to capacity
# K = int(sys.argv[2])  # information block size
# pckts = int(sys.argv[3])  # number of packets
#

pflipgap = -0.02 #gap to capacity
K = 1000 # information block size
pckts = 1000


print("pflipgap = " + str(pflipgap),flush=True)
print("K = " + str(K),flush=True)
print("Packets = " + str(pckts),flush=True)


DegBits = {}
DegBits['deg'] = np.array([3], dtype=np.int32)
DegBits['dist'] = np.array([1])

DegChk = {}
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

pflip = pflipgap+pflip_capacity

print("pflip = " + str(pflip),flush=True)


ers = []
for p in range(pckts):
    pck_ber = LDPC.SimulateBSCPacket(pflip)
    ers.append(pck_ber)
    print("Packet BER=" + str(np.mean(pck_ber)),flush=True)
    print("Ensemble BER=" + str(np.mean(ers)),flush=True)
print("BSC BER=" + str(np.mean(ers)),flush=True)

