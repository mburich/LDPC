import numpy as np
from LDPC_lib import LDPC_Code

# import sys
# if  len(sys.argv)!=4:
#     print("Not enough arguments",flush=True)
#     sys.exit()
#
# EbN0gap = float(sys.argv[1]) # gap in dBs to capacity
# K = int(sys.argv[2])  # information block size
# pckts = int(sys.argv[3])  # number of packets

EbN0gap = 1.0 # gap in dBs to capacity
K = 1000 # information block size
pckts = 1000 # number of packets

print("EbN0gap = " + str(EbN0gap),flush=True)
print("K = " + str(K),flush=True)
print("Packets = " + str(pckts),flush=True)


DegBits = {}
DegBits['deg'] = np.array([3], dtype=np.int32)
DegBits['dist'] = np.array([1])

DegChk = {}
DegChk['deg'] = np.array([6], dtype=np.int32)
DegChk['dist'] = np.array([1])

LDPC = LDPC_Code(DegBits, DegChk, K)

print("Rate=" + str(LDPC.CodeRate) + " N=" + str(LDPC.N) + " K=" + str(LDPC.K),flush=True)
Capacity_awgn = 10 * np.log10((2 ** (2 * LDPC.CodeRate) - 1) / (LDPC.CodeRate * 2))
print("Capacity_awgn=" + str(Capacity_awgn) + 'dB',flush=True)

EbN0 = EbN0gap+Capacity_awgn
N0 = (10 ** (-(EbN0) / 10)) / LDPC.CodeRate

print("EbN0 = " + str(EbN0),flush=True)
print("N0 = " + str(N0),flush=True)


ers = []
for p in range(pckts):
    pck_ber = LDPC.SimulateAWGNPacket(N0)
    ers.append(pck_ber)
    print("Packet BER=" + str(np.mean(pck_ber)),flush=True)
    print("Ensemble BER=" + str(np.mean(ers)),flush=True)
print("AWGN BER=" + str(np.mean(ers)),flush=True)

