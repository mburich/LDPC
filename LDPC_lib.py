import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt

rng = default_rng(0)

class Node(object):
    def __init__(self):
        self.ListNodes = []
        self.RxBuf = {}

    @staticmethod
    def connect(node1, node2):
        if Node.areconnected(node1, node2):
            return
        node1.ListNodes.append(node2)
        node2.ListNodes.append(node1)

        node1.RxBuf[node2] = 0
        node2.RxBuf[node1] = 0

    @staticmethod
    def areconnected(node1, node2):
        if node1 in node2.ListNodes or node2 in node1.ListNodes:
            print("Already connected",flush=True)
            return True
        else:
            return False

    def cleanmsgs(self):
        for n in self.ListNodes:
            self.RxBuf[n] = 0

    def degree(self):
        return len(self.ListNodes)

    def getmsg(self, node):
        return self.RxBuf[node]

    def sendmsg(self, node, msg):
        node.RxBuf[self] = msg


class Bit(Node):
    def __init__(self):
        Node.__init__(self)
        self.LLR = None
        self.trueval = None

    def computemsg(self, edge):
        msg = np.array(self.LLR)
        for l in self.ListNodes:
            if edge != l:
                msg += self.getmsg(l)

        return msg

    def decision(self):
        msg = np.array(self.LLR)
        for l in self.ListNodes:
            msg += self.getmsg(l)
        val = 0 if msg > 0 else 1

        # if self.trueval == None:
        #     a = 2

        if self.trueval != val:
            return 1
        else:
            return 0


class Chk(Node):
    def __init__(self):
        Node.__init__(self)

    def computemsg(self, edge):
        msg = 1.0
        for l in self.ListNodes:
            if edge != l:
                s = np.tanh(self.getmsg(l) / 2.0)
                if np.abs(s) > 0.9999999999999999:
                    s = np.sign(s) * 0.9999999999999999  # highest precision
                msg = msg * s

        if np.abs(msg) > 0.9999999999999999:
            msg = np.sign(msg) * 0.9999999999999999  # highest precision
        msg = 2.0 * np.arctanh(msg)

        return msg


# if np.abs(s) > 0.9999999999999999:
#     s = np.sign(s) * 0.9999999999999999  # highest precision
# s = 2.0 * np.arctanh(s)

class LDPC_Code:

    @staticmethod
    def matmul(A, B):
        if A.shape[1] != B.shape[0]:
            return np.array([])
        C = np.zeros((A.shape[0], B.shape[1]))
        for i in range(C.shape[0]):
            for j in range(C.shape[1]):
                C[i, j] = np.sum(A[i, :] * B[:, j]) % 2
        return C

    @staticmethod
    def gausselimination(H):
        colswap = []
        print("Building Hsys",flush=True)

        def swaprow(H, i, j):
            st = np.array(H[i, :])
            H[i, :] = H[j, :]
            H[j, :] = st
            # print("R"+str(i)+"<->R"+str(j))
            # print(H)

        def swapcol(H, i, j):
            st = np.array(H[:, i])
            H[:, i] = H[:, j]
            H[:, j] = st
            colswap.append((i, j))
            # print("C"+str(i)+"<->C"+str(j))
            # print(H)

        Hori = np.array(H)
        H = np.array(H)

        for row in range(H.shape[0]):  # find canonical form for row
            if row % 100 == 0:
                print("Pivot=" + str(row),flush=True)

            for col_t in range(row, H.shape[1]):  # look for columns candidate to be swapped i.e. zeros above row
                idx = np.where(H[:, col_t] == 1)[0]  # full matrix
                # idx = find(H[:,col_t]==1)[0] # sparse matrix
                idx = idx[idx >= row]
                if len(idx) != 0:
                    if col_t != row:
                        swapcol(H, row, col_t)
                    break

            if len(idx) == 0:
                print("singular?",flush=True)
                return (np.array([]), np.array([]), np.array([]))

            if idx[0] != row:
                swaprow(H, idx[0], row)
            for k in range(H.shape[0]):  # put zeros under pivot
                if H[k, row] == 1 and k != row:
                    H[k, :] = np.abs(H[k, :] - H[row, :])
                    # print("R"+str(k)+"+R"+str(col)+"->R"+str(k))
                    # print(H)

        for k in range(H.shape[0]-1,-1,-1):
            idx = H.shape[1]-H.shape[0]+k
            swapcol(H, idx, k)
        # for k in range(H.shape[0]):
        #     idx = k + H.shape[1] - H.shape[0]
        #     swapcol(H, idx, k)

        colswap_copy = colswap.copy()
        for k in colswap_copy:
            swapcol(Hori, k[0], k[1])
        return (H, Hori, colswap_copy)  # return systematic and original matrix with appropriate column shift

    @staticmethod
    def listdegs(Profile, Size):
        DegList = np.zeros(0)
        for d, p in zip(Profile['deg'], Profile['dist']):
            n = round(p * Size)
            DegList = np.concatenate((DegList, d * np.ones(n, dtype=np.int))) if len(DegList) != 0 else d * np.ones(n,
                                                                                                                    dtype=np.int)
        rng.shuffle(DegList[:Size])
        return DegList[:Size]

    def __init__(self, DegBits, DegChk, K=10000):
        self.G = None
        self.H = None
        self.tanhtable = None
        self.DegBits = DegBits
        self.DegChk = DegChk
        self.AvgDegBits = np.sum(DegBits['deg'] * DegBits['dist'])
        self.AvgDegChk = np.sum(DegChk['deg'] * DegChk['dist'])
        self.CodeRate = (1 - self.AvgDegBits / self.AvgDegChk)

        self.K = K
        self.N = round(self.K / self.CodeRate)
        self.CodeRate = self.K / self.N

        self.NumChk = self.N - self.K

    @staticmethod
    def xpos_to_idx(x,max_pdf,x_res):
        return (np.round((x+max_pdf)/(x_res))).astype(int)

    @staticmethod
    def buildtanhtable(x_pdf,max_pdf,x_res):
        print("Building Tanh Table")
        t = (np.tanh(x_pdf / 2.0))
        x, y = np.meshgrid(range(len(x_pdf)), range(len(x_pdf)))
        z = t[x]*t[y]
        idx_z = tuple([np.abs(z) > 0.9999999999999999])
        z[idx_z] = np.sign(z[idx_z])*0.9999999999999999
        tanhtable = LDPC_Code.xpos_to_idx(2.0*np.arctanh(z),max_pdf,x_res)
        pos_idx = np.unique(tanhtable)
        pos_dict = {}
        for i in range(tanhtable.shape[0]):
            for j in range(tanhtable.shape[1]):
                xy = tanhtable [i,j]
                if not xy in pos_dict.keys():
                    pos_dict[xy] = [[i,j]]
                else:
                    pos_dict[xy].append([i, j])
        for i in pos_dict.keys():
            pos_dict[i] = np.array(pos_dict[i])
        # for i in pos_idx:
        #     pos_dict[i] = np.where(tanhtable == i)
        print("Finished Tanh Table")

        return tanhtable,pos_dict

    def densityevolution(self,channel='AWGN',N0=None,pflip=None,max_pdf=50,pts=100,max_ite=100):

        def normalizepdf(x):
            return x/(np.sum(x))




        def slowbuildtanhtable(x_pdf):
            print("Building Slow Tanh Table")
            c = np.zeros((len(x_pdf),len(x_pdf)),dtype=np.int)
            for ka in range(len(x_pdf)):
                for kb in range(ka,len(x_pdf)):
                    s = (np.tanh(x_pdf[ka] / 2.0) * np.tanh(x_pdf[kb] / 2.0))
                    if np.abs(s) > 0.9999999999999999:
                        s = np.sign(s) * 0.9999999999999999  # highest precision
                    s = 2.0 * np.arctanh(s)
                    idx = LDPC_Code.xpos_to_idx(s,max_pdf,x_res)
                    c[ka, kb] = idx
                    c[kb, ka] = idx
            return c

        if self.AvgDegBits*self.N!=self.AvgDegChk*self.NumChk:
            print('Invalid code')
            return -1


        x_pdf = np.linspace(-max_pdf,max_pdf,pts)
        x_res = x_pdf[1]-x_pdf[0]



        if self.tanhtable is None:
            self.tanhtable = LDPC_Code.buildtanhtable(x_pdf,max_pdf,x_res)

        def slowtanhdensity(a,b):
            c = 0*x_pdf
            for ka in range(len(x_pdf)):
                for kb in range(len(x_pdf)):
                    prob = a[ka]*b[kb]
                    if prob<1e-16:
                        continue
                    idx = self.tanhtable[0][ka,kb]
                    c[idx] = prob+c[idx]
            return normalizepdf(c)

        def tanhdensity(a,b):
            pos_dict = self.tanhtable[1]
            c = 0*x_pdf
            prob_vec = a.reshape(-1,1)@b.reshape(1,-1)
            for ka in pos_dict.keys():
                idx = ka
                cords = pos_dict[ka]
                c[idx] = np.sum(prob_vec[cords[:,0],cords[:,1]])

                # x_cord = pos_dict[ka][0]
                # y_cord = pos_dict[ka][1]
                # c[idx] = np.sum(prob_vec[x_cord,y_cord])
            return normalizepdf(c)

        if channel=='AWGN':
        #ch_msg = 2y/n_var -> y = x+n -> x = -1 always [all zeros] -> mean(ch_msg) = +2/n_var -> var(ch_msg) = 4/n_var
            ch_pdf = np.exp(-((x_pdf - 4.0 / N0) ** 2) / (16 / N0))
            ch_pdf = normalizepdf(ch_pdf)


        elif channel=='BSC':
        # chmsg = np.log((1.0 - pflip)/pflip) if rcv =1 , np.log(pflip / (1.0 - pflip)) if 0
        # sending all-zeros -> w.p. pflip we get np.log(pflip / (1.0 - pflip)) -> w.p. (1-p) we get -np.log(pflip / (1.0 - pflip))
            ch_pdf = 0*x_pdf
            xp = np.log((1-pflip) / pflip)
            ch_pdf[LDPC_Code.xpos_to_idx(xp,max_pdf,x_res)] = 1-pflip
            ch_pdf[LDPC_Code.xpos_to_idx(-xp,max_pdf,x_res)] = pflip
        else:
            return

        var_to_chk = np.array(ch_pdf)
        # plt.figure(1)
        for ite in range(max_ite):
            # chk to var
            chk_ite = 2
            z = np.array(var_to_chk)
            # plt.plot(x_pdf,var_to_chk)
            # plt.show()

            chk_to_var = 0 * x_pdf
            for chk_dist in range(len(self.DegChk['deg'])):
                while chk_ite != self.DegChk['deg'][chk_dist]:
                    z = tanhdensity(var_to_chk,z)
                    chk_ite+=1
                chk_to_var += np.array(self.DegChk['dist'][chk_dist] * z)
            chk_to_var = normalizepdf(chk_to_var)

            # plt.plot(x_pdf,chk_to_var)
            # plt.show()
            # var to chk
            bits_ite = 1
            z = np.array(ch_pdf)
            var_to_chk = 0*x_pdf
            for var_dist in range(len(self.DegBits['deg'])):
                while bits_ite != self.DegBits['deg'][var_dist]:
                    z = normalizepdf(np.convolve(chk_to_var,z,mode='same'))
                    bits_ite+=1
                var_to_chk += np.array(self.DegBits['dist'][var_dist] * z)
            var_to_chk = normalizepdf(var_to_chk)


            decision = normalizepdf(np.convolve(var_to_chk,chk_to_var,mode='same'))
            ber = np.sum(decision[x_pdf<0])

            print("BER = "+str(ber)+" iteration="+str(ite))

            if ber<1e-6:
                print("Converged at iteration=" + str(ite))
                return 0
        print("Didn't converge at iteration=" + str(ite))
        return ber

    def null_check(self):
        null_check = LDPC_Code.matmul(self.G, self.Hadj.transpose())
        if np.all(null_check == 0):
            print("Null check pass")
            return 1
        else:
            print("Null check problem")
            return 0

    def buildcode(self):
        print(" Rate=" + str(self.CodeRate) + " N=" + str(self.N) + " K=" + str(self.K),flush=True)

        Hsys = np.array([])
        while len(Hsys) == 0:
            H = LDPC_Code.buildH(self.DegBits, self.DegChk, self.N, self.N - self.K)
            while len(H) == 0:
                H = LDPC_Code.buildH(self.DegBits, self.DegChk, self.N, self.N - self.K)
            (Hsys, Hadj, colswap) = LDPC_Code.gausselimination(H)

        P = np.transpose(Hsys[:, :self.K])
        G = np.hstack((np.eye(self.K), P))

        self.H = H
        self.Hadj = Hadj
        self.Hsys = Hsys
        self.G = G

        print("Parity Check Matrix",flush=True)
        print(self.H,flush=True)
        print("Parity Check Matrix - Adj",flush=True)
        print(self.Hadj,flush=True)
        print("Parity Check Matrix - Echelon form",flush=True)
        print(self.Hsys,flush=True)

        bit_deg = np.unique(np.sum(self.Hadj, 0))
        print("Bit deg = " + str(bit_deg),flush=True)
        chk_deg = np.unique(np.sum(self.Hadj, 1))
        print("Chk deg = " + str(chk_deg),flush=True)



        (self.ChkNodes, self.BitNodes) = self.buildgraph(Hadj)

    def cleanmsgs(self):
        for bit in self.BitNodes:
            bit.cleanmsgs()

        for chk in self.ChkNodes:
            chk.cleanmsgs()

    def encode(self, bits, channel):
        self.cleanmsgs()

        if bits.shape[0] != self.G.shape[0]:
            return np.array([])

        # for c in range(bits.shape[0]):
        #     self.BitNodes[c].trueval = bits[c]

        y = LDPC_Code.matmul(self.G.transpose(), bits)
        # print(LDPC_Code.matmul(self.Hsys,y))
        print("CW check="+str(np.sum(LDPC_Code.matmul(self.Hadj,y)))) #check if codeword

        for c in range(y.shape[0]):
            self.BitNodes[c].trueval = y[c]

        for chk in range(len(self.ChkNodes)):
            st=0
            for L in range(len(self.ChkNodes[chk].ListNodes)):
                st+=self.ChkNodes[chk].ListNodes[L].trueval
            if st%2!=0:
                print("Invalid chk node",flush=True)
                return np.nan
        print("Graph cw ok",flush=True)

        syndrome = np.sum(LDPC_Code.matmul(self.Hadj, y))
        if syndrome != 0.0:
            print("Invalid CW",flush=True)
            return np.nan
        print("Valid cw",flush=True)


        if channel == 'BSC':
            return y
        elif channel == 'AWGN':
            return 2 * y - 1  # energy 1
        return -1

    def decode(self, codeword, channel, n_variance=None, pflip=None, maxite=400):
        if channel == 'BSC':
            if pflip != None:
                for c in range(codeword.shape[0]):
                    if codeword[c] == 1:
                        self.BitNodes[c].LLR = -np.log((1-pflip) / pflip)
                    else:
                        self.BitNodes[c].LLR = +np.log((1-pflip) / pflip)
            else:
                print("insert pflip",flush=True)
                return

        if channel == 'AWGN':
            if n_variance != None:
                for c in range(codeword.shape[0]):
                    self.BitNodes[c].LLR = -2.0*codeword[c]/(n_variance)
            else:
                print("insert awgn",flush=True)
                return

        self.cleanmsgs()

        rpt = np.zeros(len(self.BitNodes))
        rpt_prev = np.zeros(len(self.BitNodes))

        for ite in range(maxite):
            for bit in self.BitNodes:
                for edge in bit.ListNodes:
                    msg = bit.computemsg(edge)
                    bit.sendmsg(edge, msg)
            for chk in self.ChkNodes:
                for edge in chk.ListNodes:
                    msg = chk.computemsg(edge)
                    chk.sendmsg(edge, msg)
            er = 0
            bidx = 0
            for bit in self.BitNodes:

                bt = bit.decision()
                if ite == 0:
                    rpt_prev[bidx] = bt
                else:
                    if (rpt_prev[bidx] == bt):
                        rpt[bidx] += 1
                    else:
                        rpt_prev[bidx] = bt
                        rpt[bidx] = 0
                bidx += 1
                er += bt

            print("Ite=" + str(ite) + " with " + str(er) + " errors",flush=True)
            if min(rpt) > 10:
                return er / len(self.BitNodes)

        return er / len(self.BitNodes)

    @staticmethod
    def buildH(DegBits, DegChk, CW_Size, NumChk, max_trials=100000):
        print("Building H",flush=True)
        ListDegBits = LDPC_Code.listdegs(DegBits, CW_Size)
        ListDegBits_cur = 0 * ListDegBits
        ListDegChk = LDPC_Code.listdegs(DegChk, NumChk)
        AvgDegbits = np.sum(DegBits['deg'] * DegBits['dist'])
        AvgDegChk = np.sum(DegChk['deg'] * DegChk['dist'])

        numconnections = CW_Size * AvgDegbits
        connections = 0
        H = np.zeros((NumChk, CW_Size), dtype=np.int8)

        for chk in range(len(ListDegChk)):
            for edges in range(ListDegChk[chk]):

                bit_v = np.where((H[chk, :] != 1) & (ListDegBits_cur != ListDegBits))[0]
                if len(bit_v) == 0:
                    return np.array([])
                rnd_idx = rng.integers(len(bit_v))

                bit = bit_v[rnd_idx]
                H[chk, bit] = 1
                print("Chk=" + str(chk) + " Bit=" + str(bit),flush=True)
                ListDegBits_cur[bit] = ListDegBits_cur[bit] + 1

                connections += 1
                if connections == numconnections:
                    if np.mean(np.sum(H,1))!=AvgDegChk:
                        values, counts = np.unique(np.sum(H, 1), return_counts=True)
                        DegChk['deg'] = values
                        DegChk['dist'] = counts/np.sum(counts)
                        print('Updating Check Node profile to Deg='+str(DegChk['deg'])+" Dist="+str(DegChk['dist']))

                    return H
        return H

    def buildgraph(self, H):
        ChkNodes = [Chk() for i in range(self.NumChk)]
        BitNodes = [Bit() for i in range(self.N)]

        for chk in range(H.shape[0]):
            for bit in range(H.shape[1]):
                if H[chk, bit] == 1:
                    Node.connect(BitNodes[bit], ChkNodes[chk])

        return (ChkNodes, BitNodes)

    def SimulateAWGNPacket(self,N0,maxite=100):
        if self.G is None:
            self.buildcode() # we need to build G and H first

        x = ((rng.uniform(0, 1, self.K) > 0.5) * 1).reshape(-1, 1)
        y_awgn = self.encode(x, channel='AWGN')
        n_awgn = rng.normal(0, np.sqrt(N0 / 2.0), size=y_awgn.shape[0]).reshape(-1, 1)
        r_awgn = y_awgn + n_awgn
        pck_ber = self.decode(r_awgn, channel='AWGN', n_variance=N0 / 2,maxite=maxite)
        return pck_ber

    def SimulateBSCPacket(self,pflip,maxite=100):
        if self.G is None:
            self.buildcode() # we need to build G and H first

        x = ((rng.uniform(0, 1, self.K) > 0.5) * 1).reshape(-1, 1)
        y_bec = self.encode(x,channel='BSC')
        n_bec =  ((rng.uniform(0,1,y_bec.shape[0])<pflip)*1).reshape(-1,1)
        r_bec = (y_bec+n_bec)%2
        pck_ber = self.decode(r_bec,channel='BSC',pflip=pflip,maxite=maxite)
        return pck_ber
