import numpy as np
import os
import torch
from torch import nn
import matplotlib
import matplotlib.pyplot as plt
import scipy.io
from scipy.optimize import fsolve

#import file which contains the model for the neural network used to reconstruct the eddy viscosity coefficients:
import TrainNet

from sklearn.model_selection import train_test_split
from scipy.linalg import svdvals
from sklearn import preprocessing
from scipy.interpolate import Rbf
from scipy.interpolate import griddata
from sklearn.utils.extmath import randomized_svd
from matplotlib.ticker import MaxNLocator
import math
from sklearn import *
import cvxpy as cp
plt.rcParams.update({
    "text.usetex": True,
    "font.family": 'STIXGeneral',
    'mathtext.fontset': 'stix'})
matplotlib.rcParams.update({'font.size': 16})

torch.manual_seed(0)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(0)

def readFile(fileName):
    '''
    Function that converts a file of path 'fileName' into a list.
    '''
    fileObj = open(fileName, "r") #opens the file in read mode
    words = fileObj.read().splitlines() #puts the file into an array
    words.pop(0); words.pop(0)
    words = [float(i) for i in words]
    fileObj.close()
    return words

def read_txt(file):
    '''
    Function that converts a text file of path 'file' into a numpy array.
    '''
    f = open(file)
    triplets = f.read().split()
    for i in range(0,len(triplets)):
        triplets[i]=triplets[i].split(',')
    A = np.array(triplets, dtype=np.float32)
    return A

def read_foam_vec(folder):
    '''
    Function that reads a vector field from a Foam file and converts it
    in a numpy array.
    '''
    file = open(folder, 'r')
    cont = file.read().splitlines()
    cont = cont[cont.index('(')+1: cont.index(')')]
    vec = np.zeros(shape=(len(cont), 3))
    for i, line in enumerate(cont):
        line = line[line.index('(')+1:line.index(')')]
        vec[i, :] = np.array([float(number) for number in line.split()])
    file.close()
    return vec

def read_foam_scal(folder):
    '''
    Function that reads a scalar field from a Foam file and converts it
    in a numpy array.
    '''
    file = open(folder, 'r')
    cont = file.read().splitlines()
    cont = cont[cont.index('(')+1: cont.index(')')]
    scal = np.zeros(shape=(len(cont), ))
    for i, line in enumerate(cont):
        scal[i] = float(line)
    file.close()
    return scal


def read_vec_modes(folder, n_points, Nu_modes, field):
    '''
    Function that reads the POD modes of a vector field from the folder created from ITHACA-FV.
    '''
    U_modes = np.zeros((n_points*3, Nu_modes))
    for mode in range(1, Nu_modes+1):
        U_mode = read_foam_vec(os.path.join(folder, str(mode), field))
        U_mode = U_mode.flatten()
        U_modes[:, mode-1] = U_mode
    return U_modes

def read_p_modes(folder, n_points, Np_modes):
    '''
    Function that reads the POD modes of a scalar field from the folder created from ITHACA-FV.
    '''
    p_modes = np.zeros((n_points, Np_modes))
    for mode in range(1, Np_modes+1):
        p_mode = read_foam_scal(os.path.join(folder, str(mode), 'p'))
        p_modes[:, mode-1] = p_mode
    return p_modes



class PPE_ROM():
    def __init__(self, data_folder, _Nu, _Nsup, _Np, _Nnut, _Ru):
        '''
        Class in which a reduced Poisson dynamical system can be solved, with
        or without correction terms and/or turbulence modelling.
        Input parameters: _Nu: number of velocity modes
                          _Nsup: number of supremizer modes
                          _Np: number of pressure modes
                          _Nnut: number of eddy viscosity modes
                          _Ru: number of singular values to retain for the construction
                          of the correction term.
        '''
        #Sizes
        self.data_folder = data_folder
        self.Nu = _Nu
        self.Nsup = _Nsup
        self.Np = _Np
        self.Nnut = _Nnut
        self.M = 2001
        self.M1 = 500
        self.Nu_tot = 50
        self.Nsup_tot = 0
        self.Np_tot = 50
        self.Nnut_tot = 50
        self.Ru = _Ru
        self.timeOnline = []
        self.dtOn = []
        self.ntOn = []
        self.sol = []
        self.Net = []
        self.X_train = []
        self.X_test = []
        self.Y_train = []
        self.Y_test = []
        self.inp_np_train = []
        self.inp_np_test = []
        self.out_np_train = []
        self.out_np_test = []
        self.scaling = []
        self.trainset = []
        self.testset = []
        self.trainloader = []
        self.testloader = []
        self.loss_fn = []
        self.lossplottrain = []
        self.lossplottest = []
        self.tplot = []
        self.device = torch.device('cpu')
        self.gsol = []
        self.BC = 5
        self.tau = 1000
        self.errU_proj = np.zeros(self.M)
        self.errU_full = np.zeros(self.M)
        self.errP_proj = np.zeros(self.M)
        self.errP_full = np.zeros(self.M)

        folder = os.path.join(data_folder, "Matrices")
        B = np.load(os.path.join(folder, "B.npy"))
        Bt = np.load(os.path.join(folder, "bt.npy"))
        M = np.load(os.path.join(folder, "M.npy"))
        BC1 = np.load(os.path.join(folder, "BC1.npy"))
        BC2 = np.load(os.path.join(folder, "BC2.npy"))
        BC3 = np.load(os.path.join(folder, "BC3.npy"))
        C = np.load(os.path.join(folder, "C.npy"))
        K = np.load(os.path.join(folder, "K.npy"))
        G = np.load(os.path.join(folder, "G.npy"))
        D = np.load(os.path.join(folder, "D.npy"))
        P = np.load(os.path.join(folder, "P.npy"))
        Ct1 = np.load(os.path.join(folder, "ct1.npy"))
        Ct2 = np.load(os.path.join(folder, "ct2.npy"))
        Ct1Ave = np.load(os.path.join(folder, "ct1Ave.npy"))
        Ct2Ave = np.load(os.path.join(folder, "ct2Ave.npy"))
        Ct1PPE = np.load(os.path.join(folder, "ct1PPE.npy"))
        Ct2PPE = np.load(os.path.join(folder, "ct2PPE.npy"))
        Ct1PPEAve = np.load(os.path.join(folder, "ct1PPEAve.npy"))
        Ct2PPEAve = np.load(os.path.join(folder, "ct2PPEAve.npy"))
        coeffU = np.load(os.path.join(folder, "coeefs.npy")) #100 modes x 5001 time steps
        coeffP = np.load(os.path.join(folder, "coeefsP.npy")) #50 modes x 5001 time steps
        coeffNut = np.load(os.path.join(folder, "coeefsNut.npy"))
        bcMat = read_txt(os.path.join(folder, "bcVelMat/bcVelMat0_mat.txt"))
        bcVec = read_txt(os.path.join(folder, "bcVelVec/bcVelVec0_mat.txt"))
        CtTot = Ct1 + Ct2
        CtPPETot = Ct1PPE + Ct2PPE
        CtAveTot = Ct1Ave + Ct2Ave
        CtPPEAveTot = Ct1PPEAve + Ct2PPEAve
        #nu 
        self.nu = 1e-4

        bcMat = bcMat.reshape((self.Nu_tot + self.Nsup_tot,self.Nu_tot + self.Nsup_tot))
        bcVec = bcVec.reshape((self.Nu_tot + self.Nsup_tot))
        self.CumEig_U = readFile(os.path.join(data_folder,"POD/CumEigenvalues_U"))
        self.CumEig_Usup = readFile(os.path.join(data_folder, "POD/CumEigenvalues_Usup"))
        self.CumEig_P = readFile(os.path.join(data_folder, "POD/CumEigenvalues_p"))
        Eig_U = readFile(os.path.join(data_folder, "POD/Eigenvalues_U"))
        coo = np.load("coords.npy")
        modU = read_vec_modes(os.path.join(data_folder, 'POD'), 11644, self.Nu_tot, "U") #34932x50
        modSup = read_vec_modes(os.path.join(data_folder, 'supremizer'), 11644, self.Nsup_tot, "Usup")  #34932x50
        modP = read_p_modes(os.path.join(data_folder, 'POD'), 11644, self.Np_tot) #11644x50
        ## Import useful data from FOM simulation
#        self.coo = coo[:,0:2]
        self.mod_Ux = modU[0:11644,:]
        self.mod_Uy = modU[11644:23288,:]
        self.mod_Uz = modU[23288:34932,:]
        self.modsup_Ux = modSup[0:11644,:]
        self.modsup_Uy = modSup[11644:23288,:]
        self.modsup_Uz = modSup[23288:34932,:]
        self.mod_P = modP
        self.coeffU = coeffU
        self.coeffP = coeffP
        self.EigU = Eig_U
        self.C = C

        ## Slice Matrices according to dimension
        self.B_ok = np.zeros([self.Nu+self.Nsup,self.Nu+self.Nsup])
        self.B_ok[0:self.Nu,0:self.Nu] = B[0:self.Nu,0:self.Nu]
        self.B_ok[self.Nu:self.Nu+self.Nsup,self.Nu:self.Nu+self.Nsup] = B[self.Nu_tot:self.Nu_tot+self.Nsup,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.B_ok[self.Nu:self.Nu+self.Nsup,0:self.Nu] = B[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nu]
        self.B_ok[0:self.Nu,self.Nu:self.Nu+self.Nsup] = B[0:self.Nu,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.Bt_ok = np.zeros([self.Nu+self.Nsup,self.Nu+self.Nsup])
        self.Bt_ok[0:self.Nu,0:self.Nu] = Bt[0:self.Nu,0:self.Nu]
        self.Bt_ok[self.Nu:self.Nu+self.Nsup,self.Nu:self.Nu+self.Nsup] = Bt[self.Nu_tot:self.Nu_tot+self.Nsup,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.Bt_ok[self.Nu:self.Nu+self.Nsup,0:self.Nu] = Bt[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nu]
        self.Bt_ok[0:self.Nu,self.Nu:self.Nu+self.Nsup] = Bt[0:self.Nu,self.Nu_tot:self.Nu_tot+self.Nsup]

        self.M_ok = np.zeros([self.Nu+self.Nsup,self.Nu+self.Nsup])
        self.M_ok[0:self.Nu,0:self.Nu] = M[0:self.Nu,0:self.Nu]
        self.M_ok[self.Nu:self.Nu+self.Nsup,self.Nu:self.Nu+self.Nsup] = M[self.Nu_tot:self.Nu_tot+self.Nsup,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.M_ok[self.Nu:self.Nu+self.Nsup,0:self.Nu] = M[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nu]
        self.M_ok[0:self.Nu,self.Nu:self.Nu+self.Nsup] = M[0:self.Nu,self.Nu_tot:self.Nu_tot+self.Nsup]

        self.K_ok = np.zeros([self.Nu+self.Nsup,self.Np])
        self.K_ok[0:self.Nu,0:self.Np] = K[0:self.Nu,0:self.Np]
        self.K_ok[self.Nu:self.Nu+self.Nsup,0:self.Np] = K[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Np]

        self.D_ok = np.zeros([self.Np,self.Np])
        self.D_ok = D[0:self.Np,0:self.Np]

        self.C_ok = np.zeros([self.Nu+self.Nsup,self.Nu+self.Nsup,self.Nu+self.Nsup])
        self.C_ok[0:self.Nu,0:self.Nu,0:self.Nu] = C[0:self.Nu,0:self.Nu,0:self.Nu]
        self.C_ok[self.Nu:self.Nu+self.Nsup,self.Nu:self.Nu+self.Nsup,self.Nu:self.Nu+self.Nsup] = C[self.Nu_tot:self.Nu_tot+self.Nsup,self.Nu_tot:self.Nu_tot+self.Nsup,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.C_ok[0:self.Nu,0:self.Nu,self.Nu:self.Nu+self.Nsup] = C[0:self.Nu,0:self.Nu,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.C_ok[0:self.Nu,self.Nu:self.Nu+self.Nsup,self.Nu:self.Nu+self.Nsup] = C[0:self.Nu,self.Nu_tot:self.Nu_tot+self.Nsup,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.C_ok[0:self.Nu,self.Nu:self.Nu+self.Nsup,0:self.Nu] = C[0:self.Nu,self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nu]
        self.C_ok[self.Nu:self.Nu+self.Nsup,0:self.Nu,0:self.Nu] = C[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nu,0:self.Nu]
        self.C_ok[self.Nu:self.Nu+self.Nsup,self.Nu:self.Nu+self.Nsup,0:self.Nu] = C[self.Nu_tot:self.Nu_tot+self.Nsup,self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nu]
        self.C_ok[self.Nu:self.Nu+self.Nsup,0:self.Nu,self.Nu:self.Nu+self.Nsup] = C[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nu,self.Nu_tot:self.Nu_tot+self.Nsup]

        self.G_ok = np.zeros([self.Np,self.Nu+self.Nsup,self.Nu+self.Nsup])
        self.G_ok[0:self.Np,0:self.Nu, 0:self.Nu] = G[0:self.Np,0:self.Nu, 0:self.Nu]
        self.G_ok[0:self.Np,self.Nu:self.Nu+self.Nsup, self.Nu:self.Nu+self.Nsup] = G[0:self.Np,self.Nu_tot:self.Nu_tot+self.Nsup, self.Nu_tot:self.Nu_tot+self.Nsup]

        self.CtTot_Ok = np.zeros([self.Nu+self.Nsup,self.Nnut,self.Nu+self.Nsup])
        self.CtTot_Ok[0:self.Nu,0:self.Nnut,0:self.Nu] = CtTot[0:self.Nu,0:self.Nnut,0:self.Nu]
        self.CtTot_Ok[self.Nu:self.Nu+self.Nsup,0:self.Nnut,self.Nu:self.Nu+self.Nsup] = CtTot[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nnut,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.CtTot_Ok[0:self.Nu,0:self.Nnut,self.Nu:self.Nu+self.Nsup] = CtTot[0:self.Nu,0:self.Nnut,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.CtTot_Ok[self.Nu:self.Nu+self.Nsup,0:self.Nnut,0:self.Nu] = CtTot[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nnut,0:self.Nu]
        self.CtTot_Ok[self.Nu:self.Nu+self.Nsup,0:self.Nnut,self.Nu:self.Nu+self.Nsup] = CtTot[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nnut,self.Nu_tot:self.Nu_tot+self.Nsup]

        self.CtPPETot_Ok = np.zeros([self.Nu+self.Nsup,self.Nnut,self.Nu+self.Nsup])
        self.CtPPETot_Ok[0:self.Nu,0:self.Nnut,0:self.Nu] = CtPPETot[0:self.Nu,0:self.Nnut,0:self.Nu]
        self.CtPPETot_Ok[self.Nu:self.Nu+self.Nsup,0:self.Nnut,self.Nu:self.Nu+self.Nsup] = CtPPETot[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nnut,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.CtPPETot_Ok[0:self.Nu,0:self.Nnut,self.Nu:self.Nu+self.Nsup] = CtPPETot[0:self.Nu,0:self.Nnut,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.CtPPETot_Ok[self.Nu:self.Nu+self.Nsup,0:self.Nnut,0:self.Nu] = CtPPETot[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nnut,0:self.Nu]
        self.CtPPETot_Ok[self.Nu:self.Nu+self.Nsup,0:self.Nnut,self.Nu:self.Nu+self.Nsup] = CtPPETot[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nnut,self.Nu_tot:self.Nu_tot+self.Nsup]


        self.P_ok = np.zeros([self.Np,self.Nu+self.Nsup])
        self.P_ok[0:self.Np,0:self.Nu] = P[0:self.Np,0:self.Nu]
        self.P_ok[0:self.Np,self.Nu:self.Nu+self.Nsup] = P[0:self.Np,self.Nu_tot:self.Nu_tot+self.Nsup]

        self.CtTot_Ok = np.zeros([self.Nu+self.Nsup,self.Nnut,self.Nu+self.Nsup])
        self.CtTot_Ok[0:self.Nu,0:self.Nnut,0:self.Nu] = CtTot[0:self.Nu,0:self.Nnut,0:self.Nu]
        self.CtTot_Ok[self.Nu:self.Nu+self.Nsup,0:self.Nnut,self.Nu:self.Nu+self.Nsup] = CtTot[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nnut,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.CtTot_Ok[0:self.Nu,0:self.Nnut,self.Nu:self.Nu+self.Nsup] = CtTot[0:self.Nu,0:self.Nnut,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.CtTot_Ok[self.Nu:self.Nu+self.Nsup,0:self.Nnut,0:self.Nu] = CtTot[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nnut,0:self.Nu]
        self.CtTot_Ok[self.Nu:self.Nu+self.Nsup,0:self.Nnut,self.Nu:self.Nu+self.Nsup] = CtTot[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nnut,self.Nu_tot:self.Nu_tot+self.Nsup]

        self.coeffU_ok = np.zeros([self.Nu+self.Nsup,coeffU.shape[1]])
        self.coeffNut_ok = coeffNut[0:self.Nnut,:]
        self.coeffP_ok = coeffP[0:self.Np,:]
        self.coeffU_ok[0:self.Nu,:] = coeffU[0:self.Nu,:]
        self.coeffU_ok[self.Nu:self.Nu+self.Nsup,:] = coeffU[self.Nu_tot:self.Nu_tot+self.Nsup,:]

        self.bcMat_ok = np.zeros([self.Nu+self.Nsup,self.Nu+self.Nsup])
        self.bcMat_ok[0:self.Nu,0:self.Nu] = bcMat[0:self.Nu,0:self.Nu]
        self.bcMat_ok[self.Nu:self.Nu+self.Nsup,self.Nu:self.Nu+self.Nsup] = bcMat[self.Nu_tot:self.Nu_tot+self.Nsup,self.Nu_tot:self.Nu_tot+self.Nsup]
        self.bcMat_ok[self.Nu:self.Nu+self.Nsup,0:self.Nu] = bcMat[self.Nu_tot:self.Nu_tot+self.Nsup,0:self.Nu]
        self.bcMat_ok[0:self.Nu,self.Nu:self.Nu+self.Nsup] = bcMat[0:self.Nu,self.Nu_tot:self.Nu_tot+self.Nsup]

        self.bcVec_ok = np.zeros([self.Nu+self.Nsup])
        self.bcVec_ok[0:self.Nu] = bcVec[0:self.Nu]
        self.bcVec_ok[self.Nu:self.Nu+self.Nsup] = bcVec[self.Nu_tot:self.Nu_tot+self.Nsup]

        #Times
        self.dtOff = 0.004
        self.ntOff = 5000
        self.timeOffline = np.linspace(20,20+self.ntOff*self.dtOff,self.ntOff+1)

    def correction(self):
        '''
        Function that calculates the correction terms starting from data
        of the first self.M1 snapshots.
        '''
        term1 = np.zeros((self.Nu_tot,self.M1))
        term2 = np.zeros((self.Nu,self.M1))
        self.tau_ex_U = np.zeros((self.Nu,self.M1))
        Xhat_Ut = np.zeros((self.Nu,self.M1))

        for j in range(self.M1):
            term1[:,j] = - np.transpose(coeffU[0:self.Nu_tot,j])@C[0:self.Nu_tot,0:self.Nu_tot,0:self.Nu_tot]@(coeffU[0:self.Nu_tot,j])
            term2[:,j] = - np.transpose(self.coeffU_ok[0:self.Nu,j])@self.C_ok[0:self.Nu,0:self.Nu,0:self.Nu]@(self.coeffU_ok[0:self.Nu,j])
            self.tau_ex_U[:,j] = term1[0:self.Nu,j]-term2[:,j]
            Xhat_Ut[:,j] = self.coeffU_ok[0:self.Nu,j] #rxM--> Xhat is Mxr

        self.tau_ex_U = self.tau_ex_U.T
        Xhat_U = Xhat_Ut.T


        ## Definition of the term tau exact for pressure (in second equation) for term Db and a.T@G@a
        ##model the sum of these two terms with Dtilde_pg@ab+ab.T@Gtilde_pg@a
        term1_pg = np.zeros((self.Np_tot,self.M1))
        term2_pg = np.zeros((self.Np,self.M1))
        self.tau_ex_pg = np.zeros((self.Np,self.M1))
        Xhat_t = np.zeros((self.Np,self.M1))

        for j in range(self.M1):
            term1_pg[:,j] = D[0:self.Np_tot,0:self.Np_tot]@(coeffP[0:self.Np_tot,j]) + np.transpose(coeffU[0:self.Nu_tot,j])@G[0:self.Np_tot,0:self.Nu_tot,0:self.Nu_tot]@(coeffU[0:self.Nu_tot,j])
            term2_pg[:,j] = self.D_ok[0:self.Np,0:self.Np]@(self.coeffP_ok[0:self.Np,j]) + np.transpose(self.coeffU_ok[0:self.Nu,j])@self.G_ok[0:self.Np,0:self.Nu,0:self.Nu]@(self.coeffU_ok[0:self.Nu,j])
            self.tau_ex_pg[:,j] = term1_pg[0:self.Np,j]-term2_pg[:,j]
            Xhat_t[:,j] = self.coeffP_ok[0:self.Np,j]
        self.tau_ex_pg = self.tau_ex_pg.T
        Xhat = Xhat_t.T
        ## Construction of the least square problem for pressure with linear term and quadratic term w.r.t. ab 

        ab_hat = np.concatenate((Xhat_U, Xhat), axis=1)
        Dmat_ab = np.zeros((self.M1,self.Np+self.Nu))
        Dmat_ab = ab_hat  #initialization of matrix D

        for i in range(1,self.Nu+self.Np+1,1):
            ab_hati = np.zeros((self.M1,i))
            for j in range(self.M1):
                ab_hati[j,0:i] = np.dot(ab_hat[j,i-1], ab_hat[j,0:i])
            Dmat_ab = np.append(Dmat_ab,ab_hati,axis=1)

        elements4 = np.linspace(1,self.Nu+self.Np,self.Nu+self.Np)
        rows4 = sum(elements4)+ self.Np + self.Nu
        self.rows4 = int(rows4)


        ##Definition of a exact corrective term for both first and second equation
        self.tau_ex_tot = np.zeros((self.M1,self.Nu + self.Np))
        self.tau_ex_tot = np.concatenate((self.tau_ex_U,self.tau_ex_pg), axis=1)

        ##-----------Constuction of a unique global least square problem for all corrective terms-------------------##
        O_tot = np.zeros((self.rows4,self.Np+self.Nu))

        ## Truncated SVD applied to matrix D
        U_tot, Sigma_tot, VT_tot = randomized_svd(Dmat_ab, n_components=self.R_tot, random_state=None)
        Sigma_tot = np.diag(Sigma_tot)
        for i in range(self.Np+self.Nu):
            o_itot = (VT_tot.T@(np.linalg.inv(Sigma_tot))@U_tot.T)@self.tau_ex_tot[:,i]
            O_tot[:,i] = o_itot

        self.Omat_tot = O_tot.T
        self.Jtilde_A = self.Omat_tot[:,0:self.Np+self.Nu]
        self.Jtilde_B = self.Omat_tot[:,self.Np+self.Nu:self.rows4]
        return self.Jtilde_A, self.Jtilde_B

    def U_50modes(self):
        '''
        Function which computes the full order velocity field retaining 50 modes,
        at each time step and for each node of the grid.
        '''
        U_x50 = np.zeros((11644,self.M))
        U_y50 = np.zeros((11644,self.M))
        U_z50 = np.zeros((11644,self.M))
        self.U_norm_50 = np.zeros((11644,self.M))
        for i in range(self.M):
            for j in range(self.Nu_tot):
                U_x50[:,i] = U_x50[:,i]+self.coeffU[j,i]*self.mod_Ux[:,j]
                U_y50[:,i] = U_y50[:,i]+self.coeffU[j,i]*self.mod_Uy[:,j]
                U_z50[:,i] = U_z50[:,i]+self.coeffU[j,i]*self.mod_Uz[:,j]
            for m in range(self.Nsup_tot):
                U_x50[:,i] = U_x50[:,i]+self.coeffU[self.Nu_tot+m,i]*self.modsup_Ux[:,m]
                U_y50[:,i] = U_y50[:,i]+self.coeffU[self.Nu_tot+m,i]*self.modsup_Uy[:,m]
                U_z50[:,i] = U_z50[:,i]+self.coeffU[self.Nu_tot+m,i]*self.modsup_Uz[:,m]
            self.U_norm_50[:,i] = np.sqrt(U_x50[:,i]**2+U_y50[:,i]**2+U_z50[:,i]**2)
        return self.U_norm_50

    def U_10modes(self):
        '''
        Function which computes the full order velocity field retaining 10 modes,
        at each time step and for each node of the grid.
        '''
        U_x10 = np.zeros((11644,self.M))
        U_y10 = np.zeros((11644,self.M))
        U_z10 = np.zeros((11644,self.M))
        self.U_norm_10 = np.zeros((11644,self.M))
        for i in range(self.M):
            for j in range(self.Nu):
                U_x10[:,i] = U_x10[:,i]+self.coeffU[j,i]*self.mod_Ux[:,j]
                U_y10[:,i] = U_y10[:,i]+self.coeffU[j,i]*self.mod_Uy[:,j]
                U_z10[:,i] = U_z10[:,i]+self.coeffU[j,i]*self.mod_Uz[:,j]
            for m in range(self.Nsup):
                U_x10[:,i] = U_x10[:,i]+self.coeffU[self.Nu_tot+m,i]*self.modsup_Ux[:,m]
                U_y10[:,i] = U_y10[:,i]+self.coeffU[self.Nu_tot+m,i]*self.modsup_Uy[:,m]
                U_z10[:,i] = U_z10[:,i]+self.coeffU[self.Nu_tot+m,i]*self.modsup_Uz[:,m]
            self.U_norm_10[:,i] = np.sqrt(U_x10[:,i]**2+U_y10[:,i]**2+U_z10[:,i]**2)
        return self.U_norm_10


    def P_50modes(self):
        '''
        Function which computes the full order pressure field retaining 50 modes,
        at each time step and for each node of the grid.
        '''
        self.P_50 = np.zeros((11644,self.M))
        for i in range(self.M):
            for j in range(self.Np_tot):
                self.P_50[:,i] = self.P_50[:,i]+self.coeffP[j,i]*self.mod_P[:,j]
        return self.P_50

    def P_10modes(self):
        '''
        Function which computes the full order pressure field retaining 10 modes,
        at each time step and for each node of the grid.
        '''
        self.P_10 = np.zeros((11644,self.M))
        for i in range(self.M):
            for j in range(self.Np):
                self.P_10[:,i] = self.P_10[:,i]+self.coeffP[j,i]*self.mod_P[:,j]
        return self.P_10


    def solveOnline(self, _dtOn, _ntOn):
        self.dtOn = _dtOn
        self.ntOn = _ntOn
        self.timeOnline = np.linspace(20,20+self.ntOn*self.dtOn,self.ntOn)
        self.a_0 = self.coeffU_ok[:,0]
        self.a = self.a_0
        self.b_0 = self.coeffP_ok[:,0]
        self.ab_0 = np.concatenate((self.a_0,self.b_0))
        self.a_old = self.a_0
        self.b_old = self.b_0

        ab = self.ab_0
        a = self.a_0
        self.i = 0
        for i,t in enumerate(self.timeOnline):
            self.i = i

            if i == 0:
                self.sol = ab.reshape(-1,1)
                ab = fsolve(self.residual0, ab)
            if i == 1:
                self.sol = ab.reshape(-1,1)
                ab = fsolve(self.residual0, ab)

            else:
                self.sol = np.hstack((self.sol,ab.reshape(-1,1)))
                ab = fsolve(self.residual, ab)
            self.a_old = ab[0:self.Nu+self.Nsup]
            self.b_old = ab[self.Nu+self.Nsup:self.Np+self.Nu+self.Nsup]

    # Solve the standard online problem
    def solveOnline_standard(self, _dtOn, _ntOn):
        self.dtOn = _dtOn
        self.ntOn = _ntOn
        self.timeOnline = np.linspace(20,20+self.ntOn*self.dtOn,self.ntOn)
        self.a_0 = self.coeffU_ok[:,0]
        self.a = self.a_0
        self.b_0 = self.coeffP_ok[:,0]
        self.ab_0 = np.concatenate((self.a_0,self.b_0))
        self.a_old = self.a_0
        self.b_old = self.b_0

        ab = self.ab_0
        a = self.a_0
        self.i = 0
        for i,t in enumerate(self.timeOnline):
            self.i = i
            if i == 0:
                self.sol = ab.reshape(-1,1)
                ab = fsolve(self.residual0_stand, ab)
            if i == 1:
                self.sol = ab.reshape(-1,1)
                ab = fsolve(self.residual0_stand, ab)
            else:
                self.sol = np.hstack((self.sol,ab.reshape(-1,1)))
                ab = fsolve(self.residual_stand, ab)
            self.a_old = ab[0:self.Nu+self.Nsup]
            self.b_old = ab[self.Nu+self.Nsup:self.Np+self.Nu+self.Nsup]


    # Solve the online problem with turbulence modelling and including the velocity correction
    def solveOnlineT(self, _dtOn, _ntOn):
        self.dtOn = _dtOn
        self.ntOn = _ntOn
        self.timeOnline = np.linspace(20,20+self.ntOn*self.dtOn,self.ntOn)
        self.a_0 = self.coeffU_ok[:,0]
        self.a = self.a_0
        self.b_0 = self.coeffP_ok[:,0]
        self.g_0 = self.coeffNut_ok[:,0]
        self.ab_0 = np.concatenate((self.a_0,self.b_0))
        self.a_old = self.a_0
        self.b_old = self.b_0

        ab = self.ab_0
        a = self.a_0
        self.i = 0
        for i,t in enumerate(self.timeOnline):
            self.i = i
            if i == 0:
                self.sol = ab.reshape(-1,1)
               # self.gsol = self.g_w(ab).reshape(-1,1)
                ab = fsolve(self.residualT0, ab)
            if i == 1:
                self.sol = ab.reshape(-1,1)
               # self.gsol = self.g_w(ab).reshape(-1,1)
                ab = fsolve(self.residualT0, ab)
            else:
                self.sol = np.hstack((self.sol,ab.reshape(-1,1)))
               # self.gsol = np.hstack((self.gsol,self.g_w(ab).reshape(-1,1)))
                ab = fsolve(self.residualT, ab)
            self.a_old = ab[0:self.Nu+self.Nsup]
            self.b_old = ab[self.Nu+self.Nsup:self.Np+self.Nu+self.Nsup]


    # Residual without turbulence modelling and with correction for initial time steps
    def residual0(self,ab):
        J_A, J_B = self.correction()
        res0 = ab*0
        a = ab[0:self.Nu+self.Nsup]
        a_nu = a[0:self.Nu]
        b = ab[self.Nu+self.Nsup:self.Nu+self.Nsup+self.Np]
        adot0 = (a-self.a_old)/self.dtOn
        Jterm_A = np.zeros(self.Np+self.Nu)
        Jterm_A[:] = J_A@ab

        Jterm_B = np.zeros(self.Np+self.Nu)
        k = 0
        for i in range(1,self.Nu+self.Np+1,1):
            ab_i = np.zeros(i)
            ab_i[0:i] = np.dot(ab[i-1],ab[0:i])
            Jtilde_Bi = J_B[:,k:k+i]
            k = k+i
            Jterm_B[:] = Jterm_B[:]+Jtilde_Bi@ab_i

        pen = self.BC * self.bcVec_ok - self.bcMat_ok@a
        res0[0:self.Nu+self.Nsup] = -self.M_ok@adot + self.nu*(self.B_ok+self.Bt_ok)@a - a.T@self.C_ok@a - self.K_ok@b  + pen*self.tau + Jterm_A[0:self.Nu] + Jterm_B[0:self.Nu]
        res0[-self.Np:] = self.D_ok@b + a.T@self.G_ok@a -self.nu*self.N_ok@a + Jterm_A[self.Nu:self.Nu+self.Np] + Jterm_B[self.Nu:self.Nu+self.Np]
        return res0

    # Residual without turbulence modelling and with correction 
    def residual(self,ab):
        J_A, J_B = self.correction()
        res0 = ab*0
        a = ab[0:self.Nu+self.Nsup]
        a_nu = a[0:self.Nu]
        b = ab[self.Nu+self.Nsup:self.Nu+self.Nsup+self.Np]
        adot0 = (3*a - 4*self.a_old + self.sol[0:self.Nu+self.Nsup,self.i-2])/(2*self.dtOn)
        Jterm_A = np.zeros(self.Np+self.Nu)
        Jterm_A[:] = J_A@ab

        Jterm_B = np.zeros(self.Np+self.Nu)
        k = 0
        for i in range(1,self.Nu+self.Np+1,1):
            ab_i = np.zeros(i)
            ab_i[0:i] = np.dot(ab[i-1],ab[0:i])
            Jtilde_Bi = J_B[:,k:k+i]
            k = k+i
            Jterm_B[:] = Jterm_B[:]+Jtilde_Bi@ab_i

        pen = self.BC * self.bcVec_ok - self.bcMat_ok@a
        res0[0:self.Nu+self.Nsup] = -self.M_ok@adot + self.nu*(self.B_ok+self.Bt_ok)@a - a.T@self.C_ok@a - self.K_ok@b  + pen*self.tau + Jterm_A[0:self.Nu] + Jterm_B[0:self.Nu]
        res0[-self.Np:] = self.D_ok@b + a.T@self.G_ok@a -self.nu*self.N_ok@a + Jterm_A[self.Nu:self.Nu+self.Np] + Jterm_B[self.Nu:self.Nu+self.Np]
        return res0

    # Residual with turbulence modelling and correction for initial time steps
    def residualT0(self,ab):
        J_A, J_B = self.correction()
        res = ab*0
        a = ab[0:self.Nu+self.Nsup]
        a_nu = a[0:self.Nu]
        b = ab[self.Nu+self.Nsup:self.Nu+self.Nsup+self.Np]
        adot = (a - self.a_old)/self.dtOn
        g = self.g_w(ab)
        Jterm_A = np.zeros(self.Np+self.Nu)
        Jterm_A[:] = J_A@ab

        Jterm_B = np.zeros(self.Np+self.Nu)
        k = 0
        for i in range(1,self.Nu+self.Np+1,1):
            ab_i = np.zeros(i)
            ab_i[0:i] = np.dot(ab[i-1],ab[0:i])
            Jtilde_Bi = J_B[:,k:k+i]
            k = k+i
            Jterm_B[:] = Jterm_B[:]+Jtilde_Bi@ab_i

        pen = self.BC * self.bcVec_ok - self.bcMat_ok@a
        res[0:self.Nu+self.Nsup] = -self.M_ok@adot + self.nu*(self.B_ok+self.Bt_ok)@a - a.T@self.C_ok@a - self.K_ok@b + g.T@(self.CtTot_Ok)@a + pen*self.tau + Jterm_A[0:self.Nu] + Jterm_B[0:self.Nu]
        res[-self.Np:] = self.D_ok@b + a.T@self.G_ok@a -g.T@self.CtPPETot_Ok@a -self.nu*self.N_ok@a + Jterm_A[self.Nu:self.Nu+self.Np] + Jterm_B[self.Nu:self.Nu+self.Np]
        return res

    # Residual with turbulence modelling and correction 
    def residualT(self,ab):
        J_A, J_B = self.correction()
        res = ab*0
        a = ab[0:self.Nu+self.Nsup]
        a_nu = a[0:self.Nu]
        b = ab[self.Nu+self.Nsup:self.Nu+self.Nsup+self.Np]
        adot = (3*a - 4*self.a_old + self.sol[0:self.Nu+self.Nsup,self.i-2])/(2*self.dtOn)
        g = self.g_w(ab)
        Jterm_A = np.zeros(self.Np+self.Nu)
        Jterm_A[:] = J_A@ab

        Jterm_B = np.zeros(self.Np+self.Nu)
        k = 0
        for i in range(1,self.Nu+self.Np+1,1):
            ab_i = np.zeros(i)
            ab_i[0:i] = np.dot(ab[i-1],ab[0:i])
            Jtilde_Bi = J_B[:,k:k+i]
            k = k+i
            Jterm_B[:] = Jterm_B[:]+Jtilde_Bi@ab_i

        pen = self.BC * self.bcVec_ok - self.bcMat_ok@a
        res[0:self.Nu+self.Nsup] = -self.M_ok@adot + self.nu*(self.B_ok+self.Bt_ok)@a - a.T@self.C_ok@a - self.K_ok@b + g.T@(self.CtTot_Ok)@a + pen*self.tau + Jterm_A[0:self.Nu] + Jterm_B[0:self.Nu]
        res[-self.Np:] = self.D_ok@b + a.T@self.G_ok@a -g.T@self.CtPPETot_Ok@a -self.nu*self.N_ok@a + Jterm_A[self.Nu:self.Nu+self.Np] + Jterm_B[self.Nu:self.Nu+self.Np]
        return res

    # Residual without turbulence modelling and correction for initial time steps
    def residual0_stand(self,ab):
        res = ab*0
        a = ab[0:self.Nu+self.Nsup]
        a_nu = a[0:self.Nu]
        ag = ab[0:self.Nu]
        b = ab[self.Nu+self.Nsup:self.Nu+self.Nsup+self.Np]
        adot = (a - self.a_old)/self.dtOn
        pen = self.BC * self.bcVec_ok - self.bcMat_ok@a
        res[0:Nu+Nsup] = -self.M_ok@adot + self.nu*(self.B_ok+self.Bt_ok)@a - a.T@self.C_ok@a - self.K_ok@b + pen*self.tau
        res[-Np:] = self.D_ok@b + a.T@self.G_ok@a -self.nu*self.N_ok@a
        return res

    # Residual without turbulence modelling and correction 
    def residual_stand(self,ab):
        res = ab*0
        a = ab[0:self.Nu+self.Nsup]
        a_nu = a[0:self.Nu]
        ag = ab[0:self.Nu]
        b = ab[self.Nu+self.Nsup:self.Nu+self.Nsup+self.Np]
        adot = (3*a - 4*self.a_old + self.sol[0:self.Nu+self.Nsup,self.i-2])/(2*self.dtOn)
        pen = self.BC * self.bcVec_ok - self.bcMat_ok@a
        res[0:Nu+Nsup] = -self.M_ok@adot + self.nu*(self.B_ok+self.Bt_ok)@a - a.T@self.C_ok@a - self.K_ok@b + pen*self.tau
        res[-Np:] = self.D_ok@b + a.T@self.G_ok@a -self.nu*self.N_ok@a
        return res
    # Compute coefficients the nut POD expansion by training a network with a weighted loss function
    def g_w(self,ab):
        ag = ab[0:self.Nu]
        atr = self.scaling["scaler_inp"].transform(ag.reshape(1, -1))[0,:]
        at = torch.from_numpy(atr).float()
        gt = self.Net_w.forward(at)
        g = gt.detach().numpy()
        g = self.scaling["scaler_out"].inverse_transform(g.reshape(1, -1))[0,:]
        return g

    def weighted_mse_loss(self, inp, target, weight):
        weight = weight/sum(weight)
        loss = torch.sum(weight * (inp - target) ** 2)
        return loss

    # Train the NET for a -> g mapping with weighted loss
    def train_w(self,layers, epochs, batch_size, learning_rate, act):
        # Input data for training
        self.X_train, self.X_test, self.Y_train, self.Y_test = train_test_split(self.coeffU_ok[0:self.Nu,:].T, self.coeffNut_ok.T, test_size=0.33, random_state=42)
        self.scaling = {"scaler_inp": preprocessing.MinMaxScaler(),
           "scaler_out": preprocessing.MinMaxScaler()}
        self.inp_np_train = self.scaling["scaler_inp"].fit_transform(self.X_train)
        self.out_np_train = self.scaling["scaler_out"].fit_transform(self.Y_train)
        self.inp_np_test = self.scaling["scaler_inp"].transform(self.X_test)
        self.out_np_test = self.scaling["scaler_out"].transform(self.Y_test)

        self.trainset = TrainNet.NutDataset(self.inp_np_train, self.out_np_train)
        self.testset = TrainNet.NutDataset(self.inp_np_test, self.out_np_test)
        self.trainloader = torch.utils.data.DataLoader(self.trainset, batch_size=batch_size, shuffle=True)
        self.testloader = torch.utils.data.DataLoader(self.testset, batch_size=self.inp_np_test.shape[0], shuffle=True)

        # Create Net
        self.Net_w = TrainNet.NetVelocityNut(self.Nu,self.Nnut,layers,act)
        # Define Loss
        self.loss_fn = torch.nn.MSELoss()
        self.weights = self.EigU[0:self.Nu]
        self.weights = torch.tensor(self.weights, dtype=torch.float32)
        self.modelfile = self.Net_w.hash()+"_w.pt"
        optimizer = torch.optim.Adam(self.Net_w.parameters(), lr=learning_rate)
        if os.path.isfile(self.modelfile):
            self.Net_w = torch.jit.load(self.modelfile)
        else:
            for t in range(epochs):
                batch_losses = []
                for inputs, labels in self.trainloader:
                    inputs, labels =  inputs.to(self.device), labels.to(self.device)
                    optimizer.zero_grad()
                    # forward + backward + optimize
                    outputs = self.Net_w(inputs)
                    loss = self.weighted_mse_loss(outputs, labels, self.weights)
                    loss.backward()
                    optimizer.step()
                    batch_losses.append(loss.item())
                loss = np.mean(batch_losses)

                # evaluate accuracy on test set
                batch_test_losses = []
                self.Net_w.eval()
                for inputs_test, labels_test in self.testloader:
                    inputs_test, labels_test =  inputs_test.to(self.device), labels_test.to(self.device)
                    outputs_test = self.Net_w(inputs_test)
                    test_loss = self.weighted_mse_loss(outputs_test, labels_test, self.weights)
                    batch_test_losses.append(test_loss.item())
                test_loss = np.mean(batch_test_losses)
                if t % 100 == 99:
                    print(t, "loss on train" , loss)
                    print(t, "loss on test" , test_loss)
                    self.tplot.append(t)
                    self.lossplottrain.append(loss)
                    self.lossplottest.append(test_loss)
                self.Net_w.train()

            # Save the model and scaling
            m = torch.jit.script(self.Net_w)
            np.save(self.Net_w.hash()+"_scal_inp_min_w.npy",self.scaling["scaler_inp"].min_)
            np.save(self.Net_w.hash()+"_scal_inp_sca_w.npy",self.scaling["scaler_inp"].scale_)
            np.save(self.Net_w.hash()+"_scal_out_min_w.npy",self.scaling["scaler_out"].min_)
            np.save(self.Net_w.hash()+"_scal_out_sca_w.npy",self.scaling["scaler_out"].scale_)
            np.save(self.Net_w.hash()+"_trainLoss_w.npy",self.lossplottrain)
            np.save(self.Net_w.hash()+"_testLoss_w.npy",self.lossplottest)
            m.save(self.Net_w.hash()+"_w.pt")

## Plot errors of the solution w.r.t. projection and w.r.t. full order solution
    def errors(self):
        U_50 = self.U_50modes()
        U_10 = self.U_10modes()
        U_sol = np.zeros((11644,self.M))
        U_solx = np.zeros((11644,self.M))
        U_soly = np.zeros((11644,self.M))
        U_solz = np.zeros((11644,self.M))
        P_sol = np.zeros((11644,self.M))
        P_10 = self.P_10modes()
        P_50 = self.P_50modes()
        for j in range(self.M-1):
            for i in range(self.Nu):
                U_solx[:,j] = U_solx[:,j]+self.sol[i,j]*self.mod_Ux[:,i]
                U_soly[:,j] = U_soly[:,j]+self.sol[i,j]*self.mod_Uy[:,i]
                U_solz[:,j] = U_solz[:,j]+self.sol[i,j]*self.mod_Uz[:,i]
            for m in range(self.Nsup):
                U_solx[:,j] = U_solx[:,j]+self.sol[m+self.Nu,j]*self.modsup_Ux[:,m]
                U_soly[:,j] = U_soly[:,j]+self.sol[m+self.Nu,j]*self.modsup_Uy[:,m]
                U_solz[:,j] = U_solz[:,j]+self.sol[m+self.Nu,j]*self.modsup_Uz[:,m]
            for n in range(self.Np):
                P_sol[:,j] = P_sol[:,j]+self.sol[self.Nu+self.Nsup+n,j]*self.mod_P[:,n]
            U_sol[:,j] = np.sqrt(U_solx[:,j]**2+U_soly[:,j]**2+U_solz[:,j]**2)
            self.errU_proj[j] = np.linalg.norm(U_sol[:,j]-U_10modes[:,j])/(np.linalg.norm(U_10[:,j]))*100
            self.errU_full[j] = np.linalg.norm(U_sol[:,j]-U_50modes[:,j])/(np.linalg.norm(U_50[:,j]))*100
            self.errP_proj[j] = np.linalg.norm(P_sol[:,j]-P_10[:,j])/(np.linalg.norm(P_10[:,j]))*100
            self.errP_full[j] = np.linalg.norm(P_sol[:,j]-P_50[:,j])/(np.linalg.norm(P_50[:,j]))*100

    def Ufield(self):
        U_solx = np.zeros((11644,self.M))
        U_soly = np.zeros((11644,self.M))
        U_solz = np.zeros((11644,self.M))
        U_sol = np.zeros((34932,self.M-1))
        for j in range(self.M-1):
            for i in range(self.Nu):
                U_solx[:,j] = U_solx[:,j]+self.sol[i,j]*self.mod_Ux[:,i]
                U_soly[:,j] = U_soly[:,j]+self.sol[i,j]*self.mod_Uy[:,i]
                U_solz[:,j] = U_solz[:,j]+self.sol[i,j]*self.mod_Uz[:,i]
            for n in range(self.Nsup):
                U_solx[:,j] = U_solx[:,j]+self.sol[self.Nu+n,j]*self.modsup_Ux[:,n]
                U_soly[:,j] = U_soly[:,j]+self.sol[self.Nu+n,j]*self.modsup_Uy[:,n]
                U_solz[:,j] = U_solz[:,j]+self.sol[self.Nu+n,j]*self.modsup_Uz[:,n]
        U_sol[0:11644,:] = U_solx[:,0:2000]
        U_sol[11644:23288,:] = U_soly[:,0:2000]
        U_sol[23288:34932,:] = U_solz[:,0:2000]
        return U_sol
    
    def Pfield(self):
        P_sol = np.zeros((11644,self.M-1))
        for j in range(self.M-1):
            for n in range(self.Np):
                P_sol[:,j] = P_sol[:,j]+self.sol[self.Nu+self.Nsup+n,j]*self.mod_P[:,n]
        return P_sol
    
    ## Define the variable epsilon_U, a measure of the error of the reduced velocity solution w.r.t. the projection
    def epsilon_U(self):
        Ufield = self.Ufield
        errU = np.zeros((self.M, 1))
        for j in range(self.M-1):
            U_sol[:,j] = np.sqrt(U_sol[0:11644,j]**2+U_sol[11644:23288,j]**2+U_sol[23288:34932,j]**2)
            errU[j] = np.linalg.norm(U_sol[:,j]-self.U_norm_10[:,j])
        eps_U = sum(errU)
        return eps_U

    ## Define the percentage error of the velocity solution w.r.t. the full order solution
    def U_proj(self):
        U_x10 = np.zeros((11644,self.M))
        U_y10 = np.zeros((11644,self.M))
        U_z10 = np.zeros((11644,self.M))
        U_norm_proj = np.zeros((11644,self.M))
        errU_p = np.zeros(self.M)
        U_50 = self.U_50modes()
        for i in range(self.M):
            for j in range(self.Nu):
                U_x10[:,i] = U_x10[:,i]+self.coeffU[j,i]*self.mod_Ux[:,j]
                U_y10[:,i] = U_y10[:,i]+self.coeffU[j,i]*self.mod_Uy[:,j]
                U_z10[:,i] = U_z10[:,i]+self.coeffU[j,i]*self.mod_Uz[:,j]
            for m in range(self.Nsup):
                U_x10[:,i]=U_x10[:,i]+self.coeffU[self.Nu_tot+m,i]*self.modsup_Ux[:,m]
                U_y10[:,i]=U_y10[:,i]+self.coeffU[self.Nu_tot+m,i]*self.modsup_Uy[:,m]
                U_z10[:,i]=U_z10[:,i]+self.coeffU[self.Nu_tot+m,i]*self.modsup_Uz[:,m]
            U_norm_proj[:,i] = np.sqrt(U_x10[:,i]**2+U_y10[:,i]**2+U_z10[:,i]**2)
            errU_p[i] = np.linalg.norm(U_norm_proj[:,i]-U_50[:,i])/(np.linalg.norm(U_50[:,i]))*100
        return errU_p

    def P_proj(self):
        errP_p=np.zeros(self.M)
        P_proj=np.zeros((11644,self.M))
        P_50 = self.P_50modes()
        for i in range(self.M):
            for j in range(self.Np):
                P_proj[:,i]=P_proj[:,i]+self.coeffP[j,i]*self.mod_P[:,j]
            errP_p[i] = np.linalg.norm(P_proj[:,i]-P_50[:,i])/(np.linalg.norm(P_50[:,i]))*100
        return errP_p
