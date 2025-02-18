import numpy as np
from files import *
import torch
from torch import nn
from torch.autograd import Variable
import matplotlib.pyplot as plt
import os
from sklearn.preprocessing import MinMaxScaler
from sklearn import preprocessing
from sklearn.metrics import r2_score
from torch.utils.data import Dataset
from scipy.io import mmread
import io

torch.manual_seed(0)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(0)

device = torch.device('cpu')
device
batch_size = 500

class NutDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.from_numpy(X).type(torch.float32)
        self.y = torch.from_numpy(y).type(torch.float32)
 
    def __len__(self):
        return len(self.X)
 
    def __getitem__(self, idx):
        return [self.X[idx], self.y[idx]]

class Net:
    def __init__(self,Nu,Nnut,Epochs = 10000):
        self.Nu = Nu
        self.Nnut = Nnut
        self.Epochs = Epochs

    def read(self):
        NOffSnap = np.loadtxt("ITHACAoutput/Offline/snaps",dtype=int)
        NOnSnap = np.loadtxt("ITHACAoutput/checkOff/snaps",dtype=int)
        
        ## Read the coefficients train
        # U
        inp_np_train_U = np.load("ITHACAoutput/NN/coeffs/coeffL2UTrain.npy")
        # P
        inp_np_train_P = np.load("ITHACAoutput/NN/coeffs/coeffL2PTrain.npy")
        # Nut
        out_np_train = np.load("ITHACAoutput/NN/coeffs/coeffL2NutTrain.npy")
        # Read Parameters from file train
        pars_train = np.loadtxt("parsOff_mat.txt")
        #NOffSnap = np.load("train/NOffSnap.npy")
        pars_train_np = []
        # Fill the train Parameters
        for k,j in enumerate(NOffSnap): # k is the index of NOffSnaps, j is the value
            for i in range(j):
                pars_train_np.append(pars_train[k])
        pars_train_np = np.asarray(pars_train_np) # convert into an array

        # self.weights = Variable(torch.from_numpy(np.transpose(mmread("ITHACAoutput/POD/Eigenvalues_nut")[0:self.Nnut])).type(torch.float32),requires_grad=True)        
        self.weights = Variable(torch.from_numpy((mmread("ITHACAoutput/POD/Eigenvalues_nut")[0:self.Nnut])).type(torch.float32),requires_grad=True)        
        # Read the coefficients test
        # U
        inp_np_test_U = np.load("ITHACAoutput/NN/coeffs/coeffL2UTest.npy")
        # P
        inp_np_test_P = np.load("ITHACAoutput/NN/coeffs/coeffL2PTest.npy")
        # Nut
        out_np_test = np.load("ITHACAoutput/NN/coeffs/coeffL2NutTest.npy")
        # Read Parameters from file test
        pars_test = np.loadtxt("parsOn_mat.txt")

        #NOnSnap = np.load("test/NOnSnap.npy")
        pars_test_np = []
        # Fill the train Parameters
        for k,j in enumerate(NOnSnap):
            for i in range(j):
                pars_test_np.append(pars_test[k])
        pars_test_np = np.asarray(pars_test_np)
        
        # Prepare dataset with and without angles
        #self.inp_np_train_a = np.append(np.transpose(np.expand_dims(pars_train_np,axis=0)),inp_np_train_U[:,0:self.Nu], axis = 1)
        #self.inp_np_test_a = np.append(np.transpose(np.expand_dims(pars_test_np,axis=0)),inp_np_test_U[:,0:self.Nu], axis = 1)
        print(pars_train_np.shape)
        print(inp_np_train_U.shape)
        self.inp_np_train_a = np.append(pars_train_np,inp_np_train_U[:,0:self.Nu], axis = 1)
        print(self.inp_np_train_a.shape)
        self.inp_np_test_a = np.append(pars_test_np,inp_np_test_U[:,0:self.Nu], axis = 1)
        print(self.inp_np_test_a.shape)
        self.inp_np_train_noa = inp_np_train_U[:,0:self.Nu]
        self.inp_np_test_noa = inp_np_test_U[:,0:self.Nu]
        self.out_np_train = out_np_train[:,0:self.Nnut]
        self.out_np_test = out_np_test[:,0:self.Nnut]

    def weighted_loss(self, x, labels):
        loss = (x - labels)**2
        loss = torch.mm(loss,self.weights)
        # loss = torch.square(loss)
        # pct_var = (x - labels)**2
        # out = pct_var * self.weights.expand_as(labels)
        loss = loss.sum()
        return loss


    def Net1(self,inp_np_train, inp_np_test, out_np_train, out_np_test, epochs, learning_rate = 0.001, weight_decay = 1e-7, batch_size = 500):
        # Create NN Network 1
        Nin = inp_np_train.shape[1]
        Nout = out_np_train.shape[1]
        model = torch.nn.Sequential(
            #torch.nn.Linear(Nin, 64),
            #torch.nn.Tanh(),
            #torch.nn.Linear(64, 256),
            #torch.nn.Tanh(),
	    #torch.nn.Linear(256, 64),
            #torch.nn.Tanh(),
            torch.nn.Linear(Nin, 256),
            torch.nn.Tanh(),
            torch.nn.Linear(256, 64),
            torch.nn.Tanh(),
            torch.nn.Linear(64, Nout),
        )
        
        scaling = {"scaler_inp": preprocessing.MinMaxScaler(),
                   "scaler_out": preprocessing.MinMaxScaler()}
        inp_np_train = scaling["scaler_inp"].fit_transform(inp_np_train)
        out_np_train = scaling["scaler_out"].fit_transform(out_np_train)
        inp_np_test = scaling["scaler_inp"].transform(inp_np_test)
        out_np_test = scaling["scaler_out"].transform(out_np_test)
        trainset = NutDataset(inp_np_train, out_np_train)
        testset = NutDataset(inp_np_test, out_np_test)
        
        trainloader = torch.utils.data.DataLoader(trainset, batch_size=batch_size, shuffle=True)
        testloader = torch.utils.data.DataLoader(testset, batch_size=inp_np_test.shape[0], shuffle=True)
        
        # Loss Functions
        loss_fn = torch.nn.MSELoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=100, verbose=True)
        
        tplot = []
        lossplottrain = []
        lossplottest = []
        
        for t in range(epochs):
            # Before the backward pass, use the optimizer object to zero all of the
            # gradients for the variables it will update (which are the learnable
            # weights of the model). This is because by default, gradients are
            # accumulated in buffers( i.e, not overwritten) whenever .backward()
            # is called. Checkout docs of torch.autograd.backward for more details.
            batch_losses = []
            for inputs, labels in trainloader:
                inputs, labels =  inputs.to(device), labels.to(device)
                optimizer.zero_grad()
                # forward + backward + optimize
                outputs = model(inputs)
                loss = loss_fn(outputs, labels)
                loss2 = self.weighted_loss(outputs,labels)
                loss.backward()
                optimizer.step()
                batch_losses.append(loss.item())
            loss = np.mean(batch_losses)
            
            # evaluate accuracy on test set
            batch_test_losses = []
            model.eval()
            for inputs_test, labels_test in testloader:
                inputs_test, labels_test =  inputs_test.to(device), labels_test.to(device)
                outputs_test = model(inputs_test)
                test_loss = loss_fn(outputs_test, labels_test)
                test_loss2 = self.weighted_loss(outputs_test, labels_test)
                batch_test_losses.append(test_loss.item())
            test_loss = np.mean(batch_test_losses)
            if t % 100 == 99:
                print(t, "loss on train" , loss)
                print(t, "loss on test" , test_loss)
                tplot.append(t)
                lossplottrain.append(loss)
                lossplottest.append(test_loss)
            # if(loss < 0.0003):
            #     scheduler.step(test_loss)
            model.train()

        return tplot,lossplottrain,lossplottest,model,scaling

    def train(self):
        self.t_plot, self.lossplottrain, self.lossplottest,self.model,self.scaling = self.Net1(self.inp_np_train_a, self.inp_np_test_a, self.out_np_train, self.out_np_test, 
                                                              self.Epochs, 1e-4, 1e-7, 5000)
    def plot_loss(self):
        plt.plot(self.t_plot, self.lossplottrain, label="train")
        plt.plot(self.t_plot, self.lossplottest, label="test")
        plt.legend()
        plt.show()

    def save(self):
        m = torch.jit.script(self.model)
        np.save("ITHACAoutput/NN/minAnglesInp_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_inp"].min_[:,None])
        np.save("ITHACAoutput/NN/scaleAnglesInp_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_inp"].scale_[:,None])
        np.save("ITHACAoutput/NN/minOut_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_out"].min_[:,None])        
        np.save("ITHACAoutput/NN/scaleOut_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_out"].scale_[:,None])
        m.save("ITHACAoutput/NN/Net_"+str(self.Nu) + "_" +str(self.Nnut)+".pt")
        np.save("ITHACAoutput/NN/trainLoss_"+str(self.Nu) + "_" +str(self.Nnut)+ ".npy",self.lossplottrain)
        np.save("ITHACAoutput/NN/testLoss_"+str(self.Nu) + "_" +str(self.Nnut)+ ".npy",self.lossplottest)

# sed_variable("NmodesUproj", "system/ITHACAdict", 10)

# Netok = Net(10,10,500)
# Netok.read()
# Netok.train()
# Netok.plot_loss()
# Netok.save()
