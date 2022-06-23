import torch
from torch import nn
from torch.autograd import Variable
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn import preprocessing
from sklearn.metrics import r2_score
from torch.utils.data import Dataset


class NetVelocityNut(nn.Module):
    def __init__(self, in_features, nb_classes,  
        hidden_sizes, act):
        super(NetVelocityNut, self).__init__()
        self.actName = act
        self.act = act()
        self.layers = []
        self.hidden_sizes = hidden_sizes
        self.in_features = in_features
        self.nb_classes = nb_classes
        self.layers.append(nn.Linear(in_features, hidden_sizes[0]))
        self.layers.append(act())
        for i, k in enumerate(hidden_sizes[:-1]):
             self.layers.append(nn.Linear(hidden_sizes[i], hidden_sizes[i+1]))
             self.layers.append(act())
        self.layers.append(nn.Linear(hidden_sizes[-1], nb_classes))
        self.net = nn.Sequential(*self.layers)

    def forward(self, x):
        x = self.net.forward(x)
        return x

    def hash(self):
        name = str(self.in_features)
        for i in self.hidden_sizes:
            name+= "_"+str(i)
        name+= "_"+str(self.nb_classes)
        name+= ("_"+str(self.act))[:-2]
        return name


class NutDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.from_numpy(X).type(torch.float32)
        self.y = torch.from_numpy(y).type(torch.float32)
 
    def __len__(self):
        return len(self.X)
 
    def __getitem__(self, idx):
        return [self.X[idx], self.y[idx]]

#a = NetVelocityNut(5,5,[10,10,10],act=nn.ReLU)
#print(a)
