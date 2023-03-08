import torch
import psd
import argparse
import network_cuda
import matplotlib.pyplot as plt
import numpy as np
from inference import inference, prep, to_numpy

'''

Routine that generates a timeseries and then computes the PSD between
the generated timeseries and the test data. Multiple parameters including
the smoothing hyperparameter sigma can be changed.

'''
parser = argparse.ArgumentParser(description='routine parameters')
parser.add_argument('--model_path', type=str, default="",
                    help='model path')
parser.add_argument('--nhid', type=int, default=15,
		    		help='input size as well as output size')
parser.add_argument('--lorenz_type', type=str, default="lorenz63",
                    help='lorenz type')
parser.add_argument('--bootstrapT', type=int, default=1,
		    		help='how long of a sequence should the model get')
parser.add_argument('--T', type=int, default=1000,
		    		help='prediction length')
parser.add_argument('--title', type=str, default="", help='title for plot')

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
args = parser.parse_args()
model_path = args.model_path
ninp = 3 if args.lorenz_type == "lorenz63" else 20
nhid = args.nhid
nout = ninp
lorenz_type = args.lorenz_type
bootstrapT = args.bootstrapT
title = args.title

# load model from file and into the device
model = network_cuda.LEM(ninp, nhid, nout).to(device)
if model_path != "":
	model.load_state_dict(torch.load(model_path))

test_data = np.load('data/{}_test.npy'.format(lorenz_type))

T_start = 1000
T_pred = args.T

# generate the timeseries
gt = test_data[T_start+bootstrapT:T_start+T_pred+bootstrapT,:]
x0 = test_data[T_start:T_start+bootstrapT]
y = to_numpy(inference(model, prep(x0).to(device), T_pred)).reshape(-1, nout)

# Plot for reference
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(gt[0,0],
	   		gt[0,1],
			gt[0,2],
			color='g',
			zorder=1,
			label='Prediction start')
ax.plot(gt[:, 0],
		gt[:, 1],
		gt[:, 2],
		label='test data',
		color='b', zorder=0)
ax.plot(y[:,0], y[:,1], y[:,2], label='prediction', color='r', zorder=0)
ax.legend()
ax.set_xlabel(r'$x$', fontsize=16)
ax.set_ylabel(r'$y$', fontsize=16)
ax.set_zlabel(r'$z$', fontsize=16)
plt.title(title, fontsize=16)
plt.show()