import torch
import psd
import argparse
import network_cuda
import matplotlib.pyplot as plt
import numpy as np
from inference import inference_adaptive, prep, to_numpy

'''

Routine that generates a timeseries and then computes the PSD between
the generated timeseries and the test data. Multiple parameters including
the smoothing hyperparameter sigma can be changed.

'''
parser = argparse.ArgumentParser(description='routine parameters')
parser.add_argument('--model_path', type=str, default="",
                    help='model path')
parser.add_argument('--nhid', type=int, default=16,
                    help='hidden size')
parser.add_argument('--lorenz_type', type=str, default="lorenz63",
                    help='lorenz type')
parser.add_argument('--bootstrapT', type=int, default=1,
		    		help='how long of a sequence should the model get')
parser.add_argument('--T', type=int, default=1000,
		    		help='prediction length')
parser.add_argument('--sigma', type=int, default=20)
parser.add_argument('--generate', type=bool, default=True, help='generate the data for psd or not?')

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
args = parser.parse_args()
model_path = args.model_path
ninp = 3 if args.lorenz_type == "lorenz63" else 20
nhid = args.nhid
nout = ninp
lorenz_type = args.lorenz_type
bootstrapT = args.bootstrapT
generate = args.generate
psd.SMOOTHING_SIGMA = args.sigma

test_data = np.load('data/{}_test.npy'.format(lorenz_type))

if generate:
	# load model from file and into the device
	model = network_cuda.LEM(ninp, nhid, nout).to(device)
	if model_path != "":
		model.load_state_dict(torch.load(model_path))

	T_pred = args.T

	# generate the timeseries
	x0 = np.random.normal(0, 1, (1, ninp))
	y = to_numpy(inference_adaptive(model, prep(x0).to(device), T_pred + 2 * bootstrapT, bootstrapT))[2*bootstrapT:].reshape(-1, ninp)

	# Plot for reference
	fig = plt.figure(figsize=(6, 6))
	ax = fig.add_subplot(111, projection='3d')
	ax.plot(y[:,0], y[:,1], y[:,2], label='prediction', color='r', zorder=0)
	ax.legend()
	ax.set_xlabel(r'$x$', fontsize=16)
	ax.set_ylabel(r'$y$', fontsize=16)
	ax.set_zlabel(r'$z$', fontsize=16)
	plt.show()

	#save y
	np.save('data/{}_gen.npy'.format(lorenz_type), y)

# Calculate power-spectrum-distance
print(psd.power_spectrum_error(test_data.reshape(-1, 1, ninp), y.reshape(-1, 1, ninp)))