from pathlib import Path
import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset
from torch import Tensor
import os
from matplotlib import pyplot as plt

def save_model_and_params(model, valid_losses, args, timetaken, final_test_loss):
	path = Path('models')
	path.mkdir(parents=True, exist_ok=True)
	i = 0
	while True:
		if not (path / str(i)).exists():
			break
		i += 1
	path = path / str(i)
	path.mkdir(parents=True, exist_ok=True)
	torch.save(model.state_dict(), path / 'model.pt')
	np.savetxt(path / 'valid_losses.txt', np.array(valid_losses))
	with open(path / 'args.txt', 'w') as f:
		f.write(str(args))

	with open(path / 'timetaken.txt', 'w') as f:
		f.write(str(timetaken))
	
	with open(path / 'final_test_loss.txt', 'w') as f:
		f.write(str(final_test_loss))

def slice_data(data, seq_len, stride=1):
    stack_x = []
    stack_y = []
    for i in range(0, len(data) - seq_len, stride):
        stack_x.append(data[i:i+seq_len, :])
        stack_y.append(data[i+1:i+seq_len+1, :])
    return np.stack(stack_x), np.stack(stack_y)

def create_dataloader(data, seq_len, batch_size, stride=1):
    x, y = slice_data(data, seq_len, stride)
    x = Tensor(x).float()
    y = Tensor(y).float()
    return DataLoader(TensorDataset(x, y), batch_size=batch_size, shuffle=True)

def parse_args(path):
	args_path = os.path.join(path, 'args.txt')
	args_dict = {}
	with open(args_path, 'r') as f:
		args = f.read()
		args = args.replace('Namespace(', '')[:-1].split(', ')
		for arg in args:
			key, val = arg.split('=', 1)
			args_dict[key] = val
	return args_dict

def load_losses_and_args(path):
    losses = np.loadtxt(path + 'valid_losses.txt')
    return parse_args(path), losses

# plot losses on left axis and lr on right axis

def plot_training_routine(paths, title):
    fig, ax = plt.subplots()

    val_losses = []
    lens = []
    prev_args = None
    for i, path in enumerate(paths):
        args, losses = load_losses_and_args(path)
        prev = len(val_losses)
        ax.plot(np.arange(prev, prev + len(losses)), losses, color='g', label='Validation loss' if i == 0 else None)
        lens.append(prev)
        val_losses += list(losses)
        i = 0
        for arg in args:
            if ((prev_args is None or (arg in args and arg in prev_args and args[arg] != prev_args[arg]))
                and arg in ['batch', 'lr', 'nhid', 'seq_len', 'stride']):
                plt.text(prev+1, np.max(val_losses)*np.exp(-0.1*i), '{}={}'.format(arg, args[arg]), fontsize=10)
                i += 1
        prev_args = args

    max = 1.1 * np.max(val_losses)
    min = 0.9 * np.min(0.05)
    plt.ylim(min, max)
    plt.hlines(0.05, 0, len(val_losses), color='black', alpha=0.8, label='Optimal loss')
    for i, prev in enumerate(lens):
        plt.vlines(prev, min, max, color='red', alpha=0.3, label='Training restart' if i == 0 else None)
    ax.set_yscale('log')
    ax.set_xlabel('Epoch', fontsize=14)
    ax.set_ylabel('Validation loss', fontsize=14)
    # legend outside of plot
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=14)
    plt.title(title, fontsize=14)
    plt.show()