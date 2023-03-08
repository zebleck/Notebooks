import torch

# keeps the shape of the input to the model constant
def inference(model, x0, T):
	model.eval()
	x = x0
	xs = []
	for i in range(T):
		x_prev = x
		x = model(x)
		xs.append(x[-1, :, :])
		x[:-1, :, :] = x_prev[1:, :, :]

		if i % 200 == 0:
			print('Progress: {}/{}'.format(i, T))
	return torch.stack(xs)

# Extends the input length until it reaches seq_len even if the length of the input is lower than that.
# Leads to more stable results since usually a model is trained with a fixed input length (usually >1).
def inference_adaptive(model, x0, T, seq_len):
	assert seq_len > x0.shape[0]
	model.eval()
	x = x0
	xs = []
	for i in range(T):
		x_prev = x
		x = model(x)
		xs.append(x[-1, :, :])
		x[:-1, :, :] = x_prev[1:, :, :]
		# adapt x
		if x.shape[0] < seq_len:
			x = torch.cat((x_prev, x), dim=0)
		if x.shape[0] > seq_len:
			x = x[-seq_len:, :, :]

		if i % 200 == 0:
			print('Progress: {}/{}'.format(i, T))
	return torch.stack(xs)


def prep(x):
	# if dim = 1, reshape to (1, -1, 1) tensor
	if len(x.shape) == 1:
		x = x.reshape(1, -1, 1)
	elif len(x.shape) == 2:
		x = x.reshape(1, -1, x.shape[1])

	x = torch.tensor(x).permute(1, 0, 2).float()
	return x

def to_numpy(x):
	return x.detach().cpu().numpy()