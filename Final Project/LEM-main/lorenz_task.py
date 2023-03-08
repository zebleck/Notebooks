# %load lorenz_task.py
from torch import nn, optim
import numpy as np
import torch.nn.utils
import argparse
import network_cuda
from util import save_model_and_params, slice_data, create_dataloader
import time

parser = argparse.ArgumentParser(description='training parameters')

parser.add_argument('--nhid', type=int, default=16,
                    help='hidden size')
parser.add_argument('--epochs', type=int, default=400,
                    help='max epochs')
parser.add_argument('--device', type=str, default=torch.device('cuda' if torch.cuda.is_available() else 'cpu'),
                    help='computing device')
parser.add_argument('--batch', type=int, default=32,
                    help='batch size')
parser.add_argument('--lr', type=float, default=0.00904,
                    help='learning rate')
parser.add_argument('--seed', type=int, default=1234,
                    help='random seed')
parser.add_argument('--seq_len', type=int, default=128,
                    help='length of the sequence')
parser.add_argument('--stride', type=int, default=5,
                    help='stride of the sequence')
parser.add_argument('--model_path', type=str, default="",
                    help='model path')
parser.add_argument('--lorenz_type', type=str, default="lorenz63",
                    help='lorenz type')

args = parser.parse_args()
print(args)

nhid=args.nhid
epochs=args.epochs
device=args.device
batch=args.batch
lr=args.lr
seed=args.seed
seq_len=args.seq_len
stride=args.stride
model_path=args.model_path
ninp=3 if args.lorenz_type == "lorenz63" else 20
lorenz_type=args.lorenz_type

train_val_ratio = 0.8
train_data = np.load('data/{}_on0.05_train.npy'.format(lorenz_type))
print(train_data.shape)
val_data = train_data[int(train_val_ratio * len(train_data)):]
train_data = train_data[:int(train_val_ratio * len(train_data))]
test_data = np.load('data/{}_test.npy'.format(lorenz_type))

trainloader = create_dataloader(train_data, seq_len, batch, stride=stride)
validloader = create_dataloader(val_data, seq_len, 128, stride=stride)
testloader = create_dataloader(test_data, seq_len, 128, stride=stride)

print('Number of train batches:', len(trainloader))
print('Shape of train batch:', list(next(iter(trainloader))[0].shape))
print('Number of validation batches:', len(validloader))
print('Shape of validation batch:', list(next(iter(validloader))[0].shape))
print('Number of test batches:', len(testloader))
print('Shape of test batch:', list(next(iter(testloader))[0].shape))

nout = ninp

torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
torch.manual_seed(seed)
np.random.seed(seed)


model = network_cuda.LEM(ninp, nhid, nout).to(device)

if model_path != "":
    model.load_state_dict(torch.load(model_path))

objective = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=lr)

start = time.time()

print('Training...')

def eval(dataloader):
    model.eval()
    with torch.no_grad():
        for x,y in dataloader:
            y = y.permute(1, 0, 2)
            x = x.permute(1, 0, 2)
            out = model(x.to(device))
            loss = torch.sqrt(objective(out,y.to(device))).item()
    return loss

valid_losses = []

best_loss = 10000
for epoch in range(epochs):
    model.train()
    for x,y in trainloader:
        y = y.permute(1, 0, 2)
        x = x.permute(1, 0, 2)
        optimizer.zero_grad()
        out = model(x.to(device))
        loss = objective(out, y.to(device))
        loss.backward()
        optimizer.step()
    valid_loss = eval(validloader)
    valid_losses.append(valid_loss)
    test_loss = eval(testloader)
    if (valid_loss < best_loss):
        best_loss = valid_loss
        final_test_loss = test_loss

    log = 'Epoch: {:03d}, Val loss: {:.6f}'
    print(log.format(epoch+1, valid_loss))

print('Final test loss: ' + str(final_test_loss))

print('Time taken: ' + str(time.time() - start) + ' seconds')

save_model_and_params(model, valid_losses, args, str(time.time() - start), final_test_loss)