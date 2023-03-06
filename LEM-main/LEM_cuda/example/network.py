import math
import torch
from torch import nn
import lem_cuda

class LEMFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, inputs, weights, weights_lin_z, bias,
                bias_lin_z, initial_y_state, initial_z_state, dt):
        all_y_states, all_z_states, all_X, all_X2, \
        all_multi_scales, all_lin_new_z_state = lem_cuda.forward(inputs, weights,
                                                                 weights_lin_z, bias, bias_lin_z,
                                                                 initial_y_state, initial_z_state, dt)
        ctx.save_for_backward(all_X, all_X2, all_multi_scales,
                              all_lin_new_z_state, weights, weights_lin_z, bias,
                              bias_lin_z, initial_y_state, initial_z_state, dt)
        return all_y_states, all_z_states

    @staticmethod
    def backward(ctx, grad_y_states, grad_z_states):
        outputs = lem_cuda.backward(grad_y_states.contiguous(), grad_z_states.contiguous(), *ctx.saved_tensors)
        d_inputs, d_weights, d_weights_lin_z, d_bias, d_bias_lin_z, d_initial_y_state, d_initial_z_state = outputs
        return None, d_weights, d_weights_lin_z, d_bias, d_bias_lin_z, d_initial_y_state, d_initial_z_state, None


class LEMcuda(torch.nn.Module):
    def __init__(self, ninp, nhid, dt):
        super(LEMcuda, self).__init__()
        self.ninp = ninp
        self.nhid = nhid
        self.weights = torch.nn.Parameter(torch.empty(3 * nhid, ninp + nhid))
        self.weights_lin_z = torch.nn.Parameter(torch.empty(nhid, ninp + nhid))
        self.bias = torch.nn.Parameter(torch.empty(3 * nhid))
        self.bias_lin_z = torch.nn.Parameter(torch.empty(nhid))
        self.dt = torch.tensor(dt).reshape(1, 1).cuda()

        self.reset_parameters()

    def reset_parameters(self):
        stdv = 1.0 / math.sqrt(self.nhid)
        for weight in self.parameters():
            weight.data.uniform_(-stdv, +stdv)

    def forward(self, input, states=None):
        if states is None:
            y = input.data.new(input.size(1), self.nhid).zero_()
            z = input.data.new(input.size(1), self.nhid).zero_()
            states = (y,z)
        return LEMFunction.apply(input.contiguous(), self.weights,
                                 self.weights_lin_z, self.bias,
                                 self.bias_lin_z, *states, self.dt)

class LEM(torch.nn.Module):
    def __init__(self, ninp, nhid, nout, dt=1.):
        super(LEM, self).__init__()
        self.ninp = ninp
        self.nhid = nhid
        self.rnn = LEMcuda(ninp,nhid,dt)
        self.classifier = nn.Linear(nhid,nout)

        self.init_weights()

    def init_weights(self):
        for name, param in self.named_parameters():
            if 'classifier' in name and 'weight' in name:
                nn.init.kaiming_normal_(param.data)

    def forward(self, input):
        all_y, all_z = self.rnn(input)
        out = self.classifier(all_y)
        return out