import torch.nn as nn

# keras中维度为 batch，height，width，channels
# pytorch中维度为 batch channels height width
# 原论文的特征维度为（BS，6，27，2）
# pytorch中输入维度应该为 (BS, 2, 6, 27)

def get_activation(name="silu", inplace=False):
    if name == "silu":
        module = nn.SiLU(inplace=inplace)
    elif name == "relu":
        module = nn.ReLU(inplace=inplace)
    elif name == "lrelu":
        module = nn.LeakyReLU(0.1, inplace=inplace)
    elif name == "softmax":
        module = nn.Softmax(dim=1)
    else:
        raise AttributeError("Unsupported act type: {}".format(name))
    return module

class PredPHI(nn.Module):
    def __init__(self, shape, act, fact, classes, dropout):
        super().__init__()
        
        self.activation = get_activation(act)
        self.factivation = get_activation(fact)
        
        self.conv1 = nn.Sequential(
            nn.Conv2d(in_channels=shape[0], out_channels=32, kernel_size=3, padding=1, stride=1),
            self.activation,
            nn.BatchNorm2d(num_features=32),
        )
        
        self.pool1 = nn.Sequential(
            nn.MaxPool2d(kernel_size=3),
            nn.Dropout(dropout)
        )
        
        self.FC1 = nn.Sequential(
            nn.Flatten(),
            nn.Linear(576,1024),
            self.activation,
            nn.BatchNorm1d(num_features=1024),
            nn.Dropout(dropout)  
        )
        
        self.FC2 = nn.Sequential(
            nn.Linear(1024, classes),
            self.factivation
        )
    
    def forward(self, x):
        # 输入维  torch.Size([3, 2, 6, 26])
        x = self.conv1(x)
        # 输入维  torch.Size([3, 32, 6, 26])
        x = self.pool1(x)
        # 输入维torch.Size([3, 32, 2, 8])
        x = self.FC1(x)
        # torch.Size([3, 1024])
        x = self.FC2(x)
        # 输出结果 torch.Size([3, 2])
        return x