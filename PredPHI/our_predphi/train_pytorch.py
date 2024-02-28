import sys
import os
import numpy as np
import torch as t
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import torch.nn as nn

from pytorch_model import PredPHI
from pconfig import *


data = t.ones(3, 2, 6, 27)

model = PredPHI(input_shape, act, fact, classes, dropout)
res = model(data)

print(res)