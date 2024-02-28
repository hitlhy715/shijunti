import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from config import *
from pconfig import *

from train_keras import keras_predphi

model = keras_predphi(
    datapath = fed_all_predphi_protein2,
    epochs=150,
    init_lr=1e-4,
    batch_size=32,
    feature_dim = (6,27,2) ,
    k_splits=12,
    n_class=2,
    model_path=predphi_model_out,
    loss="binary_crossentropy",
    score_thres=0.5
)

# model.kfolder()
# model.construct_model()
model.test(fed_all_predphi_protein2)