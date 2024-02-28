from keras.models import Sequential
from keras.layers.normalization.batch_normalization_v1 import BatchNormalization
from keras.layers.convolutional import Conv2D
from keras.layers.convolutional import MaxPooling2D
from keras.layers.core import Activation
from keras.layers.core import Flatten
from keras.layers.core import Dropout
from keras.layers.core import Dense
from keras import backend as K

class PredPHI:
    # width:6   height:26  depth:2
    def build(self,width, height, depth, classes, finalAct="softmax"):
        model = Sequential()
        inputShape = (width, height, depth)
        chanDim = -1
        
        if(K.image_data_format == "channels_first"):
            inputShape = (depth, height, width)  
            chanDim = 1
        
        print(inputShape)
        model.add(Conv2D(32, (3,3), padding="same", input_shape=inputShape))
        model.add(Activation("relu"))
        model.add(BatchNormalization(axis=chanDim))
        model.add(MaxPooling2D(pool_size=(3,3)))
        model.add(Dropout(0.5))
              
        model.add(Flatten())
        model.add(Dense(1024))
        model.add(Activation("relu"))
        model.add(BatchNormalization())
        model.add(Dropout(0.5))
        
        model.add(Dense(classes))
        model.add(Activation(finalAct))
        
        return model