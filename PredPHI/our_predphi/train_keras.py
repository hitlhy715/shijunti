# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
from sklearn.model_selection import KFold
from sklearn.metrics import f1_score, accuracy_score, recall_score, precision_score, confusion_matrix, roc_auc_score, matthews_corrcoef,roc_curve
import keras.applications
import matplotlib
matplotlib.use("Agg")
from keras.optimizers import Adam
from sklearn.preprocessing import LabelBinarizer
from keras.callbacks import Callback
import tensorflow as tf
from keras.callbacks import EarlyStopping
from keras_model import PredPHI
from keras.models import load_model

import warnings
warnings.filterwarnings("ignore")

class Metrics(Callback):
    def __init__(self, validation_data):
        super(Callback, self).__init__()
        self.validation_data = validation_data
    
    def on_train_begin(self, logs={}):
        self.val_f1s = []
        self.val_recalls = []
        self.val_precisions = []
 
    def on_epoch_end(self, epoch, logs={}):
        val_predict = (np.asarray(self.model.predict(self.validation_data[0]))).round()
        val_targ = self.validation_data[1]
        _val_f1 = f1_score(val_targ, val_predict,average='micro')
        _val_recall = recall_score(val_targ, val_predict,average=None)
        _val_precision = precision_score(val_targ, val_predict,average=None)
        self.val_f1s.append(_val_f1)
        self.val_recalls.append(_val_recall)
        self.val_precisions.append(_val_precision)
        return _val_f1
    
class keras_predphi:
    def __init__(self, datapath, epochs, init_lr, batch_size, feature_dim, k_splits, n_class, model_path, loss, score_thres):
        self.datapath = datapath
        self.epochs = epochs
        self.init_lr = init_lr
        self.batch_size = batch_size
        self.feature_dim = feature_dim
        self.k_splits = k_splits
        self.n_class = n_class
        self.model_path = model_path
        self.loss = loss
        self.score_thres = score_thres
        
    def kfolder(self):
        data=pd.read_csv(self.datapath ,header=None,sep='\t').values
        print(data.shape)
        ##kfold result
        kf = KFold(self.k_splits, random_state=1, shuffle=True)
        result_predphi=[]
        predphi_pred=[]
        test_y_all=[]
        for train_index, test_index in kf.split(data): 
            training=data[train_index,:]
            test=data[test_index,:]   
            feature_training,label_training=self.obtainfeature(training)
            feature_test,label_test=self.obtainfeature(test)
            feature_training2=np.array(feature_training).transpose(0,2,3,1)
            print(feature_training2.shape)
            
            label_training2 = keras.utils.to_categorical(label_training)
            feature_test2=np.array(feature_test).transpose(0,2,3,1)
            test_y_all=test_y_all+label_test.tolist()
            model = PredPHI().build(width=self.feature_dim[0], height=self.feature_dim[1], depth=self.feature_dim[2], classes=self.n_class)
            # for i in range(len(model.layers)):
            #     print(i, model.layers[i].input_shape, model.layers[i].output_shape)
                
            opt = Adam(lr=self.init_lr, decay=self.init_lr/self.epochs)
            model.compile(loss=self.loss, optimizer=opt, metrics=['acc'])
            f1=Metrics((feature_training2, label_training2))
            model.fit(feature_training2, label_training2, batch_size=self.batch_size,epochs=self.epochs,validation_data=(feature_training2, label_training2), verbose=1,
                        callbacks=[EarlyStopping(monitor='val_loss', patience=20, verbose=2, mode='auto'),f1])
            # res = model.predict(feature_test2)
            # print(res)
            # print(res.shape)
            # print(model.predict)
            result_predphi.append(self.scores(label_test,model.predict(feature_test2)[:,1]))
            predphi_pred=predphi_pred+model.predict(feature_test2)[:,1].tolist()  
        np.savetxt('result/pred-10fold.csv',np.array([test_y_all,predphi_pred]).T)
        result=np.mean(result_predphi,axis=0)
        with open('result/result_score-10fold.csv','w') as fout:
            fout.write('sen,spe,pre,f1,mcc,acc,auc,tn,fp,fn,tp\n')
            for jj in result:
                fout.write(str(jj)+',')
            fout.write('\n')

    def construct_model(self):
        ##construct model  
        data=pd.read_csv(self.datapath,sep='\t',header=None).values
        feature,label = self.obtainfeature(data)
        print(np.array(feature).shape)
        print(np.array(label).shape)
        feature2=np.array(feature).transpose(0,2,3,1)
        label2 = keras.utils.to_categorical(label)
        model = PredPHI().build(width=self.feature_dim[0], height=self.feature_dim[1], depth=self.feature_dim[2], classes=self.n_class)
        opt = Adam(lr=self.init_lr, decay=self.init_lr/self.epochs)
        model.compile(loss=self.loss, optimizer=opt, metrics=['acc'])
        f1=Metrics((feature2, label2))
        model.fit(feature2, label2, batch_size=self.batch_size,epochs=self.epochs,validation_data=(feature2, label2), verbose=1,
                    callbacks=[EarlyStopping(monitor='val_loss', patience=20, verbose=2, mode='auto'),f1])
        model.save(self.model_path)
        
    def obtainfeature(self, data):
        feature=data[:,:-1]
        label=data[:,-1]
        feature_phage=feature[:,:int(feature.shape[1]/2)]
        feature_host=feature[:,int(feature.shape[1]/2):]
        
        ##obtain phage features
        # 拼接顺序： AAC + CHONS + weight
        feature_phage_CHONS=[]
        feature_phage_weight=[]
        feature_phage_AAC=[]
        for i in feature_phage:
            aa=i[:21]  #obtain AAC
            cc=i[21:26]  #obtain CHONS
            ww=[i[26]]    #obtain weight
           
            for ii in range(27,len(i),27):   ##obtain mean,max,min...
                aa=np.concatenate((aa,i[ii:ii+21]))
                cc=np.concatenate((cc,i[ii+21:ii+26]))
                ww.append(i[ii+26])
                
            feature_phage_CHONS.append(cc.tolist())
            feature_phage_weight.append(ww)
            feature_phage_AAC.append(aa.tolist())
        # print(np.array(feature_phage_CHONS).shape)
        # print(np.array(feature_phage_weight).shape)
        # print(np.array(feature_phage_AAC).shape)
        
        ##obtain host features
        feature_host_CHONS=[]
        feature_host_weight=[]
        feature_host_AAC=[]
        for i in feature_host:
            aa=i[:21]  #obtain AAC
            cc=i[21:26]  #obtain CHONS
            ww=[i[26]]    #obtain weight
            for ii in range(27,len(i),27):
                aa=np.concatenate((aa,i[ii:ii+21]))
                cc=np.concatenate((cc,i[ii+21:ii+26]))
                ww.append(i[ii+26])
            feature_host_CHONS.append(cc.tolist())
            feature_host_weight.append(ww)
            feature_host_AAC.append(aa.tolist())
            
        ##combine phage and host features
        feature_CHONS=np.concatenate((feature_phage_CHONS,feature_host_CHONS),axis=1)
        feature_AAC=np.concatenate((feature_phage_AAC,feature_host_AAC),axis=1)  
        feature_weight = np.concatenate((feature_phage_weight, feature_host_weight), axis=1)
        # print(feature_CHONS.shape)
        # print(feature_AAC.shape)
        # print(feature_weight.shape)
        
        feature_CHONS_AAC=np.concatenate((feature_CHONS,feature_AAC, feature_weight),axis=1)
        feature_CHONS_AAC_new=[]
        for i in range(len(feature_CHONS_AAC)):
            feature_CHONS_AAC_new.append([feature_CHONS_AAC[i,:int(feature_CHONS_AAC.shape[1]/2)].reshape((6,-1)),
                                            feature_CHONS_AAC[i,int(feature_CHONS_AAC.shape[1]/2):].reshape((6,-1))])
        return feature_CHONS_AAC_new,label
    
    def scores(self, y_test, y_pred):           
        y_predlabel = [(0. if item < self.score_thres else 1.) for item in y_pred]
        tn, fp, fn, tp = confusion_matrix(y_test, y_predlabel).flatten()
        SPE = tn*1./(tn+fp)
        MCC = matthews_corrcoef(y_test, y_predlabel)
        fpr,tpr,threshold = roc_curve(y_test, y_predlabel)
        sen, spe, pre, f1, mcc, acc, auc, tn, fp, fn, tp = np.array([recall_score(y_test, y_predlabel), SPE, precision_score(y_test, y_predlabel,average='macro'), 
                                                                    f1_score(y_test, y_predlabel), MCC, accuracy_score(y_test, y_predlabel), 
                                                                    roc_auc_score(y_test, y_pred), tn, fp, fn, tp])
        return sen, spe, pre, f1, mcc, acc,auc,tn,fp,fn,tp  

    def list_of_groups(init_list, childern_list_len):
        list_of_group = zip(*(iter(init_list),) *childern_list_len)
        end_list = [list(i) for i in list_of_group]
        count = len(init_list) % childern_list_len
        end_list.append(init_list[-count:]) if count !=0 else end_list
        return end_list
    
    def test(self, test_path):
        test_kmeans=pd.read_csv(test_path,header=None,sep='\t').values
        print(test_kmeans.shape)
        feature_kmeans_CHONS_AAC_new, label_new_kmeans =self.obtainfeature(test_kmeans)
        print(os.listdir('.'))
        model = load_model(self.model_path) 
        test_X=np.array(feature_kmeans_CHONS_AAC_new).transpose(0,2,3,1)
        print(self.scores(label_new_kmeans, model.predict(test_X)[:,1]))
