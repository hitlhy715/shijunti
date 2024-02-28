# -*- coding = utf-8 -*-
# 创建时间：2023/5/5
# 最后更新时间：2023/5/5
# 作者:lhy

import os
import sys
import numpy as np
import math
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from fe_predphi_protein import fe_predphi_protein
from config import *
from pconfig import *
import os
        
class obtain_feature:
    def __init__(self, type, decode, interact, src_phage, src_bac, fea_phage, fea_bac, mix, minp, maxp, minmax, norm, all_standard):
        self.type = type  # faa/fasta 代表了 碱基对还是protein
        self.decode = decode
        self.interact = interact
        self.src_phage = src_phage
        self.src_bac = src_bac
        self.fea_phage = fea_phage
        self.fea_bac = fea_bac
        self.mix = mix
        self.minp = minp
        self.maxp = maxp
        self.minmax = minmax
        self.norm = norm
        self.all_standard = all_standard
        if self.decode == 'predphi':
            self.decoder = fe_predphi_protein
        
    def getlists(self, path):
        return os.listdir(path)

    def deal_filename(self, path):
        for i in self.getlists(path):
            new_name = i.split('_')[0] + '.faa'
            os.rename(os.path.join(path, i), os.path.join(path, new_name))
        
        
    # 小数的string类型转int类型           
    def dstring2int(self, j):
        if len(j.split('e')) == 2:
            # print(j)
            j2 = j.split('e')[0]
            p = int(j.split('e')[1])
            inter = j2.split('\n')[0].split('.')[0]
            dot = j2.split('\n')[0].split('.')[1]
            nu = int(inter) + math.pow(0.1, len(dot)) * int(dot)
            nu = nu * math.pow(0.1, abs(p))
            # print(nu)
        else:
            inter = j.split('\n')[0].split('.')[0]
            dot = j.split('\n')[0].split('.')[1]
            nu = int(inter) + math.pow(0.1, len(dot)) * int(dot)
        return nu
    
    def getname(self,name):
        if len(name.split('.')) == 2:
            return name.split('.')[0]
        if len(name.split('.')) == 3:
            return name.split('.')[0] + '.'+ name.split('.')[1]
        if len(name.split('.')) == 4:
            return name.split('.')[0] + '.'+ name.split('.')[1] + '.'+ name.split('.')[2] 
        if len(name.split('.')) == 5:
            return name.split('.')[0] + '.'+ name.split('.')[1] + '.'+ name.split('.')[2] + '.'+ name.split('.')[3] 
        
    # 细菌宿主 氨基酸序列 特征提取
    def create_bac_features(self):
        os.makedirs(self.fea_bac, exist_ok=True)
        for f in self.getlists(self.src_bac):
            fpath = os.path.join(self.src_bac, f)
            respath = os.path.join(self.fea_bac, (self.getname(f) + '.csv'))
            self.decoder(fpath,respath)
            
    # 噬菌体 氨基酸序列 特征提取    
    def create_phage_features(self):
        os.makedirs(self.fea_phage, exist_ok=True)
        for f in self.getlists(self.src_phage):
            fpath = os.path.join(self.src_phage, f)
            respath = os.path.join(self.fea_phage, (self.getname(f) + '.csv'))
            self.decoder(fpath,respath)
    
    # 将细菌宿主与噬菌体特征混合 生成m*n个文件
    def create_mix_features(self):
        os.makedirs(self.mix, exist_ok=True)
        for i in self.getlists(self.fea_bac):
            for j in self.getlists(self.fea_phage):
                resp = self.getname(i) + '#' + self.getname(j) + '.csv'
                bacp = open(os.path.join(self.fea_bac, i),'r')
                phagep = open(os.path.join(self.fea_phage, j), 'r')
                
                bacs = bacp.readlines()
                phages = phagep.readlines()
                
                bacp.close()
                phagep.close()
                
                bacs.extend(phages)
                with open(os.path.join(self.mix, resp), 'w') as m:
                    m.writelines(bacs)
               
    # 统计最大值与最小值     
    def generate_minmax(self):
        os.makedirs(self.minmax, exist_ok=True)
        if self.decode == 'predphi':
            maxs = np.zeros((predphi_AAC_nums+predphi_che_nums+predphi_mw_nums)*predphi_type_nums*2)
            mins = np.zeros((predphi_AAC_nums+predphi_che_nums+predphi_mw_nums)*predphi_type_nums*2) + 999999999
        for i in self.getlists(self.mix):
            with open(os.path.join(self.mix, i), 'r') as mf:
                ds = mf.readlines()
            for jidx in range(len(ds)):
                j = ds[jidx].split('\n')[0]
                nu = self.dstring2int(j)
            
                if nu > maxs[jidx]:
                    maxs[jidx] = nu
                if nu < mins[jidx]:
                    mins[jidx] = nu
        with open(self.minp, 'w') as minf:
            for ms in mins:
                minf.write(str(ms) + '\n')
        with open(self.maxp, 'w') as maxf:
            for ms in maxs:
                maxf.write(str(ms) + '\n')
    
    # 最大最小值归一
    def maxminNormalize(self, data, max, min):
        rdata = np.zeros(len(data))
        for idx in range(len(data)):
            if float(max[idx]) - float(min[idx]) == 0:
                rdata[idx] = 0
            elif float(max[idx]) - float(data[idx]) <= 1e-15:
                rdata[idx] = 1
            elif float(data[idx]) - float(min[idx]) <= 1e-15:
                rdata[idx] = 0
            else:
                assert float(max[idx]) >= float(min[idx]), max[idx] + ';' + min[idx]
                assert float(data[idx]) >= float(min[idx]), data[idx] + ';' + min[idx]
                rdata[idx] = ((float(data[idx]) - float(min[idx])) / (float(max[idx]) - float(min[idx])))
                assert rdata[idx] >=0 and rdata[idx] <= 1, max[idx] + ';' + min[idx] + ';' + data[idx]
        return rdata
    
    # 标准化混合后的所有特征
    def norm_feature(self):
        os.makedirs(self.norm, exist_ok=True)
        for i in self.getlists(self.mix):
            with open(os.path.join(self.mix, i), 'r') as od:
                ds = od.readlines()
            with open(self.maxp, 'r') as maxd:
                maxs = maxd.readlines()
            with open(self.minp, 'r') as mind:
                mins = mind.readlines()
                
            ds   = [_.split('\n')[0] for _ in ds]
            maxs = [_.split('\n')[0] for _ in maxs]
            mins = [_.split('\n')[0] for _ in mins]
            # minmax归一化方法
            norms = self.maxminNormalize(ds, maxs, mins)
            
            with open(os.path.join(self.norm, i), 'w') as rd:
                for i in norms:
                    rd.write(str(i) + '\n')
    
    # 将格式调整为prePHI源代码要求的输入格式，n行，每行324+1维    
    def integrate(self):
        all_data = []
        for i in self.getlists(self.norm):
            with open(os.path.join(self.norm, i), 'r') as od:
                ds = od.readlines()
                ds = [_.split('\n')[0] for _ in ds]
                ds = '\t'.join(ds)
                all_data.append([i ,ds])
        labels = {}
        with open(self.interact, 'r') as hd:
            labelss = hd.readlines()
            labelss = [_.split('\n')[0] for _ in labelss]
            labelss = [_.split('\t') for _ in labelss]
            new_labelss  = []
            for l in labelss:
                new_l0 = l[0].split(' ')
                new_l1 = l[1].split(' ')
                new_labelss.append(['_'.join(new_l0), '_'.join(new_l1),l[2]])
                
            labelss = [[_[1] + '#' + _[0], _[2]] for _ in new_labelss]
        for line in labelss:
            labels[line[0]] = line[1]
        print(labels)
            
        with open(self.all_standard, 'w') as pd:
            for line in all_data:
                name = self.getname(line[0])
                print(name,line[0])
                pd.write(line[1] + '\t' + labels[name] + '\n')
            
if __name__ == '__main__':
    obtainer = obtain_feature(type='faa', 
                              decode='predphi',
                              interact=allinteract_file1,
                              src_phage=phagepro_faa,
                              src_bac=bacpro_faa,
                              fea_phage=fed_phage_predphi_protein1,
                              fea_bac=fed_bac_predphi_protein1,
                              mix = fed_mix_predphi_protein1,
                              minp=fed_min_predphi_protein1,
                              maxp=fed_max_predphi_protein1,
                              minmax=fed_minmax_predphi_protein1,
                              norm=fed_norm_predphi_protein1,
                              all_standard=fed_all_predphi_protein1)
    obtainer.create_bac_features()
    obtainer.create_phage_features()
    obtainer.create_mix_features()
    obtainer.generate_minmax()
    obtainer.norm_feature()
    obtainer.integrate()
    
    # obtainer = obtain_feature(type='faa', 
    #                           decode='predphi',
    #                           interact=allinteract_file2,
    #                           src_phage=phagepro_faa2,
    #                           src_bac=bacpro_faa2,
    #                           fea_phage=fed_phage_predphi_protein2,
    #                           fea_bac=fed_bac_predphi_protein2,
    #                           mix = fed_mix_predphi_protein2,
    #                           minp=fed_min_predphi_protein2,
    #                           maxp=fed_max_predphi_protein2,
    #                           minmax=fed_minmax_predphi_protein2,
    #                           norm=fed_norm_predphi_protein2,
    #                           all_standard=fed_all_predphi_protein2
    #                           )
    # obtainer.create_bac_features()
    # obtainer.create_phage_features()
    # obtainer.create_mix_features()
    # obtainer.generate_minmax()
    # obtainer.norm_feature()
    # obtainer.integrate()

                

    
        
