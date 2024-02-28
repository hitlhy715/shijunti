# -*- coding = utf-8 -*-
# 创建时间：2023/5/5
# 最后更新时间：2023/5/5
# 作者:lhy

import os
import sys
import numpy as np
from collections import Counter
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from utils.fileOperater import readTxtInfo
from utils.infoStatistic import getProteinKinds
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# 输入：
#   filepath: 待提取的氨基酸序列文件路径，内部包含分段的氨基酸序列
#   respath:  目的文件路径，将处理后的特征存储在该文件路径
# 返回：无
# 功能： 实现特征提取，针对氨基酸序列提取，算法为PrePHI论文提供的提取编码方法

def fe_predphi_protein(faapath, respath):
    # 读取氨基酸序列信息，以node进行拆分
    sequences = readTxtInfo(faapath)

    # 统计当前文件中存在的所有氨基酸类型和数量
    kinds = getProteinKinds(sequences)
    print(kinds)

    features = []

    for seq in sequences:
        if seq.endswith('*'):
            seq = seq[:-1]
        # 先提取AAC特征
        AACfeatures = AAC(seq)
        # 提取化学元素丰度
        chemicalFeatures = chemicals(seq)
        # 提取分子量特征
        mwfeature = mw(seq)
        
        # 拼接三种特征，写入文件中
        features.append(AACfeatures+chemicalFeatures+mwfeature)
    
    assert len(features) == len(sequences), str(len(features))+str(len(seq))
    # print(features[0])
    # print(len(features))
    # print(len(features[0]))
    
    # 计算六种特征，拼成27*6=162维
    phifeatures = np.array(features)
    phifeatures = [phifeatures.mean(axis=0).tolist(),phifeatures.max(axis=0).tolist(),phifeatures.min(axis=0).tolist(),
                 phifeatures.std(axis=0).tolist(),phifeatures.var(axis=0).tolist(),np.median(phifeatures, axis=0).tolist()]
    # print(phifeatures)
    # print(len(phifeatures))
    # print(len(phifeatures[0]))
    
    # 写入目标文件中
    # print(respath)
    with open(respath, 'w') as output:
        for p in phifeatures:
            for i in range(len(p)):
                output.write(str(p[i]) + '\n')
    
# 提取AAC特征
def AAC(seqs):
    AA = 'ACDEFGHIKLMNPQRSTVWY*'
    #AA = 'ARNDCQEGHILKMFPSTWYV'
    encodings = []
    header = ['#']
    for i in AA:
        header.append(i)
    encodings.append(header)
    sequence = seqs
    
    #####核心代码#########
    count = Counter(sequence)
    for key in count:
        count[key] = float(count[key])/len(sequence)
    code = []
    for aa in AA:
        code.append(count[aa])
    #####################
    encodings.append(code)
    return code

Chemi_stats = {'A':{'C': 3, 'H': 7, 'O': 2, 'N': 1, 'S': 0},
                   'C':{'C': 3, 'H': 7, 'O': 2, 'N': 1, 'S': 1},
                   'D':{'C': 4, 'H': 7, 'O': 4, 'N': 1, 'S': 0},
                   'E':{'C': 5, 'H': 9, 'O': 4, 'N': 1, 'S': 0},
                   'F':{'C': 9, 'H': 11,'O': 2, 'N': 1, 'S': 0},
                   'G':{'C': 2, 'H': 5, 'O': 2, 'N': 1, 'S': 0},
                   'H':{'C': 6, 'H': 9, 'O': 2, 'N': 3, 'S': 0},
                   'I':{'C': 6, 'H': 13,'O': 2, 'N': 1, 'S': 0},
                   'K':{'C': 6, 'H': 14,'O': 2, 'N': 2, 'S': 0},
                   'L':{'C': 6, 'H': 13,'O': 2, 'N': 1, 'S': 0},
                   'M':{'C': 5, 'H': 11,'O': 2, 'N': 1, 'S': 1},
                   'N':{'C': 4, 'H': 8, 'O': 3, 'N': 2, 'S': 0},
                   'P':{'C': 5, 'H': 9, 'O': 2, 'N': 1, 'S': 0},
                   'Q':{'C': 5, 'H': 10,'O': 3, 'N': 2, 'S': 0},
                   'R':{'C': 6, 'H': 14,'O': 2, 'N': 4, 'S': 0},
                   'S':{'C': 3, 'H': 7, 'O': 3, 'N': 1, 'S': 0},
                   'T':{'C': 4, 'H': 9, 'O': 3, 'N': 1, 'S': 0},
                   'V':{'C': 5, 'H': 11,'O': 2, 'N': 1, 'S': 0},
                   'W':{'C': 11,'H': 12,'O': 2, 'N': 2, 'S': 0},
                   'Y':{'C': 9, 'H': 11,'O': 3, 'N': 1, 'S': 0}
                }

# 化学元素丰度
def chemicals(seqs):
    seq_new=seqs.replace('X','').replace('U','').replace('B','').replace('Z','')
    CE = 'CHONS'
    count = Counter(seq_new)
    code = []
    
    for c in CE:
        abundance_c = 0
        for key in count:
            num_c = Chemi_stats[key][c]
            abundance_c += num_c * count[key]
        code.append(abundance_c)
    return(code)

# 计算分子量
def mw(seq):
    seq_new=seq.replace('X','').replace('U','').replace('B','').replace('Z','')
    analysed_seq = ProteinAnalysis(seq_new)
    analysed_seq.monoisotopic = True
    mw = analysed_seq.molecular_weight()
    return([mw])
