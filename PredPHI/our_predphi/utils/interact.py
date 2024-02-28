# 从csv矩阵文件，生成txt格式的数据
import pandas as pd
import os
import numpy as np

def deal():
    datasets = './datasets/datasets2/phagepro_faa2'
    for item in os.listdir(datasets):
        old_path = os.path.join(datasets,item)
        if item.split('.')[1] == 'fasta_genome':
            new_name = os.path.join(datasets, item.split('.')[0] + '.faa')
            os.rename(old_path, new_name)
        elif item.split('.')[2] == 'fasta_genome':
            new_name = os.path.join(datasets, item.split('.')[0] + '.' + item.split('.')[1] + '.faa')
            os.rename(old_path, new_name)
            
def phagesmatch(bacs):
    path = './datasets/datasets2/phagepro_faa2'
    
    ds = []
    for i in os.listdir(path):
        if len(i.split('.')) == 2:
            ds.append(i.split('.')[0])
        if len(i.split('.')) == 3:
            ds.append(i.split('.')[0] + '.' + i.split('.')[1])
    ds = np.array(ds)
    
    new_bacs = []
    for item in bacs:
        ls = item.split(' ')
        new_item = '_'.join(ls)
        new_bacs.append(new_item)
    new_bacs = np.array(new_bacs)
    
    print(ds.shape)
    print(new_bacs.shape)
    
    print(ds)
    print(new_bacs)
    
    print(set(new_bacs) - set(ds))
    print(set(ds) - set(new_bacs))
    
    print(len(set(new_bacs) - set(ds)))
    print(len(set(ds) - set(new_bacs)))
    
def bacmatch(bacs):
    path = './datasets/datasets2/bacpro_faa2'
    ds = []
    for i in os.listdir(path):
        if len(i.split('.')) == 2:
            ds.append(i.split('.')[0])
        if len(i.split('.')) == 3:
            ds.append(i.split('.')[0] + '.' + i.split('.')[1])
    ds = np.array(ds)
    
    new_bacs = []
    for item in bacs:
        ls = item.split(' ')
        new_item = '_'.join(ls)
        new_bacs.append(new_item)
    new_bacs = np.array(new_bacs)
    
    print(ds.shape)
    print(new_bacs.shape)
    
    print(ds)
    print(new_bacs)
    
    print(set(new_bacs) - set(ds))
    print(set(ds) - set(new_bacs))
    
    print(len(set(new_bacs) - set(ds)))
    print(len(set(ds) - set(new_bacs)))


def generate(scores):
    res = []
    pnums = 0
    nnums = 0
    
    #行数
    row = len(scores.iloc[:,0])
    print(row)
    # 列数
    col = len(scores.columns.values)-1
    print(col)

    # 列名
    cols = scores.columns.values
    print(cols)
    
    for i in range(row):
        for j in range(col):
            row_name = scores.iloc[i,0]
            col_name = cols[j+1]
            res.append([row_name, col_name, int(scores.iloc[i, j+1])])
            if scores.iloc[i, j+1] == 1:
                pnums += 1  # 正样本
            elif scores.iloc[i, j+1] == 0:   
                nnums += 1  # 负样本
            
    print(len(res))
    print(pnums)
    print(nnums)
    
        # 将res写到txt文件中
    txt_file = './prodatasets/allInteraction2.txt'

    with open(txt_file, 'w') as al:
        for line in res:
            if len(line[0].split('\t'))>1 or len(line[1].split('\t'))>1:
                print('wrong')
            info = str(line[0])+'\t'+line[1]+'\t'+str(line[2])+'\n'
            al.write(info)


if __name__ == '__main__':
    score = './datasets/datasets2/score2.csv'

    scores = pd.read_csv(score, encoding='gbk')
    # print(scores.shape)
    # print(scores.head(3))

    # print(scores.iloc[:,0])
    # phagesmatch(np.array(scores.iloc[:,0]))

    # print(scores.columns.values[1:])
    # print(len(scores.columns.values[1:]))
    # bacmatch(scores.columns.values[1:])

    # 已经完成名字匹配，全都没问题
    generate(scores)