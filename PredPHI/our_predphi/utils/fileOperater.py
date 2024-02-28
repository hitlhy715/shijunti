


# 单个fasta文件的特征提取
def readTxtInfo(filepath):
    f = open(filepath,'r')
    dlist = []
    idx = -1
    for line in f:
        if line.startswith('>'):
            dlist.append([''])
            idx += 1
        else:
            # print(len(line))
            dlist[idx] += line.split('\n')[0]
    f.close()
       
    # print(len(dlist))
    # print(dlist[0])
    # print(''.join(dlist[0]))
    
    return [''.join(d) for d in dlist]