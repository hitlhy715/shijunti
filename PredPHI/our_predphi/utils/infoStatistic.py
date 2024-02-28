from collections import Counter


# 获取当前文件所有出现过的氨基酸种类
# 输入：数组格式，每个元素是一个node的字符串
def getProteinKinds(seqs):
    kinds = []
    for s in seqs:
        kind = list(Counter(s).keys())
        kinds.extend(kind)
    kinds = list(set(kinds))
    # print(kinds)
    # print(len(kinds))
    return ''.join(kinds)