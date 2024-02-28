# ---------------------------------特征转换时的参数----------------------------
# predphi特征提取的维度参数
predphi_AAC_nums = 21
predphi_che_nums = 5
predphi_mw_nums = 1
predphi_type_nums = 6

# 模型参数
input_shape = [2,6,27]
act = 'relu'
fact = 'softmax'
classes = 2
dropout = 0.5

base_processed = 'PredPHI/our_predphi/prodatasets/predphi/'
# ----------------------------------第一批数据生成文件路径-----------------------------
base_processed1 = base_processed + 'data1/'

# 经过prephi特征提取后的细菌文件夹
fed_bac_predphi_protein1 = base_processed1 + 'fed_bac_predphi_protein1'

# 经过predphi特征提取后的噬菌体文件夹
fed_phage_predphi_protein1 = base_processed1 + 'fed_phage_predphi_protein1'

# 经过predphi特征提取后的混合324维特征序列文件夹
fed_mix_predphi_protein1 = base_processed1 + 'fed_mix_predphi_protein1'

# predphi蛋白质序列最大值最小值文件夹
fed_minmax_predphi_protein1 = base_processed1 + 'fed_minmax_predphi_protein1'

# predphi蛋白质序列最小值文件
fed_min_predphi_protein1 = fed_minmax_predphi_protein1 + '/min.csv'

# predphi蛋白质序列最大值文件
fed_max_predphi_protein1 = fed_minmax_predphi_protein1 + '/max.csv'

# predphi蛋白质序列特征归一化后的文件夹
fed_norm_predphi_protein1 = base_processed1 + 'fed_norm_predphi_protein1'

# prephi包含所有数据的csv文件
fed_all_predphi_protein1 = base_processed1 + 'fed_all_predphi_protein1.csv'

# ----------------------------------第二批数据生成文件路径-----------------------------
base_processed2 = base_processed + 'data2/'

# 经过prephi特征提取后的细菌文件夹
fed_bac_predphi_protein2 = base_processed2 + 'fed_bac_predphi_protein2'

# 经过predphi特征提取后的噬菌体文件夹
fed_phage_predphi_protein2 = base_processed2 + 'fed_phage_predphi_protein2'

# 经过predphi特征提取后的混合324维特征序列文件夹
fed_mix_predphi_protein2 = base_processed2 + 'fed_mix_predphi_protein2'

# predphi蛋白质序列最大值最小值文件夹
fed_minmax_predphi_protein2 = base_processed2 + 'fed_minmax_predphi_protein2'

# predphi蛋白质序列最小值文件
fed_min_predphi_protein2 = fed_minmax_predphi_protein2 + '/min.csv'

# predphi蛋白质序列最大值文件
fed_max_predphi_protein2 = fed_minmax_predphi_protein2 + '/max.csv'

# predphi蛋白质序列特征归一化后的文件夹
fed_norm_predphi_protein2 = base_processed2 + 'fed_norm_predphi_protein2'

# prephi包含所有数据的csv文件
fed_all_predphi_protein2 = base_processed2 + 'fed_all_predphi_protein2.csv'

# ----------------------------------结果文件路径-----------------------------
base_result = 'PredPHI/our_predphi/result/'

# predphi生成的模型文件
predphi_model_out = base_result + 'prephi_model.h5'