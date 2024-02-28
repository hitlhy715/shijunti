# -----------------------------------原始数据路径------------------------------
# 原始数据根目录
base_datasets = './datasets/'

# -----------------------------------第一批------------------------------
# 第一批数据根目录
base_datasets1 = base_datasets + 'datasets1/'

# 包含第一批数据全部作用关系的表格
score_file = base_datasets1 + 'score.csv'

# 噬菌体fasta数据路径          ####第一批数据
phage_fasta = base_datasets1 + 'phage_fasta'

# 噬菌体vcf数据路径            ####第一批数据
phage_vcf = base_datasets1 + 'phage_vcf'

# 宿主fasta数据路径            ####第一批数据
bac_fasta = base_datasets1 + 'bac_fasta'

# 噬菌体vcf数据路径            ####第一批数据
bac_vcf = base_datasets1 + 'bac_vcf'

# 噬菌体氨基酸数据路径          ####第一批数据
bacpro_faa = base_datasets1 + 'bacprotein_faa'

# 宿主氨基酸数据路径            ####第一批数据
phagepro_faa = base_datasets1 + 'phageprotein_faa'

# -----------------------------------第二批------------------------------

# 第二批数据根目录
base_datasets2 = base_datasets + 'datasets2/'

# 第二批噬菌体核苷酸数据路径
phage_fasta2 = base_datasets2 + 'phage_fasta2'

# 第二批噬菌体氨基酸数据路径
phagepro_faa2 = base_datasets2 + 'phagepro_faa2'

# 第二批宿主核苷酸数据路径
bac_fasta2 = base_datasets2 + 'bac_fasta2'

#第二批宿主氨基酸数据路径
bacpro_faa2 = base_datasets2 + 'bacpro_faa2'

# ----------------------------------生成文件路径-----------------------------
base_processed = './prodatasets/'

# 包含第一批数据全部13029个相互作用关系的txt文件名
allinteract_file1 = base_processed + 'allInteraction1.txt'

# 包含第二批数据全部  个相互作用关系的txt文件名
allinteract_file2 = base_processed + 'allInteraction2.txt'