# 目录结构： config.py 一些目录名配置
#           datasets  原始数据
#           predphi   模型代码
#               original code 原始论文作者提供的代码
#               our predphi   修正后的代码，核心代码
#                   prodataset   中间数据，基本被删了，可以通过代码生成
#                   result       训出来的模型文件
#                   utils        功能函数，文件操作和一些类型转换
#                       DNA2Animo    dna转氨基酸
#                       fileoperater  文件操作函数
#                       infoStatistic   数据统计函数
#                       interact      数据预处理函数
#                   fe_predphi_protein 按照论文中的方式，进行氨基酸特征提取
#                   keras model  模型
#                   main 运行文件，训练和推理时可以直接用main运行，但是需要先把中间文件都按次序生成
#                   obtain feature 完整的特征提取，用来生成中间文件
#                   pconfig 模型独特的一些配置项，训练参数，中间文件路径等，主要是考虑存在其他模型时
#                   train keras 模型的训练和推理的类


# 执行顺序
# result里面是训好的，可以直接推理用
# 中间数据生成，obtainfeature中
#   主函数if __name__ == '__main__':中
#   先定义类，然后下面的函数按顺序每次执行一个即可

# 模型训练和推理，main文件中
# model.construct_model() 为训练
# model.test(fed_all_predphi_protein2) 为推理