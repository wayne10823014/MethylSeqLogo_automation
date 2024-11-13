import argparse

def get_parser():
    parser = argparse.ArgumentParser(description='MethylSeqLogo')
    
    # 必需的文件输入
    input_file_group = parser.add_argument_group('必需的文件输入')
    input_file_group.add_argument('-J', '--jaspar', required=True, help='JASPAR 文件路径')
    input_file_group.add_argument('-R', '--remap', required=True, help='ReMap 文件路径')
    input_file_group.add_argument('-M', '--methylationinfo', nargs='+', type=str, help='WGBS 在 CpG 位点的甲基化信息')
    input_file_group.add_argument('-CG', '--mcg', required=True, nargs='+', type=str, help='WGBS 在 CpG 位点的甲基化信息')
    input_file_group.add_argument('-CHG', '--mchg', required=True, nargs='+', type=str, help='WGBS 在 CHG 位点的甲基化信息')
    input_file_group.add_argument('-CHH', '--mchh', required=True, nargs='+', type=str, help='WGBS 在 CHH 位点的甲基化信息')
    
    # 与 TFBS 计算相关的参数
    tfbs_group = parser.add_argument_group('TFBS 计算参数')
    tfbs_group.add_argument('-TF', '--transcriptionfactor', metavar='TF', default='ZBTB33', help='转录因子名称')
    tfbs_group.add_argument('-S', '--species', choices=['human', 'mouse', 'arabidopsis', 'maize'], default='human', metavar='species', help='输入数据的物种')
    tfbs_group.add_argument('-C', '--celltypes', choices=['H1-hESC', 'HepG2', 'K562', 'GM12878', 'HeLa-S3', 'naive-hESC', 'primed-hESC', 'mESC', 'mMPGC_E13p5', 'mMPGC_E16p5', 'mFPGC_E13p5', 'mFPGC_E16p5', 'Leaf', 'Inflorescence', 'Root', 'shoot'], default='HepG2', metavar='cell types', help='细胞类型')
    tfbs_group.add_argument('-T', '--threshold', type=float, default=0.4, help='显示甲基化概率的阈值')
    tfbs_group.add_argument('-ROI', '--regionofinterest', default='whole_genome', help='感兴趣的基因组区域')
    # tfbs_group.add_argument('--bg', help='背景模型文件')  # 如有需要，可取消注释
    
    # 与 MethylSeqLogo 显示相关的参数
    display_group = parser.add_argument_group('MethylSeqLogo 显示参数')
    display_group.add_argument('-bot', '--beginningoftfbs', type=int, default=1, help='绘制 Logo 的 TFBS 起始位置')
    display_group.add_argument('-plen', '--plotlen', type=int, help='绘制 Logo 的 TFBS 长度')
    display_group.add_argument('-m', '--mode', choices=['Methyl', 'Regular'], default='Methyl', help='是否包含甲基化信息')
    display_group.add_argument('-l', '--logotype', choices=['Shannon', 'Kullback-Liebler', 'riverlake', 'all'], default='Kullback-Liebler', help='Logo 类型')
    
    return parser
