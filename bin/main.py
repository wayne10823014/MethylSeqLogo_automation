import os
import sys
import time
import unittest
import logging

import pandas as pd
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.SeqUtils.CheckSum import seguid

from package.arg import get_parser
from package.datapreprocess import get_tfbs, get_seq, methylread_counter, flanking_bed
from package.backgroundprocess import read_bgprob_table, promoter, isfile, neighbor
from package.calculate import *
from package.figure_setting import *


def remove_dup_seqs(seq_dict):
    """从序列字典中移除重复的序列。"""
    checksums = set()
    unique_seq_dict = {}
    for key, sequence in seq_dict.items():
        checksum = seguid(Seq(sequence))
        if checksum in checksums:
            continue
        checksums.add(checksum)
        unique_seq_dict[key] = sequence
    return unique_seq_dict


def main():
    filestart = time.time()

    # 设置日志记录配置
    logging.basicConfig(
        filename="methylseqlogo.log",
        level=logging.INFO,  # 将日志级别设置为 INFO，忽略 DEBUG 信息
        format='%(asctime)s, %(name)s - %(levelname)s: %(message)s'
    )
    logger = logging.getLogger(__name__)

    # 设置特定库的日志级别为 WARNING，避免它们输出 DEBUG 信息
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    logging.getLogger('PIL.PngImagePlugin').setLevel(logging.WARNING)

    parser = get_parser()
    args = parser.parse_args()

    dir_path = os.path.dirname(os.path.realpath(__file__))

    TF = args.transcriptionfactor
    species = args.species
    celltype = args.celltypes
    mode = args.mode
    logotype = args.logotype
    region = args.regionofinterest

    threshold = float(args.threshold)
    beginningoftfbs = args.beginningoftfbs
    plotlen = args.plotlen

    methylbed = [args.mcg, args.mchg, args.mchh]

    pd.set_option('display.float_format', lambda x: '%.2f' % x)

    logger.info("开始运行 methylseqlogo 脚本")

    if logotype in ['Shannon', 'Kullback-Liebler']:
        logger.info("绘制 %s %s %s logo，物种：%s，细胞类型：%s，区域：%s",
                    TF, mode, logotype, species, celltype, region)
    elif logotype == 'riverlake':
        logger.info("绘制 %s %s %s，物种：%s，细胞类型：%s", TF, mode, logotype, species, celltype)
    else:
        logger.info("绘制 %s %s 所有 logos，物种：%s，细胞类型：%s", TF, mode, species, celltype)

    # 确定跨度
    spanL = int(region) if region.isdigit() else 2
    spanR = spanL

    jaspar_file = args.jaspar
    remap_file = args.remap
    tfbs_bed, modify_coordinate_list, modify_coordinate, motif_len = get_tfbs(
        remap_file, jaspar_file, species, spanL, spanR)

    # 获取序列并构建 total 和 seq_dict
    total, seq_dict = get_seq(species, modify_coordinate_list, modify_coordinate, spanL, spanR, motif_len)

    # 解析序列
    try:
        # 将序列字典转换为 DataFrame，每个序列作为一行
        seqdata = pd.DataFrame([list(seq.upper()) for seq in seq_dict.values()])
    except Exception as e:
        logger.warning("解析序列时出错：%s。正在移除重复序列。", e)
        seq_dict = remove_dup_seqs(seq_dict)
        seqdata = pd.DataFrame([list(seq.upper()) for seq in seq_dict.values()])

    logger.info("序列数量：%d", len(seqdata))

    start = time.time()

    # 输出目录
    output_dir = os.path.join(dir_path, '..', 'Output1')
    os.makedirs(output_dir, exist_ok=True)

    logoname = f"{TF}_{species}_{celltype}_{region}_{mode}_{logotype}_"

    ctx_file = os.path.join(output_dir, logoname + 'ctx.csv')
    cread_file = os.path.join(output_dir, logoname + 'cread.csv')
    tread_file = os.path.join(output_dir, logoname + 'tread.csv')

    logger.debug("上下文文件路径：%s", ctx_file)

    # 读取或生成甲基化数据
    # if os.path.isfile(ctx_file) and os.path.isfile(cread_file) and os.path.isfile(tread_file):
    #     ctxdata = pd.read_csv(ctx_file, sep='\t')
    #     creaddata = pd.read_csv(cread_file, sep='\t')
    #     treaddata = pd.read_csv(tread_file, sep='\t')
    # else:
    logger.debug("无 ctx, cread, tread 数据，开始计算...")
    ctxdata, creaddata, treaddata = methylread_counter(tfbs_bed, methylbed, total)

        # 将数据保存为 CSV 文件
    ctxdata.to_csv(ctx_file, sep='\t', index=False)
    creaddata.to_csv(cread_file, sep='\t', index=False)
    treaddata.to_csv(tread_file, sep='\t', index=False)

    end = time.time()
    logger.info("甲基化读取计数完成，用时 %d 分钟", (end - start) // 60)

    # 计算基序长度
    motif_len_seqdata = seqdata.shape[1] - (spanL + spanR)
    logger.info("%s 结合基序长度为 %d bp", TF, motif_len_seqdata)

    # 确定绘图长度
    if plotlen is None:
        plotlen = motif_len_seqdata - beginningoftfbs + 1

    # 确定结束位置
    if (beginningoftfbs + plotlen > motif_len_seqdata):
        endpos = motif_len_seqdata
        logger.warning("用户定义的绘图长度超过基序长度")
        logger.info("绘制 %s 结合基序，从位置 %d 到 %d", TF, beginningoftfbs, motif_len_seqdata)
    else:
        endpos = beginningoftfbs + plotlen - 1
        logger.info("绘制 %s 结合基序，从位置 %d 到 %d", TF, beginningoftfbs, endpos)

    # 选择相关数据
    seqdata = seqdata.iloc[:, beginningoftfbs - 1 + spanL: beginningoftfbs - 1 + spanL + motif_len_seqdata]
    ctxdata = ctxdata.iloc[:, beginningoftfbs - 1:endpos]
    creaddata = creaddata.iloc[:, beginningoftfbs - 1:endpos]
    treaddata = treaddata.iloc[:, beginningoftfbs - 1:endpos]

    pseudocount = 1.0  # 设置伪计数（用于计算模块）

    # 检查物种基因组文件是否存在
    genome_file = os.path.join(dir_path, 'genome', f'{species}.fa')
    isfile(genome_file)

    # 读取或生成背景概率
    if region == 'whole_genome':
        bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
    elif region.isdigit():
        logger.info('处理邻近区域')
        bg_prob_file = os.path.join(dir_path, '..', 'Background_probability', 'neighbor',
                                    f"{species}_{TF}_{celltype}_{region}_probability.txt")
        bg_methyl_prob_file = os.path.join(dir_path, '..', 'Background_probability', 'neighbor',
                                           f"{species}_{TF}_{celltype}_{region}_methyl_probability.txt")
        if os.path.isfile(bg_prob_file) and os.path.isfile(bg_methyl_prob_file):
            logger.info('找到现有的邻近背景概率')
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region, TF)
        else:
            logger.info('生成邻近背景概率')
            flankingbed = flanking_bed(tfbs_bed, spanL, spanR)
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = neighbor(flankingbed, species, methylbed, celltype, region, TF)
    else:
        logger.info('处理启动子区域')
        bg_prob_file = os.path.join(dir_path, '..', 'Background_probability', 'promoter',
                                    f"{species}_{celltype}_{region}_probability.txt")
        bg_methyl_prob_file = os.path.join(dir_path, '..', 'Background_probability', 'promoter',
                                           f"{species}_{celltype}_{region}_methyl_probability.txt")
        if os.path.isfile(bg_prob_file) and os.path.isfile(bg_methyl_prob_file):
            logger.info('找到现有的启动子背景概率')
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
        else:
            logger.info('生成启动子背景概率')
            start_time = time.time()
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = promoter(tfbs_bed, species, methylbed, celltype, region)
            end_time = time.time()
            logger.info('启动子背景计算完成，用时 %d 分钟', (end_time - start_time) // 60)

    # 计算位置概率矩阵
    ppm = calc_ppm(seqdata, TF, species, celltype, region)

    # 计算甲基化概率
    C_ratio, G_ratio, Cmethyls, Gmethyls, Freqs_ = calc_methylprob(
        ctxdata, creaddata, treaddata, bg_mCG, bg_mCHG, bg_mCHH, plotlen)

    # 根据 logotype 生成图形
    if logotype in ['Kullback-Liebler', 'Shannon']:
        Cents = calc_methylation_entropy(
            C_ratio, G_ratio, Cmethyls, Gmethyls, bg_mCG, bg_mCHG, bg_mCHH, logotype)
        entropys = calc_totalEntropy(ppm, bgpps, Cents, logotype, plotlen, TF, species, celltype, region)
        four_base_heights = calc_perBaseEntropy(
            entropys, ppm, mode, TF, species, celltype, region)
        dippm = to4basedippm(seqdata, plotlen)
        dientropys, bg_dientropys_max, bg_dientropys_min = twomerBg(bgpps, dippm, plotlen)
        fig = set_fig(entropys, logotype, mode, plotlen)
        plotobj = SeqLogoPlot(fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls,
                              bgpps, dientropys, bg_dientropys_max, bg_dientropys_min,
                              bg_mCG, bg_mCHG, bg_mCHH, Freqs_, mode, plotlen, threshold, TF)
        plotobj.plot_logo()
        logoname_file = f"{logoname}seqlogo.png"
        output_filepath = os.path.join(output_dir, logoname_file)
        plt.savefig(output_filepath, bbox_inches='tight', dpi=600)
        logger.info('%s 已保存至 %s', logoname_file, output_dir)
    elif logotype == 'riverlake':
        dippm = to4basedippm(seqdata)
        Cents = calc_methylation_entropy(
            C_ratio, G_ratio, bg_mCG, bg_mCHG, bg_mCHH, logotype)
        entropys = calc_totalEntropy(ppm, bgpps, Cents, logotype)
        four_base_heights = calc_perBaseEntropy(entropys, ppm)
        fig = plt.figure(figsize=(plotlen + 1, 3.0))
        riverlakeobj = riverLake(fig, celltype, ppm, dippm, Cmethyls, Gmethyls,
                                 bgpps, bg_mCG, bg_mCHG, bg_mCHH, Freqs_)
        riverlakeobj.plotRiverLake()
        logoname_file = f"{logoname}seqlogo_bar7.pdf"
        output_filepath = os.path.join(output_dir, logoname_file)
        plt.savefig(output_filepath, bbox_inches='tight', dpi=600)
        logger.info('%s 已保存至 %s', logoname_file, output_dir)
    elif logotype == 'all':
        for lt in ['Kullback-Liebler', 'Shannon']:
            Cents = calc_methylation_entropy(
                C_ratio, G_ratio, bg_mCG, bg_mCHG, bg_mCHH, lt)
            entropys = calc_totalEntropy(ppm, bgpps, Cents, lt)
            four_base_heights = calc_perBaseEntropy(entropys, ppm)
            fig = set_fig(entropys)
            plotobj = SeqLogoPlot(fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls,
                                  bgpps, bg_mCG, bg_mCHG, bg_mCHH, Freqs_)
            plotobj.plot_logo()
            logoname_file = f"{logoname}{lt}_seqlogo_bar7.pdf"
            output_filepath = os.path.join(output_dir, logoname_file)
            plt.savefig(output_filepath, bbox_inches='tight', dpi=600)
            logger.info('%s 已保存至 %s', logoname_file, output_dir)
        dippm = to4basedippm(seqdata)
        dientropys = twomerBg(bgpps, dippm)
        Cents = calc_methylation_entropy(
            C_ratio, G_ratio, bg_mCG, bg_mCHG, bg_mCHH, 'riverlake')
        entropys = calc_totalEntropy(ppm, bgpps, Cents, 'riverlake')
        four_base_heights = calc_perBaseEntropy(entropys, ppm)
        fig = plt.figure(figsize=(plotlen, 3.0))
        riverlakeobj = riverLake(fig, celltype, ppm, dippm, Cmethyls, Gmethyls,
                                 bgpps, bg_mCG, bg_mCHG, bg_mCHH, Freqs_)
        riverlakeobj.plotRiverLake()
        logoname_file = f"{logoname}riverlake_seqlogo_bar7.png"
        output_filepath = os.path.join(output_dir, logoname_file)
        plt.savefig(output_filepath, bbox_inches='tight', dpi=600)
        logger.info('%s 已保存至 %s', logoname_file, output_dir)

    total_time = time.time() - filestart
    minutes = total_time // 60
    seconds = round(total_time % 60)
    logger.info("总耗时：%d 分 %d 秒", int(minutes), seconds)

    print("正在运行测试以验证输出是否正确...")

    # 运行测试以验证输出
    try:
        from package import test_methylseqlogo_outputs
    except ImportError:
        import test_methylseqlogo_outputs

    # 创建测试加载器和测试套件
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(test_methylseqlogo_outputs)

    # 运行测试
    runner = unittest.TextTestRunner()
    result = runner.run(suite)

    # 检查测试结果
    if not result.wasSuccessful():
        print("测试未通过，输出存在错误。")
        logger.info("测试未通过，输出存在错误。")
        sys.exit(1)  # 退出程序，返回错误码
    else:
        print("所有测试通过，输出正确。")
        logger.info("所有测试通过，输出正确。")
    logger.info("")
    logger.info("--------------------------------------------")
    logger.info("")


if __name__ == '__main__':
    main()
