import os
import sys
import time
import unittest
import logging

import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

from package.arg import get_parser
from package.datapreprocess import get_tfbs, get_seq, methylread_counter, flanking_bed
from package.backgroundprocess import read_bgprob_table, promoter, isfile, neighbor
from package.calculate import *
from package.figure_setting import *


def remove_dup_seqs(records):
    """從 SeqRecord 叠代器中移除重覆的序列。"""
    checksums = set()
    for record in records:
        checksum = seguid(record.seq)
        if checksum in checksums:
            continue
        checksums.add(checksum)
        yield record


def main():
    filestart = time.time()

    # 設置日志記錄配置
    logging.basicConfig(
        filename="methylseqlogo.log",
        level=logging.INFO,  # 將日志級別設置為 INFO，忽略 DEBUG 信息
        format='%(asctime)s, %(name)s - %(levelname)s: %(message)s'
    )
    logger = logging.getLogger(__name__)

    # 設置特定庫的日志級別為 WARNING，避免它們輸出 DEBUG 信息
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

    logger.info("開始運行 methylseqlogo 腳本")

    if logotype in ['Shannon', 'Kullback-Liebler']:
        logger.info("繪制 %s %s %s logo，物種：%s，細胞類型：%s，區域：%s",
                    TF, mode, logotype, species, celltype, region)
    elif logotype == 'riverlake':
        logger.info("繪制 %s %s %s，物種：%s，細胞類型：%s", TF, mode, logotype, species, celltype)
    else:
        logger.info("繪制 %s %s 所有 logos，物種：%s，細胞類型：%s", TF, mode, species, celltype)

    # 確定跨度
    spanL = int(region) if region.isdigit() else 2
    spanR = spanL

    jaspar_file = args.jaspar
    remap_file = args.remap
    tfbs_bed, tfbs = get_tfbs(remap_file, jaspar_file, species, spanL, spanR)

    # 獲取序列
    seq_file = get_seq(species, tfbs[0], tfbs[1], spanL, spanR)

    # 解析序列
    try:
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, 'fasta'))
    except Exception as e:
        logger.warning("解析序列時出錯：%s。正在移除重覆序列。", e)
        seq_records = remove_dup_seqs(SeqIO.parse(seq_file, 'fasta'))
        seq_dict = SeqIO.to_dict(seq_records)

    seqdata = pd.DataFrame.from_dict(seq_dict, orient='index').reset_index(drop=True)
    logger.info("序列數量：%d", len(seqdata))

    start = time.time()

    # 輸出目錄
    output_dir = os.path.join(dir_path, '..', 'Output1')
    os.makedirs(output_dir, exist_ok=True)

    logoname = f"{TF}_{species}_{celltype}_{region}_{mode}_{logotype}_"

    ctx_file = os.path.join(output_dir, logoname + 'ctx.csv')
    cread_file = os.path.join(output_dir, logoname + 'cread.csv')
    tread_file = os.path.join(output_dir, logoname + 'tread.csv')

    logger.debug("上下文文件路徑：%s", ctx_file)

    # 讀取或生成甲基化數據
    if os.path.isfile(ctx_file) and os.path.isfile(cread_file) and os.path.isfile(tread_file):
        ctxdata = pd.read_csv(ctx_file, sep='\t')
        creaddata = pd.read_csv(cread_file, sep='\t')
        treaddata = pd.read_csv(tread_file, sep='\t')
    else:
        print("noData")
        logger.debug("無ctx,ctead,tread資料")
        ctxdata, creaddata, treaddata = methylread_counter(tfbs_bed, methylbed)

    end = time.time()
    logger.info("甲基化讀取計數完成，用時 %d 分鐘", (end - start) // 60)

    # 計算基序長度
    motif_len = len(seqdata.columns) - (spanL + spanR)
    logger.info("%s 結合基序長度為 %d bp", TF, motif_len)

    # 確定繪圖長度
    if plotlen is None:
        plotlen = motif_len - beginningoftfbs + 1

    # 確定結束位置
    if (beginningoftfbs + plotlen > motif_len):
        endpos = motif_len
        logger.warning("用戶定義的繪圖長度超過基序長度")
        logger.info("繪制 %s 結合基序，從位置 %d 到 %d", TF, beginningoftfbs, motif_len)
    else:
        endpos = beginningoftfbs + plotlen - 1
        logger.info("繪制 %s 結合基序，從位置 %d 到 %d", TF, beginningoftfbs, endpos)

    # 選擇相關數據
    seqdata = seqdata.iloc[:, beginningoftfbs - 1 + spanL: beginningoftfbs - 1 + spanL + motif_len]
    ctxdata = ctxdata.iloc[:, beginningoftfbs - 1:endpos]
    creaddata = creaddata.iloc[:, beginningoftfbs - 1:endpos]
    treaddata = treaddata.iloc[:, beginningoftfbs - 1:endpos]

    # 將數據保存為 CSV 文件
    ctxdata.to_csv(ctx_file, sep='\t', index=False)
    creaddata.to_csv(cread_file, sep='\t', index=False)
    treaddata.to_csv(tread_file, sep='\t', index=False)

    pseudocount = 1.0  # 設置偽計數（用於計算模塊）

    # 檢查物種基因組文件是否存在
    genome_file = os.path.join(dir_path, 'genome', f'{species}.fa')
    isfile(genome_file)

    # 讀取或生成背景概率
    if region == 'whole_genome':
        bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
    elif region.isdigit():
        logger.info('處理鄰近區域')
        bg_prob_file = os.path.join(dir_path, '..', 'Background_probability', 'neighbor',
                                    f"{species}_{TF}_{celltype}_{region}_probability.txt")
        bg_methyl_prob_file = os.path.join(dir_path, '..', 'Background_probability', 'neighbor',
                                           f"{species}_{TF}_{celltype}_{region}_methyl_probability.txt")
        if os.path.isfile(bg_prob_file) and os.path.isfile(bg_methyl_prob_file):
            logger.info('找到現有的鄰近背景概率')
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region, TF)
        else:
            logger.info('生成鄰近背景概率')
            flankingbed = flanking_bed(tfbs_bed, spanL, spanR)
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = neighbor(flankingbed, species, methylbed, celltype, region, TF)
    else:
        logger.info('處理啟動子區域')
        bg_prob_file = os.path.join(dir_path, '..', 'Background_probability', 'promoter',
                                    f"{species}_{celltype}_{region}_probability.txt")
        bg_methyl_prob_file = os.path.join(dir_path, '..', 'Background_probability', 'promoter',
                                           f"{species}_{celltype}_{region}_methyl_probability.txt")
        if os.path.isfile(bg_prob_file) and os.path.isfile(bg_methyl_prob_file):
            logger.info('找到現有的啟動子背景概率')
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
        else:
            logger.info('生成啟動子背景概率')
            start_time = time.time()
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = promoter(tfbs_bed, species, methylbed, celltype, region)
            end_time = time.time()
            logger.info('啟動子背景計算完成，用時 %d 分鐘', (end_time - start_time) // 60)

    # 計算位置概率矩陣
    ppm = calc_ppm(seqdata, TF, species, celltype, region)

    # 計算甲基化概率
    C_ratio, G_ratio, Cmethyls, Gmethyls, Freqs_ = calc_methylprob(
        ctxdata, creaddata, treaddata, bg_mCG, bg_mCHG, bg_mCHH, plotlen)

    # 根據 logotype 生成圖形
    if logotype in ['Kullback-Liebler', 'Shannon']:
        Cents = calc_methylation_entropy(
            C_ratio, G_ratio, Cmethyls, Gmethyls, bg_mCG, bg_mCHG, bg_mCHH, logotype)
        entropys = calc_totalEntropy(ppm, bgpps, Cents, logotype, plotlen, TF, species, celltype, region)
        four_base_heights = calc_perBaseEntropy(
            entropys, ppm, mode, TF, species, celltype, region)
        dippm = to4basedippm(seqdata, plotlen)
        dientropys, bg_dientropys_max, bg_dientropys_min = twomerBg(bgpps, dippm, plotlen)
        fig = set_fig(entropys, logotype, mode, plotlen)
        plotobj = seqLogoPlot(fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls,
                              bgpps, dientropys, bg_dientropys_max, bg_dientropys_min,
                              bg_mCG, bg_mCHG, bg_mCHH, Freqs_, mode, plotlen, threshold, TF)
        plotobj.plotlogo()
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
            plotobj = seqLogoPlot(fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls,
                                  bgpps, bg_mCG, bg_mCHG, bg_mCHH, Freqs_)
            plotobj.plotlogo()
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
    logger.info("總耗時：%d 分 %d 秒", int(minutes), seconds)

    print("正在運行測試以驗證輸出是否正確...")

    # 運行測試以驗證輸出
    try:
        from package import test_methylseqlogo_outputs
    except ImportError:
        import test_methylseqlogo_outputs

    # 創建測試加載器和測試套件
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(test_methylseqlogo_outputs)

    # 運行測試
    runner = unittest.TextTestRunner()
    result = runner.run(suite)

    # 檢查測試結果
    if not result.wasSuccessful():
        print("測試未通過，輸出存在錯誤。")
        logger.info("測試未通過，輸出存在錯誤。")
        sys.exit(1)  # 退出程序，返回錯誤碼
    else:
        print("所有測試通過，輸出正確。")
        logger.info("所有測試通過，輸出正確。")


if __name__ == '__main__':
    main()
