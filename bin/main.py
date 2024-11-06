import os
import sys
import time
import logging
import argparse
import unittest

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
import matplotlib.pyplot as plt

from package.arg import get_parser
from package.datapreprocess import get_tfbs, get_seq, methylread_counter, flanking_bed
from package.backgroundprocess import read_bgprob_table, promoter, isfile, neighbor
from package.calculate import *
from package.figure_setting import *


def remove_dup_seqs(records):
    """移除重覆的序列記錄。"""
    checksums = set()
    for record in records:
        checksum = seguid(record.seq)
        if checksum in checksums:
            continue
        checksums.add(checksum)
        yield record


def main():
    start_time = time.time()
    parser = get_parser()
    parser.add_argument('--output_dir', type=str, default='./Output1', help='輸出目錄')
    args = parser.parse_args()

    # 調整日志配置
    logging.basicConfig(filename='methylseqlogo.log', level=logging.INFO,
                        format='%(asctime)s %(levelname)s:%(message)s')

    # 設置特定庫的日志級別
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)

    dir_path = os.path.dirname(os.path.realpath(__file__))

    # 提取命令行參數
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
    output_dir = args.output_dir

    pd.set_option('display.float_format', lambda x: '%.2f' % x)

    # 記錄繪圖信息
    if logotype in ['Shannon', 'Kullback-Liebler']:
        logging.info(f"繪制 {TF} {mode} {logotype} logo，物種：{species}，細胞類型：{celltype}，區域：{region} 背景。")
    elif logotype == 'riverlake':
        logging.info(f"繪制 {TF} {mode} {logotype}，物種：{species}，細胞類型：{celltype}。")
    else:
        logging.info(f"繪制 {TF} {mode} 所有 logos，物種：{species}，細胞類型：{celltype}。")

    # 確定跨度
    spanL = int(region) if region.isdigit() else 2
    spanR = spanL

    # 獲取 TFBS 和序列
    jaspar_file = args.jaspar
    remap_file = args.remap
    tfbs_bed, tfbs = get_tfbs(remap_file, jaspar_file, species, spanL, spanR)
    seq_file = get_seq(species, tfbs[0], tfbs[1], spanL, spanR)

    # 加載序列數據
    try:
        seqdata = SeqIO.to_dict(SeqIO.parse(seq_file, 'fasta'))
    except Exception as e:
        logging.warning("發現重覆序列，正在移除。")
        seqdata = remove_dup_seqs(SeqIO.parse(seq_file, 'fasta'))
        seqdata = SeqIO.to_dict(seqdata)

    seqdata = pd.DataFrame.from_dict(seqdata, orient='index').reset_index(drop=True)
    logging.info(f"序列數量：{len(seqdata)}")

    # 準備輸出目錄
    os.makedirs(output_dir, exist_ok=True)

    # 定義輸出文件名
    logoname = f"{TF}_{species}_{celltype}_{region}_{mode}_{logotype}_"
    ctx_file = os.path.join(output_dir, logoname + 'ctx.csv')
    cread_file = os.path.join(output_dir, logoname + 'cread.csv')
    tread_file = os.path.join(output_dir, logoname + 'tread.csv')

    # 加載或計算甲基化數據
    if os.path.isfile(ctx_file) and os.path.isfile(cread_file) and os.path.isfile(tread_file):
        ctxdata = pd.read_csv(ctx_file, sep='\t')
        creaddata = pd.read_csv(cread_file, sep='\t')
        treaddata = pd.read_csv(tread_file, sep='\t')
    else:
        ctxdata, creaddata, treaddata = methylread_counter(tfbs_bed, methylbed)
        # 保存數據以備將來使用
        ctxdata.to_csv(ctx_file, sep='\t', index=False)
        creaddata.to_csv(cread_file, sep='\t', index=False)
        treaddata.to_csv(tread_file, sep='\t', index=False)

    logging.info("甲基化數據處理完成。")

    # 計算基序長度
    motif_len = len(seqdata.columns) - (spanL + spanR)
    logging.info(f"{TF} 結合基序長度為 {motif_len} bp")

    # 調整繪圖長度
    if plotlen is None:
        plotlen = motif_len - beginningoftfbs + 1
    else:
        plotlen = int(plotlen)

    endpos = min(beginningoftfbs + plotlen - 1, motif_len)
    if (beginningoftfbs + plotlen - 1) > motif_len:
        logging.warning("用戶定義的繪圖長度超過了基序長度。")
    logging.info(f"繪制 {TF} 結合基序，從位置 {beginningoftfbs} 到位置 {endpos}")

    # 根據繪圖位置切片數據
    seqdata = seqdata.iloc[:, beginningoftfbs - 1 + spanL: beginningoftfbs - 1 + spanL + motif_len]
    ctxdata = ctxdata.iloc[:, beginningoftfbs - 1:endpos]
    creaddata = creaddata.iloc[:, beginningoftfbs - 1:endpos]
    treaddata = treaddata.iloc[:, beginningoftfbs - 1:endpos]

    # 檢查基因組文件是否存在
    genome_file = os.path.join(dir_path, 'genome', f'{species}.fa')
    isfile(genome_file)

    # 讀取背景概率
    if region == 'whole_genome':
        bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
    elif region.isdigit():
        path = os.path.join(dir_path, '..', 'Background_probability', 'neighbor',
                            f"{species}_{TF}_{celltype}_{region}_probability.txt")
        path_methyl = os.path.join(dir_path, '..', 'Background_probability', 'neighbor',
                                   f"{species}_{TF}_{celltype}_{region}_methyl_probability.txt")
        if os.path.isfile(path) and os.path.isfile(path_methyl):
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region, TF)
        else:
            flankingbed = flanking_bed(tfbs_bed, spanL, spanR)
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = neighbor(flankingbed, species, methylbed, celltype, region, TF)
    else:
        path = os.path.join(dir_path, '..', 'Background_probability', 'promoter',
                            f"{species}_{celltype}_{region}_probability.txt")
        path_methyl = os.path.join(dir_path, '..', 'Background_probability', 'promoter',
                                   f"{species}_{celltype}_{region}_methyl_probability.txt")
        if os.path.isfile(path) and os.path.isfile(path_methyl):
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
        else:
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = promoter(tfbs_bed, species, methylbed, celltype, region)

    # 執行計算
    ppm = calc_ppm(seqdata, TF, species, celltype, region)
    C_ratio, G_ratio, Cmethyls, Gmethyls, Freqs_ = calc_methylprob(ctxdata, creaddata, treaddata,
                                                                   bg_mCG, bg_mCHG, bg_mCHH, plotlen)

    # 繪圖
    if logotype in ['Kullback-Liebler', 'Shannon']:
        Cents = calc_methylation_entropy(C_ratio, G_ratio, Cmethyls, Gmethyls,
                                         bg_mCG, bg_mCHG, bg_mCHH, logotype)
        entropys = calc_totalEntropy(ppm, bgpps, Cents, logotype, plotlen, TF, species, celltype, region)
        four_base_heights = calc_perBaseEntropy(entropys, ppm, mode, TF, species, celltype, region)
        dippm = to4basedippm(seqdata, plotlen)
        dientropys, bg_dientropys_max, bg_dientropys_min = twomerBg(bgpps, dippm, plotlen)
        fig = set_fig(entropys, logotype, mode, plotlen)
        plotobj = seqLogoPlot(fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls, bgpps,
                              dientropys, bg_dientropys_max, bg_dientropys_min,
                              bg_mCG, bg_mCHG, bg_mCHH, Freqs_, mode, plotlen, threshold, TF)
        plotobj.plotlogo()
        logoname_png = f"{TF}_{species}_{celltype}_{region}_{mode}_{logotype}_seqlogo.png"
        plt.savefig(os.path.join(output_dir, logoname_png), bbox_inches='tight', dpi=600)
        logging.info(f"{logoname_png} 已保存到 {output_dir}。")
    elif logotype == 'riverlake':
        # 實現 riverlake 類型的繪圖
        pass  # 根據需要實現
    elif logotype == 'all':
        # 實現所有 logotype 的繪圖
        pass  # 根據需要實現

    total_time = time.time() - start_time
    minutes, seconds = divmod(total_time, 60)
    logging.info(f"總執行時間：{int(minutes)} 分 {int(seconds)} 秒")

    # 運行單元測試
    logging.info("正在運行測試以驗證輸出是否正確...")
    try:
        from package import test_methylseqlogo_outputs
    except ImportError:
        import test_methylseqlogo_outputs

    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(test_methylseqlogo_outputs)
    runner = unittest.TextTestRunner()
    result = runner.run(suite)

    if not result.wasSuccessful():
        logging.error("測試未通過，輸出存在錯誤。")
        sys.exit(1)
    else:
        logging.info("所有測試通過，輸出正確。")

    logging.info("")
    logging.info("--------------------------------------------")
    logging.info("")
if __name__ == '__main__':
    main()
