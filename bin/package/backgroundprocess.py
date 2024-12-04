import os
import time
import itertools
from functools import reduce
from multiprocessing import Pool, cpu_count

import pandas as pd
import numpy as np
import pybedtools

from package.datapreprocess import modify_bed

# 獲取當前腳本的目錄路徑
dir_path = os.path.dirname(os.path.realpath(__file__))

# 定義輸出目錄路徑
output_path = os.path.join(dir_path, "../../Output1/")

# 獲取CPU核心數
CPU_COUNT = cpu_count()


def isfile(path, flag=1):
    """
    檢查文件是否存在，如果存在則讀取為 DataFrame。

    參數:
    - path (str): 文件路徑。
    - flag (int): 標志位，默認為1。當 flag 為1時，如果文件不存在則打印警告信息。

    返回:
    - DataFrame，如果文件存在。
    - None，如果文件不存在。
    """
    pd.set_option('display.float_format', lambda x: '%.4f' % x)
    if flag:
        if not os.path.isfile(path):
            print(f"{path} does not exist. Job cancelled.")
            print(f"{path} 檔案不存在，工作終止。")
            return None
    if os.path.isfile(path):
        data = pd.read_table(path, sep="\t", header=0, index_col=0)
        return data
    else:
        print(f"{path} does not exist. Job cancelled.")
        print(f"{path} 檔案不存在，工作終止。")
        return None


def read_bgprob_table(species, celltype, region, TF=None):
    """
    讀取背景概率表。

    返回:
    - bgpps (DataFrame): 背景概率。
    - bg_mCG (float): 背景甲基化概率（mCG）。
    - bg_mCHG (float): 背景甲基化概率（mCHG）。
    - bg_mCHH (float): 背景甲基化概率（mCHH）。
    """
    if region == 'whole_genome':
        species_file = os.path.join(dir_path, "../../Background_probability", region,
                                    f"{species}_{region}_probability.txt")
        celltype_file = os.path.join(dir_path, "../../Background_probability", region,
                                     f"{species}_{celltype}_{region}_methyl_probability.txt")
    elif region == 'promoter':
        species_file = os.path.join(dir_path, "../../Background_probability", "promoter",
                                    f"{species}_{celltype}_{region}_probability.txt")
        celltype_file = os.path.join(dir_path, "../../Background_probability", "promoter",
                                     f"{species}_{celltype}_{region}_methyl_probability.txt")
    else : 
        species_file = os.path.join(dir_path, "../../Background_probability", "neighbor",
                                    f"{species}_{TF}_{celltype}_{region}_probability.txt")
        celltype_file = os.path.join(dir_path, "../../Background_probability", "neighbor",
                                     f"{species}_{TF}_{celltype}_{region}_methyl_probability.txt")

    species_data = isfile(species_file)
    celltype_data = isfile(celltype_file)

    if species_data is None or celltype_data is None:
        print("Error reading background probability tables.")
        return None, None, None, None

    region_key = 'neighbor' if region.isdigit() else region

    probmatrix_species = species_data[region_key].astype('float64')
    bgpps = probmatrix_species[['A', 'C', 'G', 'T', 'CpG', 'CHG', 'CHH',
                                'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
                                'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']]

    probmatrix_celltype = celltype_data[region_key].astype('float64')
    bg_mCG = probmatrix_celltype['mCG'].astype('float64')
    bg_mCHG = probmatrix_celltype['mCHG'].astype('float64')
    bg_mCHH = probmatrix_celltype['mCHH'].astype('float64')

    print("Background methylation probabilities: ")
    print(probmatrix_celltype[['mCG', 'mCHG', 'mCHH']])
    print("\n")

    return bgpps, bg_mCG, bg_mCHG, bg_mCHH


def count_nucleotides_and_motifs(seq):
    """
    統計序列中的單核苷酸和二核苷酸頻率，並統計CpG、CHG和CHH基序的出現次數。

    參數:
    - seq (str): DNA序列。

    返回:
    - 列表，包含以下內容：
        - cpg (int): CpG位點數量。
        - chg (int): CHG位點數量。
        - chh (int): CHH位點數量。
        - mononucleotide_counts (dict): A、C、G、T的計數。
        - dinucleotide_counts (dict): 二核苷酸（AA、AC、...、TT）的計數。
    """
    if seq.startswith('>') or (seq.startswith('N') and len(set(seq)) == 2):
        return 0

    mononucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    dinucleotide_counts = {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0,
                           'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0,
                           'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0,
                           'TA': 0, 'TC': 0, 'TG': 0, 'TT': 0}
    cpg, chg, chh = 0, 0, 0

    seq = seq.upper()
    seq_length = len(seq)

    for i in range(seq_length):
        base = seq[i]
        if base in "ACGT":
            mononucleotide_counts[base] += 1

            if i < seq_length - 1 and seq[i + 1] in "ACGT":
                dinucleotide = base + seq[i + 1]
                if dinucleotide in dinucleotide_counts:
                    dinucleotide_counts[dinucleotide] += 1

        # 計算CpG、CHG、CHH基序
        if 2 <= i <= seq_length - 3:
            if base == 'C':
                if seq[i + 1] == 'G':
                    cpg += 1
                elif seq[i + 2] == 'G':
                    chg += 1
                else:
                    chh += 1
            elif base == 'G':
                if seq[i - 1] == 'C':
                    cpg += 1
                elif seq[i - 2] == 'C':
                    chg += 1
                else:
                    chh += 1

    return [cpg, chg, chh, mononucleotide_counts, dinucleotide_counts]


def process_sequences(file_sequences, species, celltype, region, directory, TF):
    """
    處理序列，計算核苷酸和基序計數，並將概率保存到文件。

    參數:
    - file_sequences (list): 序列列表。
    - species (str): 物種名稱。
    - celltype (str): 細胞類型。
    - region (str): 區域。
    - directory (str): 保存輸出文件的目錄。
    - TF (str): 轉錄因子名稱。
    """
    print('~~ Processing sequences ~~')
    start_time = time.time()
    total_cpg, total_chg, total_chh = 0, 0, 0

    mononucleotide_totals = np.zeros(4)  # A, C, G, T
    dinucleotide_totals = np.zeros(16)  # AA, AC, ..., TT

    with Pool(CPU_COUNT) as pool:
        results = pool.map(count_nucleotides_and_motifs, file_sequences)

    for result in results:
        if isinstance(result, int):
            continue  # 跳過無效結果

        cpg, chg, chh, mononucleotide_counts, dinucleotide_counts = result

        total_cpg += cpg
        total_chg += chg
        total_chh += chh

        mononucleotide_totals += np.array([mononucleotide_counts['A'], mononucleotide_counts['C'],
                                           mononucleotide_counts['G'], mononucleotide_counts['T']])
        dinucleotide_totals += np.array([dinucleotide_counts[dinuc] for dinuc in sorted(dinucleotide_counts.keys())])

    total_motifs = total_cpg + total_chg + total_chh
    cpg_chg_chh_probs = [f'{total_cpg / total_motifs:.4f}', f'{total_chg / total_motifs:.4f}',
                         f'{total_chh / total_motifs:.4f}']

    total_mono = mononucleotide_totals.sum()
    mononucleotide_probs = (mononucleotide_totals / total_mono).round(4)
    for prob in mononucleotide_probs:
        cpg_chg_chh_probs.append(f'{prob:.4f}')

    dinucleotide_list = dinucleotide_totals.tolist()
    m_factors = []
    index = 0
    for i in range(0, len(dinucleotide_list), 4):
        dinuc_group = dinucleotide_list[i:i + 4]
        group_sum = sum(dinuc_group)
        if group_sum == 0:
            m_factor = 0
        else:
            m_factor = mononucleotide_probs[index] / group_sum
        m_factors.extend([m_factor] * 4)
        index += 1

    m_factors = np.array(m_factors)
    dinucleotide_probs = (dinucleotide_totals * m_factors).round(4)
    for prob in dinucleotide_probs:
        cpg_chg_chh_probs.append(f'{prob:.4f}')

    elapsed_time = time.time() - start_time
    print(f'Calculation finished in {elapsed_time / 60:.2f} minutes.')

    index_labels = ['CpG', 'CHG', 'CHH', 'A', 'C', 'G', 'T'] + sorted(dinucleotide_counts.keys())
    df = pd.DataFrame(cpg_chg_chh_probs, columns=[directory], index=index_labels)

    output_path = os.path.join(dir_path, '../../Background_probability', directory,
                               f'{species}_{TF}_{celltype}_{region}_probability.txt')
    df.to_csv(output_path, sep='\t')


def read_big_file(big_fileA, big_fileB):
    """
    生成器函數，逐行讀取兩個大文件。

    參數:
    - big_fileA (file object): 第一個文件對象。
    - big_fileB (file object): 第二個文件對象。

    產出:
    - 每次產出兩個文件對應的行，作為元組。
    """
    while True:
        data_chunk_from_fileA = big_fileA.readline()
        data_chunk_from_fileB = big_fileB.readline()
        if data_chunk_from_fileA == '' or data_chunk_from_fileB == '':
            break

        yield data_chunk_from_fileA.strip(), data_chunk_from_fileB.strip()


def calc_mlevel(fileA, fileB):
    """
    計算甲基化水平。

    參數:
    - fileA (str): 第一個文件的路徑。
    - fileB (str): 第二個文件的路徑。

    返回:
    - 甲基化概率（str），格式化為小數點後四位。
    """
    print('~~ Calculating methylation level ~~')
    total, count = 0, 0
    length = 900000000  # 讀取的最大行數
    s_time = time.time()
    with Pool(CPU_COUNT - 2) as pool:
        with open(fileA, 'r') as big_fileA, open(fileB, 'r') as big_fileB:
            for wgbs_from_fileA, wgbs_from_fileB in itertools.islice(read_big_file(big_fileA, big_fileB), length):
                wgbs_from_fileA = wgbs_from_fileA.split('\t')
                wgbs_from_fileB = wgbs_from_fileB.split('\t')
                try:
                    reads_A = int(wgbs_from_fileA[-2])
                    reads_B = int(wgbs_from_fileB[-2])
                    meth_A = int(wgbs_from_fileA[-1]) / 100
                    meth_B = int(wgbs_from_fileB[-1]) / 100

                    if reads_A < 4 and reads_B < 4:
                        continue
                    elif reads_A == 0:
                        total += meth_B
                        count += 1
                    elif reads_B == 0:
                        total += meth_A
                        count += 1
                    else:
                        read_A = min(reads_A, 4)
                        read_B = min(reads_B, 4)
                        weighted_meth = (meth_A * read_A + meth_B * read_B) / (read_A + read_B)
                        total += weighted_meth
                        count += 1
                except:
                    pass

            pool.close()
            pool.join()

    end = time.time()
    print('耗費時間：', (end - s_time) / 60, 'min')
    print('甲基化概率：', total / count, '\n')
    return f'{total / count:.4f}'

def calculate_methylation_level_whole_genome(wgbs_files, species_name, celltype, region, directory, TF):
    print('~~ Calculating methylation level for whole genome ~~')
    names = ['CG', 'CHG', 'CHH']
    total = []
    count = 0
    for files in wgbs_files:
        print(names[count])
        start = time.time()

        wgbs_file_A = files[0]
        wgbs_file_B = files[1]

        # 计算甲基化水平
        result = calc_mlevel_whole_genome(wgbs_file_A, wgbs_file_B)
        total.append(result)
        count += 1

        end = time.time()
        print(f"Processing {names[count-1]} finished, cost {(end - start) / 60:.2f} min\n")

    df = pd.DataFrame(total, columns=[directory], index=['mCG', 'mCHG', 'mCHH'])
    print('~~ Saving results ~~')
    output_path = os.path.join(dir_path, '../../Background_probability', directory,
                               f'{species_name}_{celltype}_{region}_methyl_probability.txt')
    df.to_csv(output_path, sep='\t', float_format='%.4f')
    print('~~ Results saved ~~')


def calc_mlevel_stream(wgbs_tfbs_A, wgbs_tfbs_B):
    print('~~ Calculating methylation level ~~')
    total, count = 0, 0
    s_time = time.time()
    
    # 使用 zip 迭代两个生成器
    for line_A, line_B in zip(wgbs_tfbs_A, wgbs_tfbs_B):
        wgbs_from_fileA = str(line_A).strip().split('\t')
        wgbs_from_fileB = str(line_B).strip().split('\t')
        try:
            reads_A = int(wgbs_from_fileA[-2])
            reads_B = int(wgbs_from_fileB[-2])
            meth_A = float(wgbs_from_fileA[-1]) / 100
            meth_B = float(wgbs_from_fileB[-1]) / 100

            if reads_A < 4 and reads_B < 4:
                continue
            elif reads_A == 0:
                total += meth_B
                count += 1
            elif reads_B == 0:
                total += meth_A
                count += 1
            else:
                read_A = min(reads_A, 4)
                read_B = min(reads_B, 4)
                weighted_meth = (meth_A * read_A + meth_B * read_B) / (read_A + read_B)
                total += weighted_meth
                count += 1
        except:
            pass

    end = time.time()
    print('耗费时间：', (end - s_time) / 60, 'min')
    if count == 0:
        print('没有有效的数据点。')
        return 0.0
    print('甲基化概率：', total / count, '\n')
    return f'{total / count:.4f}'


import os
import pybedtools

def create_genome_file(fasta_path, genome_file_path):
    fai_path = fasta_path + '.fai'
    if not os.path.exists(fai_path):
        raise FileNotFoundError(f"FASTA 索引文件 {fai_path} 不存在。请先使用 'samtools faidx' 创建索引文件。")
    with open(fai_path, 'r') as fai_file, open(genome_file_path, 'w') as genome_file:
        for line in fai_file:
            fields = line.strip().split('\t')
            chrom = fields[0]
            size = fields[1]
            genome_file.write(f"{chrom}\t{size}\n")

def define_promoter_regions(tss_bed, upstream=1000, downstream=200, genome_file_path=None):
    if genome_file_path is None:
        raise ValueError("必须提供 genome_file_path 参数，以避免区域超出染色体边界。")
    # 使用 'g' 参数指定基因组文件路径
    promoter_bed = tss_bed.slop(
        l=upstream,
        r=downstream,
        s=True,
        g=genome_file_path
    )
    return promoter_bed

def promoter(species, wgbs_files, celltype, TF=None):
    tss_bed_path = "/home/wayne/MethylSeqLogo_automation/bin/TSS/TSS.bed"
    tss_bed = pybedtools.BedTool(tss_bed_path)
    # 基因组文件路径
    genome_file_path = os.path.join(dir_path, '../genome', f'{species}.genome')
    # 创建基因组文件（如果尚未创建）
    fasta_path = os.path.join(dir_path, '../genome', f'{species}.fa')
    if not os.path.exists(genome_file_path):
        create_genome_file(fasta_path, genome_file_path)
    # 定义启动子区域
    promoter_bed = define_promoter_regions(tss_bed, upstream=1000, downstream=200, genome_file_path=genome_file_path)
    # 后续代码
    # 读取基因组序列
    fasta = pybedtools.BedTool(fasta_path)
    # 提取启动子序列
    tfbss_fa = promoter_bed.sequence(fi=fasta, s=True)
    with open(tfbss_fa.seqfn, 'r') as f:
        tfbss_fasta = f.readlines()
    # 处理序列以计算核苷酸和基序计数
    process_sequences(tfbss_fasta, species, celltype, 'promoter', 'promoter', TF)
    # 计算甲基化水平
    calculate_methylation_level(promoter_bed, wgbs_files, species, celltype, 'promoter', 'promoter', TF)
    # 读取背景概率表
    return read_bgprob_table(species, celltype, 'promoter', TF)


def neighbor(tfbs_bed, species, wgbs_files, celltype, region, TF):
    """
    處理鄰近區域的背景概率和甲基化水平。

    參數:
    - tfbs_bed (BedTool): TFBS的BedTool對象。
    - species (str): 物種名稱。
    - wgbs_files (list): WGBS文件列表。
    - celltype (str): 細胞類型。
    - region (str): 區域。
    - TF (str): 轉錄因子名稱。

    返回:
    - 調用read_bgprob_table的結果。
    """
    fasta = pybedtools.BedTool(os.path.join(dir_path, '../genome', f'{species}.fa'))
    tfbss_fa = tfbs_bed.sequence(fi=fasta, s=True)
    with open(tfbss_fa.seqfn, 'r') as f:
        tfbss_fasta = f.readlines()
    process_sequences(tfbss_fasta, species, celltype, region, 'neighbor', TF)
    calculate_methylation_level(tfbs_bed, wgbs_files, species, celltype, region, 'neighbor', TF)

    return read_bgprob_table(species, celltype, region, TF)

def whole(species, wgbs_files, celltype, TF=None):
    """
    处理全基因组，计算核苷酸背景概率和甲基化水平。

    参数:
    - species (str): 物种名称。
    - wgbs_files (list): WGBS文件列表。
    - celltype (str): 细胞类型。
    - TF (str): 转录因子名称（可选）。

    返回:
    - 调用 read_bgprob_table 的结果。
    """
    # 读取基因组序列
    fasta_path = os.path.join(dir_path, '../genome', f'{species}.fa')
    fasta = pybedtools.BedTool(fasta_path)
    
    # 处理全基因组序列
    process_sequences_large(fasta_path, species, celltype, 'whole_genome', 'whole_genome', TF)
    
    # 计算全基因组甲基化水平
    calculate_methylation_level_whole_genome(wgbs_files, species, celltype, 'whole_genome', 'whole_genome', TF)
    
    # 返回背景概率表和甲基化水平
    return read_bgprob_table(species, celltype, 'whole_genome', TF)


def calc_mlevel_whole_genome(wgbs_file_A, wgbs_file_B):
    print('~~ Calculating methylation level ~~')
    total, count = 0, 0
    s_time = time.time()

    # 使用生成器逐行读取文件，避免将整个文件加载到内存中
    def read_wgbs_file(file_path):
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                chrom, pos = parts[0], int(parts[1])
                reads = int(parts[-2])
                meth = float(parts[-1]) / 100
                yield (chrom, pos, reads, meth)

    # 创建迭代器
    iter_A = read_wgbs_file(wgbs_file_A)
    iter_B = read_wgbs_file(wgbs_file_B)

    # 使用合并排序的方式同时遍历两个文件
    try:
        item_A = next(iter_A)
        item_B = next(iter_B)
        while True:
            if item_A[:2] == item_B[:2]:
                # 位点匹配，进行计算
                reads_A, meth_A = item_A[2], item_A[3]
                reads_B, meth_B = item_B[2], item_B[3]
                if reads_A < 4 and reads_B < 4:
                    pass  # 跳过
                elif reads_A == 0:
                    total += meth_B
                    count += 1
                elif reads_B == 0:
                    total += meth_A
                    count += 1
                else:
                    read_A = min(reads_A, 4)
                    read_B = min(reads_B, 4)
                    weighted_meth = (meth_A * read_A + meth_B * read_B) / (read_A + read_B)
                    total += weighted_meth
                    count += 1
                item_A = next(iter_A)
                item_B = next(iter_B)
            elif item_A[:2] < item_B[:2]:
                item_A = next(iter_A)
            else:
                item_B = next(iter_B)
    except StopIteration:
        pass  # 遍历完所有数据

    end = time.time()
    print('耗费时间：', (end - s_time) / 60, 'min')
    if count == 0:
        print('没有有效的数据点。')
        return '0.0000'
    print('甲基化概率：', total / count, '\n')
    return f'{total / count:.4f}'

def calculate_methylation_level_whole_genome(wgbs_files, species_name, celltype, region, directory, TF):
    print('~~ Calculating methylation level for whole genome ~~')
    names = ['CG', 'CHG', 'CHH']
    total = []
    count = 0
    for files in wgbs_files:
        print(names[count])
        start = time.time()

        wgbs_file_A = files[0]
        wgbs_file_B = files[1]

        # 计算甲基化水平
        result = calc_mlevel_whole_genome(wgbs_file_A, wgbs_file_B)
        total.append(result)
        count += 1

        end = time.time()
        print(f"Processing {names[count-1]} finished, cost {(end - start) / 60:.2f} min\n")

    df = pd.DataFrame(total, columns=[directory], index=['mCG', 'mCHG', 'mCHH'])
    print('~~ Saving results ~~')
    output_path = os.path.join(dir_path, '../../Background_probability', directory,
                               f'{species_name}_{celltype}_{region}_methyl_probability.txt')
    df.to_csv(output_path, sep='\t', float_format='%.4f')
    print('~~ Results saved ~~')

def process_sequences_large(file_path, species, celltype, region, directory, TF):
    """
    处理大型基因组序列文件，计算核苷酸和基序计数。

    参数:
    - file_path (str): 基因组 FASTA 文件的路径。
    - species (str): 物种名称。
    - celltype (str): 细胞类型。
    - region (str): 区域。
    - directory (str): 保存输出文件的目录。
    - TF (str): 转录因子名称。
    """
    print('~~ Processing sequences ~~')
    start_time = time.time()
    total_cpg, total_chg, total_chh = 0, 0, 0

    mononucleotide_totals = np.zeros(4)  # A, C, G, T
    dinucleotide_totals = np.zeros(16)  # AA, AC, ..., TT

    # 打开基因组 FASTA 文件，逐步读取序列
    with open(file_path, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                # 处理前一个序列
                if seq:
                    result = count_nucleotides_and_motifs(seq)
                    if isinstance(result, int):
                        continue  # 跳过无效结果

                    cpg, chg, chh, mononucleotide_counts, dinucleotide_counts = result

                    total_cpg += cpg
                    total_chg += chg
                    total_chh += chh

                    mononucleotide_totals += np.array([
                        mononucleotide_counts['A'],
                        mononucleotide_counts['C'],
                        mononucleotide_counts['G'],
                        mononucleotide_counts['T']
                    ])
                    dinucleotide_totals += np.array([
                        dinucleotide_counts[dinuc] for dinuc in sorted(dinucleotide_counts.keys())
                    ])

                seq = ''  # 重置序列
            else:
                seq += line.strip().upper()
        # 处理最后一个序列
        if seq:
            result = count_nucleotides_and_motifs(seq)
            if not isinstance(result, int):
                cpg, chg, chh, mononucleotide_counts, dinucleotide_counts = result

                total_cpg += cpg
                total_chg += chg
                total_chh += chh

                mononucleotide_totals += np.array([
                    mononucleotide_counts['A'],
                    mononucleotide_counts['C'],
                    mononucleotide_counts['G'],
                    mononucleotide_counts['T']
                ])
                dinucleotide_totals += np.array([
                    dinucleotide_counts[dinuc] for dinuc in sorted(dinucleotide_counts.keys())
                ])

    total_motifs = total_cpg + total_chg + total_chh
    cpg_chg_chh_probs = [
        f'{total_cpg / total_motifs:.4f}',
        f'{total_chg / total_motifs:.4f}',
        f'{total_chh / total_motifs:.4f}'
    ]

    total_mono = mononucleotide_totals.sum()
    mononucleotide_probs = (mononucleotide_totals / total_mono).round(4)
    cpg_chg_chh_probs.extend([f'{prob:.4f}' for prob in mononucleotide_probs])

    dinucleotide_list = dinucleotide_totals.tolist()
    m_factors = []
    index = 0
    for i in range(0, len(dinucleotide_list), 4):
        dinuc_group = dinucleotide_list[i:i + 4]
        group_sum = sum(dinuc_group)
        if group_sum == 0:
            m_factor = 0
        else:
            m_factor = mononucleotide_probs[index] / group_sum
        m_factors.extend([m_factor] * 4)
        index += 1

    m_factors = np.array(m_factors)
    dinucleotide_probs = (dinucleotide_totals * m_factors).round(4)
    cpg_chg_chh_probs.extend([f'{prob:.4f}' for prob in dinucleotide_probs])

    elapsed_time = time.time() - start_time
    print(f'Calculation finished in {elapsed_time / 60:.2f} minutes.')

    index_labels = ['CpG', 'CHG', 'CHH', 'A', 'C', 'G', 'T'] + sorted(dinucleotide_counts.keys())
    df = pd.DataFrame(cpg_chg_chh_probs, columns=[directory], index=index_labels)

    output_path = os.path.join(dir_path, '../../Background_probability', directory,
                               f'{species}_{region}_probability.txt')
    df.to_csv(output_path, sep='\t')
