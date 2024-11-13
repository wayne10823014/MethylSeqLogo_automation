import argparse
import time
from multiprocessing import Pool, cpu_count
from functools import reduce
import numpy as np
import os
import itertools
import pandas as pd

def get_parser():
    """
    创建命令行参数解析器。
    """
    parser = argparse.ArgumentParser(description='背景概率计算程序')
    parser.add_argument(
        '-c', '--celltype',
        choices=[
            'H1-hESC', 'HepG2', 'K562', 'GM12878', 'HeLa-S3', 'naive-hESC', 'primed-hESC',
            'mESC', 'mMPGC_E13p5', 'mMPGC_E16p5', 'mFPGC_E13p5', 'mFPGC_E16p5',
            'Leaf', 'Inflorescence', 'Root', 'shoot'
        ],
        default='HepG2',
        metavar='cell types',
        help='WGBS 文件所属的细胞类型'
    )
    parser.add_argument(
        '-s', '--speciesname',
        default='human',
        metavar='species',
        help='输入数据的物种名称，用于 MethylSeqLogo'
    )
    parser.add_argument(
        '-fa', '--fasta',
        default='/genome/human.fa',
        metavar='species',
        help='物种的 FASTA 文件路径'
    )

    wgbs_group = parser.add_argument_group('甲基化计算参数')
    wgbs_group.add_argument('-CG', '--mcg', nargs=2, type=str, help='CpG 位点的 WGBS 文件')
    wgbs_group.add_argument('-CHG', '--mchg', nargs=2, type=str, help='CHG 位点的 WGBS 文件')
    wgbs_group.add_argument('-CHH', '--mchh', nargs=2, type=str, help='CHH 位点的 WGBS 文件')

    return parser

def read_big_file(file_a, file_b):
    """
    按行读取两个大文件，生成器形式返回每一行的数据。

    参数：
    - file_a, file_b: 文件对象
    """
    while True:
        line_a = file_a.readline()
        line_b = file_b.readline()
        if not line_a or not line_b:
            break
        yield line_a.strip(), line_b.strip()

def calc_mlevel(file_a_path, file_b_path):
    """
    计算甲基化水平。

    参数：
    - file_a_path, file_b_path: 两个文件的路径

    返回：
    - total_mlevel: 总的甲基化水平，保留四位小数的字符串
    """
    total, count = 0, 0
    max_length = 900_000_000  # 最大读取长度

    with open(file_a_path, 'r') as file_a, open(file_b_path, 'r') as file_b:
        for line_a, line_b in itertools.islice(read_big_file(file_a, file_b), max_length):
            wgbs_a = line_a.split('\t')
            wgbs_b = line_b.split('\t')
            try:
                reads_a = int(wgbs_a[-2])
                reads_b = int(wgbs_b[-2])
                methyl_a = int(wgbs_a[-1])
                methyl_b = int(wgbs_b[-1])

                if reads_a < 4 and reads_b < 4:
                    continue
                elif reads_a == 0:
                    total += methyl_b / 100
                elif reads_b == 0:
                    total += methyl_a / 100
                else:
                    reads_a = min(reads_a, 4)
                    reads_b = min(reads_b, 4)
                    weighted_methyl = ((methyl_a * reads_a + methyl_b * reads_b) / 100) / (reads_a + reads_b)
                    total += weighted_methyl
                count += 1
            except (ValueError, IndexError):
                continue

    if count == 0:
        return '0.0000'
    return f'{total / count:.4f}'

def calc(seq):
    """
    计算给定序列的单核苷酸、二核苷酸频率以及 CpG、CHG、CHH 的计数。

    参数：
    - seq: DNA 序列字符串

    返回：
    - [cpg_count, chg_count, chh_count, mononucleotide_dict, dinucleotide_dict]
    """
    if seq.startswith('>') or (seq.startswith('N') and len(set(seq)) == 2):
        return 0

    dinucleotide_dict = {k: 0 for k in ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
                                        'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']}
    mononucleotide_dict = {k: 0 for k in ['A', 'C', 'G', 'T']}
    cpg_count, chg_count, chh_count = 0, 0, 0

    seq = seq.upper()
    seq_len = len(seq)

    for i in range(seq_len):
        base = seq[i]
        if base in mononucleotide_dict:
            mononucleotide_dict[base] += 1

            if i < seq_len - 1 and seq[i + 1] in mononucleotide_dict:
                dinucleotide = base + seq[i + 1]
                if dinucleotide in dinucleotide_dict:
                    dinucleotide_dict[dinucleotide] += 1

        if 2 <= i <= seq_len - 3:
            if base == 'C':
                if seq[i + 1] == 'G':
                    cpg_count += 1
                elif seq[i + 2] == 'G':
                    chg_count += 1
                else:
                    chh_count += 1
            elif base == 'G':
                if seq[i - 1] == 'C':
                    cpg_count += 1
                elif seq[i - 2] == 'C':
                    chg_count += 1
                else:
                    chh_count += 1

    return [cpg_count, chg_count, chh_count, mononucleotide_dict, dinucleotide_dict]

def read_file(sequences):
    """
    读取序列文件，计算核苷酸和甲基化相关信息，并将结果保存为文件。

    参数：
    - sequences: 序列列表
    """
    start_time = time.time()
    cpg_total, chg_total, chh_total = 0, 0, 0
    mononucleotide_counts = np.zeros(4)
    dinucleotide_counts = np.zeros(16)

    with Pool(cpu_count()) as pool:
        results = pool.map(calc, sequences)

    for result in results:
        if isinstance(result, int):
            continue
        cpg_total += result[0]
        chg_total += result[1]
        chh_total += result[2]
        mononucleotide_counts += np.array(list(result[3].values()))
        dinucleotide_counts += np.array(list(result[4].values()))

    total_cpg_chg_chh = cpg_total + chg_total + chh_total
    if total_cpg_chg_chh == 0:
        total_cpg_chg_chh = 1  # 防止除以零

    cpg_chg_chh_probs = [
        f'{cpg_total / total_cpg_chg_chh:.4f}',
        f'{chg_total / total_cpg_chg_chh:.4f}',
        f'{chh_total / total_cpg_chg_chh:.4f}'
    ]

    mono_total = mononucleotide_counts.sum()
    if mono_total == 0:
        mono_total = 1  # 防止除以零

    mononucleotide_probs = mononucleotide_counts / mono_total
    mononucleotide_probs = [f'{prob:.4f}' for prob in mononucleotide_probs]

    cpg_chg_chh_probs.extend(mononucleotide_probs)

    dinucleotide_probs = dinucleotide_counts.copy()
    index = 0
    for i in range(0, len(dinucleotide_counts), 4):
        dinucleotide_group = dinucleotide_counts[i:i+4]
        group_sum = dinucleotide_group.sum()
        if group_sum == 0:
            group_sum = 1  # 防止除以零
        scale_factor = mononucleotide_counts[index] / group_sum
        dinucleotide_probs[i:i+4] *= scale_factor
        index += 1

    dinucleotide_probs = [f'{prob:.4f}' for prob in dinucleotide_probs]
    cpg_chg_chh_probs.extend(dinucleotide_probs)

    end_time = time.time()
    print(f'计算完成，耗时 {(end_time - start_time) / 60:.2f} 分钟')

    nucleotides = ['CpG', 'CHG', 'CHH', 'A', 'C', 'G', 'T']
    dinucleotides = [
        'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
        'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'
    ]
    index_labels = nucleotides + dinucleotides

    df = pd.DataFrame(cpg_chg_chh_probs, columns=['whole_genome'], index=index_labels)
    output_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        '..',
        'Background_probability',
        'whole_genome',
        f'{speciesname}_whole_genome_probability.txt'
    )
    df.to_csv(output_path, sep='\t')
    print(f'结果已保存到 {output_path}')

def main():
    global cpu_count, speciesname, wgbs_file, celltype

    parser = get_parser()
    args = parser.parse_args()

    speciesname = args.speciesname
    celltype = args.celltype
    wgbs_file = [args.mcg, args.mchg, args.mchh]
    fasta_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), args.fasta)

    print(f'\n文件路径：{fasta_path}')
    print('计算大约需要 10 分钟时间。\n')

    start_time = time.time()
    with open(fasta_path, 'r') as fasta_file:
        sequences = fasta_file.readlines()
    end_time = time.time()
    print(f'读取文件完成，耗时 {(end_time - start_time) / 60:.2f} 分钟\n')

    read_file(sequences)

    if wgbs_file[0]:
        print('开始计算甲基化水平...\n')

        for label, files in zip(['CHG', 'CHH', 'CG'], wgbs_file[1:] + [wgbs_file[0]]):
            if files and len(files) == 2:
                start = time.time()
                result = calc_mlevel(files[0], files[1])
                end = time.time()
                print(f'{label} 文件计算完成，耗时 {(end - start) / 60:.2f} 分钟\n')
                wgbs_result[label] = result
            else:
                wgbs_result[label] = '0.0000'

        df = pd.DataFrame.from_dict(wgbs_result, orient='index', columns=['whole_genome'])
        output_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            '..',
            'Background_probability',
            'whole_genome',
            f'{speciesname}_{celltype}_whole_genome_methyl_probability.txt'
        )
        df.to_csv(output_path, sep='\t', float_format='%.4f')
        print(f'甲基化水平结果已保存到 {output_path}')

if __name__ == '__main__':
    wgbs_result = {}
    main()
