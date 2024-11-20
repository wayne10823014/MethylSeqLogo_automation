import os
import time
from collections import defaultdict

import pandas as pd
import pyranges as pr  # 使用 PyRanges 替代 pybedtools


def init_list():
    """初始化一个空列表，用于 defaultdict。"""
    return []


def flanking_bed(tfbs_of_tf_in_celltype, span1, span2):
    """
    生成 TFBS 两侧的序列，用于邻近区域分析。

    参数：
    - tfbs_of_tf_in_celltype (PyRanges): TFBS 在特定细胞类型中的位置信息。
    - span1 (int): 左侧延伸的碱基数。
    - span2 (int): 右侧延伸的碱基数。

    返回：
    - modify_coordinate (PyRanges): 修改后的位置信息。
    """
    # 将每个区域左右延伸
    df = tfbs_of_tf_in_celltype.df.copy()
    left_flank = df.copy()
    left_flank['End'] = left_flank['Start']
    left_flank['Start'] = (left_flank['Start'] - span1).clip(lower=0)

    right_flank = df.copy()
    right_flank['Start'] = right_flank['End']
    right_flank['End'] = right_flank['End'] + span2

    modify_coordinate_df = pd.concat([left_flank, right_flank], ignore_index=True)
    modify_coordinate = pr.PyRanges(modify_coordinate_df)

    print('Modified coordinates:', len(modify_coordinate_df), len(modify_coordinate))
    return modify_coordinate


def modify_bed(tfbs_of_tf_in_celltype, span1, span2):
    """
    修改 TFBS 的坐标，延伸左右两侧，用于获取上下文序列。

    参数：
    - tfbs_of_tf_in_celltype (PyRanges): TFBS 在特定细胞类型中的位置信息。
    - span1 (int): 左侧延伸的碱基数。
    - span2 (int): 右侧延伸的碱基数。

    返回：
    - modify_coordinate_list (list): 修改后的坐标列表。
    - modify_coordinate (PyRanges): 修改后的位置信息。
    """
    df = tfbs_of_tf_in_celltype.df.copy()
    df['Start'] = (df['Start'] - span1).clip(lower=0)
    df['End'] = df['End'] + span2

    modify_coordinate_list = df.values.tolist()
    modify_coordinate = pr.PyRanges(df)

    print('Modified coordinates:', len(modify_coordinate_list), len(modify_coordinate))
    return modify_coordinate_list, modify_coordinate


def get_tfbs(chip_seq_of_celltype, tfbs_file, species, span1, span2):
    """
    获取特定细胞类型中 TF 的 TFBS（motif）。

    参数：
    - chip_seq_of_celltype (str): ChIP-seq 数据文件路径。
    - tfbs_file (str): TFBS（motif）文件路径。
    - species (str): 物种名称。
    - span1 (int): 左侧延伸的碱基数。
    - span2 (int): 右侧延伸的碱基数。

    返回：
    - tfbs_of_tf_in_celltype (PyRanges): 交集得到的 TFBS。
    - modify_coordinate_list (list): 修改后的坐标列表。
    - motif_len (int): motif 的长度。
    """
    import pyranges as pr
    from pyfaidx import Fasta
    import os

    tfbs_site = pr.read_bed(tfbs_file)
    celltype = pr.read_bed(chip_seq_of_celltype)
    tfbs_of_tf_in_celltype = tfbs_site.join(celltype)

    # 获取序列
    fasta_path = os.path.join(dir_path, '../genome', f'{species}.fa')
    fasta = Fasta(fasta_path)

    sequences = []
    for idx, row in tfbs_of_tf_in_celltype.df.iterrows():
        chrom = row['Chromosome']
        start = row['Start']
        end = row['End']
        seq = fasta[chrom][start:end]
        sequences.append(str(seq))

    motif_len = len(sequences[0]) if sequences else 0

    print('\nIntersection finished...')
    print(f'{len(tfbs_of_tf_in_celltype)} transcription factor binding sites discovered')

    modify_coordinate_list, modify_coordinate = modify_bed(tfbs_of_tf_in_celltype, span1, span2)

    return tfbs_of_tf_in_celltype, modify_coordinate_list, modify_coordinate, motif_len


def get_seq(species, seq_list, tfbs_modify_coordinate, span1, span2, motif_len):
    """
    获取序列，并构建 total 数据结构。

    参数：
    - species (str): 物种名称。
    - seq_list (list): 修改后的坐标列表。
    - tfbs_modify_coordinate (PyRanges): 修改后的位置信息。
    - span1 (int): 左侧延伸的碱基数。
    - span2 (int): 右侧延伸的碱基数。
    - motif_len (int): motif 的长度。

    返回：
    - total (defaultdict): 存储序列信息的字典。
    - seq_dict (dict): 序列字典，键为序列标识符，值为序列字符串。
    """
    total = defaultdict(init_list)
    from pyfaidx import Fasta
    fasta_path = os.path.join(dir_path, '../genome', f'{species}.fa')
    fasta = Fasta(fasta_path)

    seq_dict = {}
    count = 0
    for idx, row in tfbs_modify_coordinate.df.iterrows():
        chrom = row['Chromosome']
        start = row['Start']
        end = row['End']
        strand = row.get('Strand', '+')
        sequence = str(fasta[chrom][start:end])

        tfbsstartpos = int(seq_list[count][1]) + span1
        tfbsendpos = int(seq_list[count][2]) - span2
        key = f'>{chrom}:{tfbsstartpos}-{tfbsendpos}{strand}'

        seq_dict[key] = sequence
        total[key].append([0] * motif_len)  # ctx
        total[key].append([0] * motif_len)  # tread
        total[key].append([0] * motif_len)  # cread

        for base in range(span1, span1 + motif_len):
            if sequence[base] not in ["c", "C", "g", "G"]:
                continue
            else:
                temp = total[key][0]
                # 以下根据链的方向和碱基类型标记上下文
                if strand == '+' and sequence[base] in ["c", "C"]:
                    if sequence[base + 1] in ["g", "G"]:
                        temp[base - span1] = 'X'
                    elif sequence[base + 2] in ["g", "G"]:
                        temp[base - span1] = 'Y'
                    else:
                        temp[base - span1] = 'Z'
                elif strand == '+' and sequence[base] in ["G", "g"]:
                    if sequence[base - 1] in ["c", "C"]:
                        temp[base - span1] = 'x'
                    elif sequence[base - 2] in ["c", "C"]:
                        temp[base - span1] = 'y'
                    else:
                        temp[base - span1] = 'z'
                elif strand == '-' and sequence[base] in ["G", "g"]:
                    if sequence[base - 1] in ["c", "C"]:
                        temp[base - span1] = 'x'
                    elif sequence[base - 2] in ["c", "C"]:
                        temp[base - span1] = 'y'
                    else:
                        temp[base - span1] = 'z'
                elif strand == '-' and sequence[base] in ["c", "C"]:
                    if sequence[base + 1] in ["g", "G"]:
                        temp[base - span1] = 'X'
                    elif sequence[base + 2] in ["g", "G"]:
                        temp[base - span1] = 'Y'
                    else:
                        temp[base - span1] = 'Z'
        count += 1

    return total, seq_dict



def methylread_counter(TFBSFile, WGBSFile, total):
    print('\nStarting methylread_counter...')
    TFBSs = TFBSFile
    WGBSFile = [filepath for sublist in WGBSFile for filepath in sublist]

    for count, file in enumerate(WGBSFile, start=1):
        print(f"Processing file {count}: {file}")
        try:
            # 将 TFBS 转换为 DataFrame，以便后续处理
            tfbs_df = TFBSs.df

            # 使用迭代器逐块读取 WGBS 文件
            chunksize = 10 ** 6  # 每次读取 1,000,000 行，可根据内存情况调整
            col_names = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand',
                         'ThickStart', 'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes']

            for chunk in pd.read_csv(file, sep='\t', names=col_names, chunksize=chunksize, header=None):
                # 将 chunk 转换为 PyRanges 对象
                WGBS_chunk = pr.PyRanges(chunk)

                # 进行交集操作
                Methylation_of_TFBSs = WGBS_chunk.join(TFBSs, suffix="_TFBS")

                df = Methylation_of_TFBSs.df

                if df.empty:
                    continue

                for idx, row in df.iterrows():
                    try:
                        tfbsstartpos = row['Start_TFBS']
                        tfbsendpos = row['End_TFBS']
                        chrom = row['Chromosome']
                        base_MethylationInfo_pos = row['Start']
                        strand = row.get('Strand_TFBS', '+')

                        methylated_reads = int(row['BlockCount'])
                        unmethylated_reads = int(row['BlockSizes'])

                        capped_read = methylated_reads + unmethylated_reads
                        methyl_read = methylated_reads

                        key = f'>{chrom}:{tfbsstartpos}-{tfbsendpos}{strand}'
                        if key not in total:
                            continue

                        if strand == '+':
                            temp = total[key][1]
                            idx_pos = int(base_MethylationInfo_pos) - int(tfbsstartpos)
                            if idx_pos < 0 or idx_pos >= len(temp):
                                continue
                            temp[idx_pos] += unmethylated_reads
                            temp1 = total[key][2]
                            temp1[idx_pos] += methyl_read
                        else:
                            temp = total[key][1]
                            idx_pos = int(tfbsendpos) - int(base_MethylationInfo_pos) - 1
                            if idx_pos < 0 or idx_pos >= len(temp):
                                continue
                            temp[idx_pos] += unmethylated_reads
                            temp1 = total[key][2]
                            temp1[idx_pos] += methyl_read
                    except Exception as e:
                        print(f"Error processing row {idx} in file {file}: {e}")
                        continue

        except Exception as e:
            print(f"Error occurred while processing file {count}: {file}")
            print(f"Exception: {e}")
            continue

    # 构建并返回结果 DataFrame
    ctx = []
    tread = []
    cread = []
    for k in total:
        try:
            ctx.append(total[k][0])
            tread.append(total[k][1])
            cread.append(total[k][2])
        except IndexError:
            pass

    print("Saving ctx_dict to CSV...")
    ctx_dict = pd.DataFrame(ctx)
    cread_dict = pd.DataFrame(cread)
    tread_dict = pd.DataFrame(tread)

    return ctx_dict, cread_dict, tread_dict




# 获取当前脚本的目录路径
dir_path = os.path.dirname(os.path.realpath(__file__))
