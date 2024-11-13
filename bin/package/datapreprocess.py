import os
import time
from collections import defaultdict

import pandas as pd
import pybedtools


def init_list():
    """初始化一个空列表，用于 defaultdict。"""
    return []


def flanking_bed(tfbs_of_tf_in_celltype, span1, span2):
    """
    生成 TFBS 两侧的序列，用于邻近区域分析。

    参数：
    - tfbs_of_tf_in_celltype (BedTool): TFBS 在特定细胞类型中的位置信息。
    - span1 (int): 左侧延伸的碱基数。
    - span2 (int): 右侧延伸的碱基数。

    返回：
    - modify_coordinate (BedTool): 修改后的位置信息。
    """
    modify_coordinate_list = []
    for tfbs in tfbs_of_tf_in_celltype:
        tfbs = list(tfbs)
        # 左侧延伸
        start = str(max(0, int(tfbs[1]) - span1))
        modify_coordinate_list.append([tfbs[0], start, tfbs[1]])
        # 右侧延伸
        end = str(int(tfbs[2]) + span2)
        modify_coordinate_list.append([tfbs[0], tfbs[2], end])

    gene_df = pd.DataFrame(modify_coordinate_list)
    modify_coordinate = pybedtools.BedTool.from_dataframe(gene_df)
    print('Modified coordinates:', len(modify_coordinate_list), len(modify_coordinate))
    return modify_coordinate


def modify_bed(tfbs_of_tf_in_celltype, span1, span2):
    """
    修改 TFBS 的坐标，延伸左右两侧，用于获取上下文序列。

    参数：
    - tfbs_of_tf_in_celltype (BedTool): TFBS 在特定细胞类型中的位置信息。
    - span1 (int): 左侧延伸的碱基数。
    - span2 (int): 右侧延伸的碱基数。

    返回：
    - modify_coordinate_list (list): 修改后的坐标列表。
    - modify_coordinate (BedTool): 修改后的位置信息。
    """
    modify_coordinate_list = []
    for tfbs in tfbs_of_tf_in_celltype:
        tfbs = list(tfbs)
        tfbs[1] = str(max(0, int(tfbs[1]) - span1))
        tfbs[2] = str(int(tfbs[2]) + span2)
        modify_coordinate_list.append(tfbs)
    gene_df = pd.DataFrame(modify_coordinate_list)
    modify_coordinate = pybedtools.BedTool.from_dataframe(gene_df)
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
    - tfbs_of_tf_in_celltype (BedTool): 交集得到的 TFBS。
    - modify_coordinate_list (list): 修改后的坐标列表。
    - motif_len (int): motif 的长度。
    """
    tfbs_site = pybedtools.BedTool(tfbs_file)
    celltype = pybedtools.BedTool(chip_seq_of_celltype)
    tfbs_of_tf_in_celltype = tfbs_site.intersect(celltype, wa=True)

    fasta = pybedtools.BedTool(os.path.join(dir_path, '../genome', f'{species}.fa'))
    TFBSs_fasta = tfbs_of_tf_in_celltype.sequence(fi=fasta, s=True, fo='myc')
    with open(TFBSs_fasta.seqfn) as f:
        TFBSs_fasta_lines = f.readlines()

    motif_len = None
    for line in TFBSs_fasta_lines:
        if line.startswith(">"):
            continue
        else:
            motif_len = len(line.strip())
            break

    print('\nBedTool intersection finished...')
    print(f'{len(tfbs_of_tf_in_celltype)} transcription factor binding sites discovered')

    modify_coordinate_list, modify_coordinate = modify_bed(tfbs_of_tf_in_celltype, span1, span2)

    return tfbs_of_tf_in_celltype, modify_coordinate_list, modify_coordinate, motif_len


def get_seq(species, seq_list, tfbs_modify_coordinate, span1, span2, motif_len):
    """
    获取序列，并构建 total 数据结构。

    参数：
    - species (str): 物种名称。
    - seq_list (list): 修改后的坐标列表。
    - tfbs_modify_coordinate (BedTool): 修改后的位置信息。
    - span1 (int): 左侧延伸的碱基数。
    - span2 (int): 右侧延伸的碱基数。
    - motif_len (int): motif 的长度。

    返回：
    - total (defaultdict): 存储序列信息的字典。
    - seq_file (str): 序列文件路径。
    """
    total = defaultdict(init_list)
    fasta = pybedtools.BedTool(os.path.join(dir_path, '../genome', f'{species}.fa'))
    tfbss_fa = tfbs_modify_coordinate.sequence(fi=fasta, s=True)
    with open(tfbss_fa.seqfn) as f:
        tfbss_fasta = f.readlines()

    count = 0
    for i in range(len(tfbss_fasta)):
        if tfbss_fasta[i].startswith('>'):
            chr = seq_list[count][0]
            tfbsstartpos = int(seq_list[count][1]) + span1
            tfbsendpos = int(seq_list[count][2]) - span2
            strand = tfbss_fasta[i][-3]
            key = f'>{chr}:{tfbsstartpos}-{tfbsendpos}{strand}'
            count += 1
            continue
        else:
            total[key].append([0] * motif_len)  # ctx
            total[key].append([0] * motif_len)  # tread
            total[key].append([0] * motif_len)  # cread
            sequence = tfbss_fasta[i].strip()
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
    return total, tfbss_fa.seqfn


def methylread_counter(TFBSFile, WGBSFile, total):
    """
    计算每个 TFBS 的甲基化信息。

    参数：
    - TFBSFile (BedTool): TFBS 文件。
    - WGBSFile (list): WGBS 文件列表。
    - total (defaultdict): 存储序列信息的字典。

    返回：
    - ctx_dict (DataFrame): 上下文信息。
    - cread_dict (DataFrame): 甲基化的读计数。
    - tread_dict (DataFrame): 未甲基化的读计数。
    """
    print('\nStarting methylread_counter...')
    count = 1
    TFBSs = TFBSFile
    WGBSFile = [filepath for sublist in WGBSFile for filepath in sublist]

    for file in WGBSFile:
        WGBS = pybedtools.BedTool(file)
        start = time.time()
        Methylation_of_TFBSs = WGBS.intersect(TFBSs, wa=True, wb=True)
        end = time.time()
        print(f'Intersection finished, cost {(end - start) / 60:.2f} min')
        start = time.time()
        with open(Methylation_of_TFBSs.fn) as f:
            readMethylofTFBSsFile = f.readlines()
        print(f'File {count} processing finished...')

        for line in readMethylofTFBSsFile:
            fields = line.strip().split()
            tfbsstartpos = fields[12]
            tfbsendpos = fields[13]
            chr = fields[0]
            base_MethylationInfo_pos = fields[1]
            capped_read = int(fields[9])
            methyl_read = capped_read * int(fields[10]) // 100
            strand = fields[-1]

            key = f'>{chr}:{tfbsstartpos}-{tfbsendpos}{strand}'
            if key not in total:
                continue  # 如果键不在 total 中，跳过

            if strand == '+':
                # tread
                temp = total[key][1]
                idx = int(base_MethylationInfo_pos) - int(tfbsstartpos)
                temp[idx] += capped_read - methyl_read
                # cread
                temp1 = total[key][2]
                temp1[idx] += methyl_read
            else:
                try:
                    temp = total[key][1]
                    idx = int(tfbsstartpos) - int(base_MethylationInfo_pos) - 1
                    temp[idx] += capped_read - methyl_read
                    # cread
                    temp1 = total[key][2]
                    temp1[idx] += methyl_read
                except IndexError:
                    print("IndexError occurred.")
                    pass
        end = time.time()
        count += 1

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
