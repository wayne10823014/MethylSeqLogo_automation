import os
from Bio import SeqIO
import pandas as pd
from decimal import *
import numpy as np
import matplotlib
import sys
import unittest
import json
from jinja2 import Environment, FileSystemLoader  # 新增

# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects
from package.arg import get_parser
from package.datapreprocess import get_tfbs, get_seq, methylread_counter, flanking_bed
from package.backgroundprocess import read_bgprob_table, promoter, isfile, neighbor
from package.calculate import *
from package.figure_setting import *
import time

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid


def remove_dup_seqs(records):
    """"SeqRecord iterator to removing duplicate sequences."""
    checksums = set()
    for record in records:
        checksum = seguid(record.seq)
        if checksum in checksums:
            # print("Ignoring %s" % record.id)
            continue
        checksums.add(checksum)
        yield record


def main():
    filestart = time.time()
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

    print("\n")
    if logotype in ['Shannon', 'Kullback-Liebler']:
        print("Plotting " + TF + ' ' + mode + ' ' + logotype + ' logo of ' + species + ' ' + celltype + ' ' + region + ' background.')
    elif logotype == 'riverlake':
        print("Plotting " + TF + ' ' + mode + ' ' + logotype + ' of ' + species + ' ' + celltype + '.')
    else:
        print("Plotting " + TF + ' ' + mode + ' all logos of ' + species + ' ' + celltype + '.')
    print("\n")

    # span = 50 if region == "promoter" else 2
    spanL = int(region) if region.isdigit() else 2
    spanR = spanL

    jaspar_file = args.jaspar
    remap_file = args.remap
    tfbs_bed, tfbs = get_tfbs(remap_file, jaspar_file, species, spanL, spanR)

    seqdata = get_seq(species, tfbs[0], tfbs[1], spanL, spanR)

    try:
        seqdata = SeqIO.to_dict(SeqIO.parse(seqdata, 'fasta'))
    except:
        # print("SeqRecord iterator to removing duplicate sequences.")
        seqdata = remove_dup_seqs(SeqIO.parse(seqdata, 'fasta'))
        seqdata = SeqIO.to_dict(seqdata)

    seqdata = pd.DataFrame.from_dict(seqdata, orient='index')

    seqdata.reset_index(drop=True)
    # seqdata.to_csv('seq.csv')
    print(len(seqdata))
    start = time.time()

    # 定义输出目录，仍然是 Output1 目录
    output_dir = "/home/wayne/MethylSeqLogo_automation/Output1/"
    os.makedirs(output_dir, exist_ok=True)  # 如果目录不存在，创建它

    # 文件名的前缀
    logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + logotype + '_'

    # 使用 os.path.join 将路径和文件名连接起来
    ctx_filepath = os.path.join(output_dir, logoname + "ctx.csv")
    cread_filepath = os.path.join(output_dir, logoname + "cread.csv")
    tread_filepath = os.path.join(output_dir, logoname + "tread.csv")

    # 检查文件是否存在
    if os.path.isfile(ctx_filepath) and os.path.isfile(cread_filepath) and os.path.isfile(tread_filepath):
        ctxdata = pd.read_csv(ctx_filepath, sep='\t')
        creaddata = pd.read_csv(cread_filepath, sep='\t')
        treaddata = pd.read_csv(tread_filepath, sep='\t')
    else:
        ctxdata, creaddata, treaddata = methylread_counter(tfbs_bed, methylbed)

    end = time.time()
    print('\nmethylatedread_counter finished...,total cost', (end - start) // 60, 'min\n')

    motif_len = len(seqdata.columns) - (spanL + spanR)
    print(TF + " binding motif is " + str(motif_len) + "bp")

    if plotlen is None:
        plotlen = motif_len - beginningoftfbs + 1
    else:
        pass

    global endpos
    if (beginningoftfbs + plotlen > motif_len):
        endpos = motif_len
        print("Warning: user defined plotting length is over motif length" + '\n')
        print("Plotting " + TF + " binding motif from pos " + str(beginningoftfbs) + " to pos " + str(motif_len))
    else:
        endpos = beginningoftfbs + plotlen - 1
        print("Plotting " + TF + " binding motif from pos " + str(beginningoftfbs) + " to pos " + str(endpos))

    print("\n")
    seqdata = seqdata.iloc[:, beginningoftfbs - 1 + spanL: beginningoftfbs - 1 + spanL + motif_len]
    ctxdata = ctxdata.iloc[:, beginningoftfbs - 1:endpos]
    creaddata = creaddata.iloc[:, beginningoftfbs - 1: endpos]
    treaddata = treaddata.iloc[:, beginningoftfbs - 1: endpos]

    # 将数据存成 CSV 文件
    ctxdata.to_csv(ctx_filepath, sep='\t', index=False)
    creaddata.to_csv(cread_filepath, sep='\t', index=False)
    treaddata.to_csv(tread_filepath, sep='\t', index=False)

    global pseudocount
    pseudocount = 1.0

    # Check species.fa exists
    isfile(dir_path + '/genome/' + species + '.fa')

    if region == 'whole_genome':
        bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
    elif (region.isdigit()):
        print('~~neighbor~~')
        path = dir_path + "/../Background_probability/neighbor/" + species + '_' + TF + '_' + celltype + '_' + region + '_probability.txt'
        path1 = dir_path + "/../Background_probability/neighbor/" + species + '_' + TF + '_' + celltype + '_' + region + '_methyl_probability.txt'
        if os.path.isfile(path) and os.path.isfile(path1):
            print('~~find neighbor~~')
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region, TF)
        else:
            print('~~not find neighbor~~')
            flankingbed = flanking_bed(tfbs_bed, spanL, spanR)
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = neighbor(flankingbed, species, methylbed, celltype, region, TF)

    else:
        print('~~promoter~~')
        start = time.time()
        path = dir_path + "/../Background_probability/promoter/" + species + '_' + celltype + '_' + region + '_probability.txt'
        path1 = dir_path + "/../Background_probability/promoter/" + species + '_' + celltype + '_' + region + '_methyl_probability.txt'
        if os.path.isfile(path) and os.path.isfile(path1):
            print('~~find promoter~~')
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
        else:
            print('~~not find promoter~~')
            bgpps, bg_mCG, bg_mCHG, bg_mCHH = promoter(tfbs_bed, species, methylbed, celltype, region)
            end = time.time()
            print('promoter bg calc finished...,total cost', (end - start) // 60, 'min\n')

    ppm = calc_ppm(seqdata, TF, species, celltype, region)

    # 定义 ppm.json 的文件名
    ppmname = f"{TF}_{species}_{celltype}_{region}_ppm.json"

    # 定义保存路径
    ppm_filepath = os.path.join(output_dir, ppmname)

    # 保存 ppm.json
    ppm_list = ppm[['A', 'C', 'G', 'T']].values.tolist()
    with open(ppm_filepath, 'w') as f:
        json.dump(ppm_list, f)
    print("PPM 数据已保存到 ppm.json 文件中，路径为：", ppm_filepath)

    # 使用 Jinja2 渲染模板，生成 home.html
    template_dir = os.path.join(dir_path, 'package')  # 模板目录，即 home_template.html 所在目录
    print("dir_path:", dir_path)
    print("template_dir:", template_dir)
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('home_template.html')

    # 计算从 home.html 到 ppm.json 的相对路径
    # home.html 位于 bin/package/ 目录
    # ppm.json 位于 Output1/ 目录
    # 因此，相对路径为 '../../Output1/ppmname'

    relative_ppm_path = os.path.relpath(ppm_filepath, start=os.path.join(dir_path, 'bin', 'package'))

    # 渲染模板，传入 ppmname
    rendered_html = template.render(ppmname=relative_ppm_path)

    # 保存生成的 home.html，覆盖原有的 home.html
    home_html_path = os.path.join(dir_path, 'package', 'home.html')
    print("home_html_path:", home_html_path)
    with open(home_html_path, 'w') as f:
        f.write(rendered_html)
    print("home.html 已生成，路径为：", home_html_path)

    # 继续您的代码...

    # main()
if __name__ == '__main__':
    main()
