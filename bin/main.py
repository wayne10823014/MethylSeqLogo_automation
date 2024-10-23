import os
from Bio import SeqIO
import pandas as pd
from decimal import *
import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects
from package.arg import get_parser
from package.datapreprocess import  get_tfbs,get_seq, methylread_counter, flanking_bed
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
            # print ("Ignoring %s" % record.id)
            continue
        checksums.add(checksum)
        yield record

def main():
    filestart= time.time()
    # logging.basicConfig(filename="methylseqlogo.log", level = logging.DEBUG, format = '%(asctime)s, %(name)s-%(levelname)-s: %(message)s')
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


    print ("\n")
    if logotype in ['Shannon', 'Kullback-Liebler']:
        print ("Plotting " + TF + ' ' + mode + ' ' + logotype + ' logo of ' + species + ' ' + celltype + ' ' + region + ' background.')
    elif logotype == 'riverlake':
        print ("Plotting " + TF + ' ' + mode + ' ' + logotype + ' of ' + species + ' ' + celltype + '.')
    else:
        print ("Plotting " + TF + ' ' + mode + ' all logos of ' + species + ' ' + celltype + '.')
    print ("\n")
    


    # span = 50 if region == "promoter" else 2 
    spanL = int(region) if region.isdigit()  else 2
    spanR = spanL

    jaspar_file= args.jaspar
    remap_file = args.remap
    tfbs_bed, tfbs = get_tfbs(remap_file, jaspar_file, species, spanL, spanR)


    
    seqdata = get_seq(species, tfbs[0], tfbs[1],spanL, spanR)
    
    try:
        seqdata = SeqIO.to_dict(SeqIO.parse(seqdata, 'fasta'))
    except:
        # print("SeqRecord iterator to removing duplicate sequences.")
        seqdata = remove_dup_seqs(SeqIO.parse(seqdata, 'fasta'))
        seqdata = SeqIO.to_dict(seqdata)
    
    seqdata= pd.DataFrame.from_dict(seqdata, orient= 'index')
    
    seqdata.reset_index(drop= True)
    # seqdata.to_csv( 'seq.csv')
    print(len(seqdata))  
    start = time.time()

    output_dir = "/home/wayne/MethylSeqLogo_automation/Output1/"

    logoname = output_dir + TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + logotype + '_'

    ctx_file = logoname + 'ctx.csv'
    cread_file = logoname + 'cread.csv'
    tread_file = logoname + 'tread.csv'

    print(ctx_file)

    # 檢查文件是否存在
    if os.path.isfile(ctx_file) and os.path.isfile(cread_file) and os.path.isfile(tread_file):
        ctxdata = pd.read_csv(ctx_file, sep='\t')
        creaddata = pd.read_csv(cread_file, sep='\t')
        treaddata = pd.read_csv(tread_file, sep='\t')
    else:
        ctxdata, creaddata, treaddata = methylread_counter(tfbs_bed, methylbed)
    
    end = time.time()
    print('\nmethylatedread_counter finished...,total cost', (end-start)//60, 'min\n')


    motif_len = len(seqdata.columns) - (spanL+spanR)
    print (TF + " binding motif is " + str(motif_len) + "bp")

    if plotlen is None:
        plotlen = motif_len - beginningoftfbs + 1
    else:
        pass

    global endpos
    if (beginningoftfbs + plotlen > motif_len):
        endpos = motif_len
        print ("Warning: user defined plotting length is over motif length" + '\n')
        print ("Plotting " + TF + " binding motif from pos " + str(beginningoftfbs) + " to pos " + str(motif_len))
    else:
        endpos = beginningoftfbs + plotlen - 1
        print ("Plotting " + TF + " binding motif from pos " + str(beginningoftfbs) + " to pos " + str(endpos))

  
    print("\n")
    seqdata = seqdata.iloc[:, beginningoftfbs - 1 + spanL : beginningoftfbs - 1 + spanL + motif_len]
    # print(seqdata)
    ctxdata = ctxdata.iloc[:, beginningoftfbs - 1:endpos]
    # print(ctxdata)
    creaddata = creaddata.iloc[:, beginningoftfbs - 1 : endpos]
    treaddata = treaddata.iloc[:, beginningoftfbs - 1 : endpos]
    # print(creaddata)
    # print(treaddata)

    output_dir = "/home/wayne/MethylSeqLogo_automation/Output1"

    # 確保目錄存在，如果不存在則創建
    os.makedirs(output_dir, exist_ok=True)

    # 文件名的前綴
    logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + logotype + '_'

    # 使用 os.path.join 將路徑和文件名連接起來
    ctx_filepath = os.path.join(output_dir, logoname + "ctx.csv")
    cread_filepath = os.path.join(output_dir, logoname + "cread.csv")
    tread_filepath = os.path.join(output_dir, logoname + "tread.csv")

    # 將數據存成 CSV 文件
    ctxdata.to_csv(ctx_filepath, sep='\t', index=False)
    creaddata.to_csv(cread_filepath, sep='\t', index=False)
    treaddata.to_csv(tread_filepath, sep='\t', index=False)

    global pseudocount
    pseudocount = 1.0

      
    #Check species.fa exits
    isfile(dir_path + '/genome/' + species +'.fa')

    if region == 'whole_genome':
        bgpps, bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
    elif (region.isdigit()):
       print('~~neighbor~~')
       path = dir_path + "/../Background_probability/neighbor/"+ species + '_' + TF  +'_'+ celltype + '_' + region + '_probability.txt'
       path1 = dir_path + "/../Background_probability/neighbor/" + species + '_' + TF + '_' + celltype + '_' + region  +'_methyl_probability.txt'
       if os.path.isfile(path) and os.path.isfile(path1):
            print('~~find neighbor~~')
            bgpps,  bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region, TF)
       else:
            print('~~not find neighbor~~')
            flankingbed = flanking_bed(tfbs_bed, spanL, spanR)
            bgpps,  bg_mCG, bg_mCHG, bg_mCHH =  neighbor(flankingbed, species, methylbed, celltype, region, TF)

    else:
       print('~~promoter~~')
       start = time.time()
       path = dir_path + "/../Background_probability/promoter/"+ species + '_' + celltype + '_' + region + '_probability.txt'
       path1 = dir_path + "/../Background_probability/promoter/"+ species + '_' + celltype + '_' + region +'_methyl_probability.txt'
       if os.path.isfile(path) and os.path.isfile(path1):
            print('~~find promoter~~')
            bgpps,  bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
       else:
            print('~~not find promoter~~')
            bgpps,  bg_mCG, bg_mCHG, bg_mCHH =  promoter(tfbs_bed, species, methylbed, celltype, region)
            end = time.time()
            print('promoter bg calc finished...,total cost', (end-start)//60, 'min\n')


    ppm = calc_ppm(seqdata, TF, species, celltype, region,)
    C_ratio, G_ratio, Cmethyls, Gmethyls, Freqs_ = calc_methylprob(ctxdata, creaddata, treaddata, bg_mCG, bg_mCHG, bg_mCHH, plotlen)

    if logotype in ['Kullback-Liebler', 'Shannon']:
        Cents = calc_methylation_entropy( C_ratio, G_ratio, Cmethyls, Gmethyls, bg_mCG, bg_mCHG, bg_mCHH, logotype)
        entropys = calc_totalEntropy(ppm, bgpps, Cents, logotype, plotlen, TF, species, celltype, region)
        four_base_heights = calc_perBaseEntropy(entropys, ppm, mode, TF, species, celltype, region)
        dippm = to4basedippm(seqdata, plotlen)
        dientropys, bg_dientropys_max, bg_dientropys_min = twomerBg(bgpps, dippm, plotlen)
        fig = set_fig(entropys, logotype, mode, plotlen)
        plotobj = seqLogoPlot(fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls, bgpps, dientropys, bg_dientropys_max, bg_dientropys_min, bg_mCG, bg_mCHG, bg_mCHH, Freqs_, mode, plotlen, threshold, TF)
        plotobj.plotlogo()
        logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + logotype + '_seqlogo.png'
        plt.savefig(dir_path + '/../Output1/' + logoname, bbox_inches = 'tight', dpi = 600)
        print (logoname + ' is saved in' + dir_path + '/../Output1/.')
    elif logotype == 'riverlake':
        dippm = to4basedippm(seqdata)
        Cents = calc_methylation_entropy(C_ratio, G_ratio, bg_mCG, bg_mCHG, bg_mCHH, logotype)
        entropys = calc_totalEntropy(ppm, bgpps, Cents, logotype)
        four_base_heights = calc_perBaseEntropy(entropys, ppm)
        # fig = set_fig(entropys)
        fig = plt.figure(figsize = (plotlen+1, 3.0))
        riverlakeobj = riverLake(fig, celltype, ppm, dippm, Cmethyls, Gmethyls, bgpps, bg_mCG, bg_mCHG, bg_mCHH, Freqs_)
        riverlakeobj.plotRiverLake()
        logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + logotype + '_seqlogo_bar7.pdf'
        plt.savefig(dir_path + '/../Output1/' + logoname, bbox_inches = 'tight', dpi = 600)
        print (logoname + ' is saved in' + dir_path + '/../Output1/.')
    elif logotype == 'all':
        for i in ['Kullback-Liebler', 'Shannon']:
            # Cents = methylationEntropy(JiCs, PiCs, J_bCG, J_bCHG, J_bCHH, logotype)
            Cents = calc_methylation_entropy(C_ratio, G_ratio, bg_mCG, bg_mCHG, bg_mCHH, i)
            entropys = calc_totalEntropy(ppm, bgpps, Cents, i)
            # entropys = totalEntropy(ppm, bgpps, Cents, logotype)
            four_base_heights = calc_perBaseEntropy(entropys, ppm)
            fig = set_fig(entropys)
            plotobj = seqLogoPlot(fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls, bgpps, bg_mCG, bg_mCHG, bg_mCHH, Freqs_)
            plotobj.plotlogo()
            logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + i + '_seqlogo_bar7.pdf'
            plt.savefig(dir_path + '/../Output1/' + logoname, bbox_inches = 'tight', dpi = 600)
            print (logoname + ' is saved in ./Output1/.')
        dippm = to4basedippm(seqdata)
        dientropys = twomerBg(bgpps, dippm)
        Cents = calc_methylation_entropy(C_ratio, G_ratio, bg_mCG, bg_mCHG, bg_mCHH, 'riverlake')
        entropys = calc_totalEntropy(ppm, bgpps, Cents, 'riverlake')
        four_base_heights = calc_perBaseEntropy(entropys, ppm)
        # fig = set_fig(entropys)

        fig = plt.figure(figsize = (plotlen, 3.0))
        riverlakeobj = riverLake(fig, celltype, ppm, dippm, Cmethyls, Gmethyls, bgpps, bg_mCG, bg_mCHG, bg_mCHH, Freqs_)
        riverlakeobj.plotRiverLake()
        logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + 'riverlake' + '_seqlogo_bar7.png'
        plt.savefig(dir_path + '/../Output1/' + logoname, bbox_inches = 'tight', dpi = 600)
        print (logoname + ' is saved in ./Output1/.')	
    fileend = time.time()
    print('total cost:', (fileend-filestart)/60, 'min')
    print ("\n")

# main()
if __name__ == '__main__':
    main()
