import os
from Bio import SeqIO
import pandas as pd
from decimal import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects

from package.arg import get_parser
from package.datapreprocess import  get_tfbs,get_seq,methylread_counter
from package.backgroundprocess import read_bgprob_table
from package.calculate import *
from package.figure_setting import *
import time


def main():
    

    # logging.basicConfig(filename="methylseqlogo.log", level = logging.DEBUG, format = '%(asctime)s, %(name)s-%(levelname)-s: %(message)s')
    parser = get_parser()
    args = parser.parse_args()

    global dir_path
    dir_path = os.path.dirname(os.path.realpath(__file__))
    print (dir_path)

    global TF
    TF = args.transcriptionfactor  

    global species
    species = args.species

    global celltype
    celltype = args.celltypes

    global mode
    mode = args.mode

    global logotype
    logotype = args.logotype

    global region
    region = args.regionofinterest

    global threshold
    threshold = float(args.threshold)

    global beginningoftfbs  
    beginningoftfbs = args.beginningoftfbs

    global plotlen
    plotlen = args.plotlen

    jaspar_file = args.jaspar
    remap_file = args.remap


    # print dir_path
    pd.set_option('display.float_format', lambda x: '%.2f' % x)


    print ("\n")
    if logotype in ['Shannon', 'Kullback-Liebler']:
        print ("Plotting " + TF + ' ' + mode + ' ' + logotype + ' logo of ' + species + ' ' + celltype + ' ' + region + ' background.')
    elif logotype == 'riverlake':
        print ("Plotting " + TF + ' ' + mode + ' ' + logotype + ' of ' + species + ' ' + celltype + '.')
    else:
        print ("Plotting " + TF + ' ' + mode + ' all logos of ' + species + ' ' + celltype + '.')
    print ("\n")


 
#***************************
    
    

    

    # methylbed = args.methylationinfo
    
    # tfbs_bed,tfbs = get_tfbs(remap_file, jaspar_file,dir_path)
    
    # seqdata = get_methyl_context(dir_path, species, tfbs[0], tfbs[1])
    
    # seqdata= SeqIO.to_dict(SeqIO.parse(seqdata, 'fasta'))
    # seqdata= pd.DataFrame.from_dict(seqdata, orient= 'index')
    # seqdata.reset_index(drop= True)

    # start = time.time()
    # ctxdata,creaddata, treaddata = methylatedread_counter(tfbs_bed, methylbed)
    # end = time.time()
    # print('\nmethylatedread_counter finished...,total cost', (end-start)//60, 'min\n')
    # # print(creaddata)
#***************************

    seqdata= SeqIO.to_dict(SeqIO.parse('/home/yuling/methylseqlogov2/Input/MYC_H1-hESC_whole_genome_binding_sites_seq.fa', 'fasta'))
    seqdata= pd.DataFrame.from_dict(seqdata, orient= 'index')
    seqdata.reset_index(drop= True)
    ctxdata= pd.read_csv('/home/yuling/methylseqlogov2/Input/MYC_H1-hESC_promoter_binding_sites_ctx.txt', sep= "\t")
    creaddata= pd.read_csv('/home/yuling/methylseqlogov2/Input/MYC_H1-hESC_whole_genome_binding_sites_cread.txt', sep= "\t")
    treaddata= pd.read_csv('/home/yuling/methylseqlogov2/Input/MYC_H1-hESC_promoter_binding_sites_tread.txt', sep= "\t")    

    # Determine plotting window
    # Default setting (0, motiflen)
    global motif_len
    motif_len = len(seqdata.columns)-4+4
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
    # seqdata = seqdata.iloc[:, spanwindowleft-1+2:endpos+2]
    seqdata = seqdata.iloc[:, beginningoftfbs-1:endpos]
    ctxdata = ctxdata.iloc[:, beginningoftfbs-1:endpos]
    creaddata = creaddata.iloc[:, beginningoftfbs-1:endpos]
    treaddata = treaddata.iloc[:, beginningoftfbs-1:endpos]

    global pseudocount
    pseudocount = 1.0

    ppm = calc_ppm(seqdata, TF, species, celltype, region,)
    # JiCs, PiCs, Cmethyls, Gmethyls = MethylProbability(ctxdata, creaddata, treaddata, J_bCG, J_bCHG, J_bCHH)

    bgpps,  bg_mCG, bg_mCHG, bg_mCHH = read_bgprob_table(species, celltype, region)
    JiCs, PiCs, Cmethyls, Gmethyls, Freqs_ = calc_methylprob(ctxdata, creaddata, treaddata, bg_mCG, bg_mCHG, bg_mCHH, plotlen)

    if logotype in ['Kullback-Liebler', 'Shannon']:
        Cents = calc_methylation_entropy(JiCs, PiCs, bg_mCG, bg_mCHG, bg_mCHH, logotype)
        entropys = calc_totalEntropy(ppm, bgpps, Cents, logotype, plotlen, TF, species, celltype, region)
        four_base_heights = calc_perBaseEntropy(entropys, ppm, mode, TF, species, celltype, region)
        dippm = to4basedippm(seqdata, plotlen)
        dientropys, bg_dientropys_max, bg_dientropys_min = twomerBg(bgpps, dippm, plotlen)
        fig = set_fig(entropys)
        plotobj = seqLogoPlot(fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls, bgpps, dientropys, bg_dientropys_max, bg_dientropys_min, bg_mCG, bg_mCHG, bg_mCHH, Freqs_)
        plotobj.plotlogo()
        logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + logotype + '_seqlogo.png'
        plt.savefig(dir_path + '/../Output/' + logoname, bbox_inches = 'tight', dpi = 600)
        print (logoname + ' is saved in' + dir_path + '/../Output/.')
    # elif logotype == 'riverlake':
    #     dippm = to4basedippm(seqdata)
    #     Cents = methylationEntropy(JiCs, PiCs, J_bCG, J_bCHG, J_bCHH, logotype)
    #     entropys = totalEntropy(ppm, bgpps, Cents, logotype)
    #     four_base_heights = perBaseEntropy(entropys, ppm)
    #     # fig = set_fig(entropys)
    #     fig = plt.figure(figsize = (plotlen+1, 3.0))
    #     riverlakeobj = riverLake(fig, celltype, ppm, dippm, Cmethyls, Gmethyls, bgpps, J_bCG, J_bCHG, J_bCHH, Freqs_)
    #     riverlakeobj.plotRiverLake()
    #     logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + logotype + '_seqlogo_bar7.pdf'
    #     plt.savefig(dir_path + '/../Output/' + logoname, bbox_inches = 'tight', dpi = 600)
    #     print (logoname + ' is saved in' + dir_path + '/../Output/.')
    # elif logotype == 'all':
    #     for i in ['Kullback-Liebler', 'Shannon']:
    #         # Cents = methylationEntropy(JiCs, PiCs, J_bCG, J_bCHG, J_bCHH, logotype)
    #         Cents = methylationEntropy(JiCs, PiCs, J_bCG, J_bCHG, J_bCHH, i)
    #         entropys = totalEntropy(ppm, bgpps, Cents, i)
    #         # entropys = totalEntropy(ppm, bgpps, Cents, logotype)
    #         four_base_heights = perBaseEntropy(entropys, ppm)
    #         fig = set_fig(entropys)
    #         plotobj = seqLogoPlot(fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls, bgpps, J_bCG, J_bCHG, J_bCHH, Freqs_)
    #         plotobj.plotlogo()
    #         logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + i + '_seqlogo_bar7.pdf'
    #         plt.savefig(dir_path + '/../Output/' + logoname, bbox_inches = 'tight', dpi = 600)
    #         print (logoname + ' is saved in ./Output/.')
    #     dippm = to4basedippm(seqdata)
    #     dientropys = twomerBg(bgpps, dippm)
    #     Cents = methylationEntropy(JiCs, PiCs, J_bCG, J_bCHG, J_bCHH, 'riverlake')
    #     entropys = totalEntropy(ppm, bgpps, Cents, 'riverlake')
    #     four_base_heights = perBaseEntropy(entropys, ppm)
    #     # fig = set_fig(entropys)

    #     fig = plt.figure(figsize = (plotlen, 3.0))
    #     riverlakeobj = riverLake(fig, celltype, ppm, dippm, Cmethyls, Gmethyls, bgpps, J_bCG, J_bCHG, J_bCHH, Freqs_)
    #     riverlakeobj.plotRiverLake()
    #     logoname = TF + '_' + species + '_' + celltype + '_' + region + '_' + mode + '_' + 'riverlake' + '_seqlogo_bar7.png'
    #     plt.savefig(dir_path + '/../Output/' + logoname, bbox_inches = 'tight', dpi = 600)
    #     print (logoname + ' is saved in ./Output/.')	
    # fileend = time.time()
    # print('total cost:', (fileend-filestart)/60, 'min')
    # print ("\n")

main()
# if __name__ == '__main__':
# main()