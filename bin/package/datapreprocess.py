
import pybedtools 
from Bio import SeqIO
import pandas as pd
from collections import defaultdict
import time,os


def init_list():
    return []

global total
total = defaultdict(init_list)

global motif_len
motif_len = 0

dir_path = os.path.dirname(os.path.realpath(__file__))
path = dir_path + "/../../Output/"

def flanking_bed(tfbs_of_tf_in_celltype, span1, span2):
    
    modify_coordinate_list = []
    temp = []
    for tfbs in tfbs_of_tf_in_celltype:
        temp = []
        tfbs = list(tfbs)
        temp.append(tfbs[0])
        start = str(max(0, int(tfbs[1]) - span1))
        temp.append(start)
        temp.append(tfbs[1])
        modify_coordinate_list.append(temp)

        temp = []
        temp.append(tfbs[0])
        temp.append(tfbs[2])
        temp.append(str(int(tfbs[2]) + span2))
        modify_coordinate_list.append(temp)

    
    gene_df = pd.DataFrame(modify_coordinate_list)
    modify_coordinate = pybedtools.BedTool.from_dataframe(gene_df)
    # gene_df.to_csv(dir_path + '/temp/' + 'modify_coordinate.bed', sep = '\t', index = False,header=False)
    print('modify',len(modify_coordinate_list),len(modify_coordinate))
    return modify_coordinate

def modify_bed(tfbs_of_tf_in_celltype, span1, span2):
    """In order to get context give C base,modify coordinate of TFBSs in bed

    #### Args:
        - tfbs_of_tf_in_celltype(_BedToolObject_) : TFBSs's bed

    #### Returns:
        - modify_coordinate_list (_list_) 
        - modify_coordinate (_BedToolObject_) 
    """
    modify_coordinate_list = []
    for tfbs in tfbs_of_tf_in_celltype:
        tfbs = list(tfbs)
        tfbs[1] = str(max(0, int(tfbs[1]) - span1))
        tfbs[2] = str(int(tfbs[2]) + span2)
        modify_coordinate_list.append(tfbs)
    gene_df = pd.DataFrame(modify_coordinate_list)
    modify_coordinate = pybedtools.BedTool.from_dataframe(gene_df)
    # gene_df.to_csv(dir_path + '/temp/' + 'modify_coordinate.bed', sep = '\t', index = False,header=False)
    print('modify',len(modify_coordinate_list),len(modify_coordinate))
    return modify_coordinate_list, modify_coordinate

def get_tfbs(chip_seq_of_celltype, tfbs, species, span1,span2):
    global motif_len

    """ Find TFBss( motif ) of TF in celltype \n
    (two bed file must use the same genome assembly)

    #### Args:
      - chip_seq_of_celltype (_.bed file_) : chip_seq data for TF in celltype
      - tfbs (_.bed file_) : TFBSs( motif ) of TF 

    #### Returns:
        - tfbs_of_tf_in_celltype (_BedToolObject_) : intersect output -> TFBSs( motif ) of TF in celltype
        - modify_bedfile(tfbs_of_tf_in_celltype) (_func_) :
 
    """

 
    tfbs_site = pybedtools.BedTool(tfbs)  
    celltype = pybedtools.BedTool(chip_seq_of_celltype)
    tfbs_of_tf_in_celltype = tfbs_site.intersect(celltype, wa = True)

    fasta = pybedtools.BedTool(dir_path + '/../genome/'+species+'.fa')
    TFBSs_fasta  = tfbs_of_tf_in_celltype.sequence(fi=fasta,s=True,fo='myc')
    TFBSs_fasta = open(TFBSs_fasta.seqfn).readlines()
    txt = []
    for i in TFBSs_fasta:
        if i.startswith(">"):
            key = i.replace('\n','')
            keys = key.split('(')
            key = keys[0] 
            continue

        else:
            tfbs = list(i)
            txt.append(tfbs[:-1])
            break
    motif_len = len(txt[0])
    
    
    print('\n'+'bedtool intersection finish...')
    print(len(tfbs_of_tf_in_celltype),' transcription factor binding sites discovered')

    return tfbs_of_tf_in_celltype, modify_bed(tfbs_of_tf_in_celltype, span1, span2)


    
def get_seq(species, seq, tfbs_modify_coordinate, span1, span2):   
    tfbss = pybedtools.BedTool(tfbs_modify_coordinate)
    fasta = pybedtools.BedTool(dir_path + '/../genome/'+species+'.fa')
    tfbss_fa  = tfbss.sequence(fi=fasta,s=True)
    tfbss_fasta = open(tfbss_fa.seqfn).readlines()
    count = 0
    for rows in range(len(tfbss_fasta)):
      
        if tfbss_fasta[rows].startswith('>'):
           
            chr = seq[count][0]
            tfbsstartpos = int(seq[count][1]) + span1
            tfbsendpos  = int(seq[count][2]) - span2
            strand = tfbss_fasta[rows][-3]
            key = '>'+str(chr)+':'+str(tfbsstartpos)+'-'+str(tfbsendpos)+strand
            count += 1 
            continue

        else:
            
            total[key].append([0]*motif_len)
            total[key].append([0]*motif_len)
            total[key].append([0]*motif_len)
            for base in range(span1, span1+motif_len, 1):
                if tfbss_fasta[rows][base] not in ["c", "C", "g", "G"]:
                    continue
                else:
                    temp = total[key][0]
                    if strand == '+' and tfbss_fasta[rows][base] in ["c", "C"]:
                        if tfbss_fasta[rows][base+1] in ["g","G"] :
                            temp[base-span1] = 'X'
                        #CHG condition | C base
                        elif tfbss_fasta[rows][base+2] in ["g","G"]:
                            temp[base-span1] = 'Y'
                        #CHH condition | C base    
                        else:
                            temp[base-span1] = 'Z'
                    elif strand == '+' and tfbss_fasta[rows][base] in ["G", "g"]:
                        #CG condition | G base
                        if tfbss_fasta[rows][base-1] in ["c", "C"] :
                            
                            temp[base-span1] = 'x'
                        elif tfbss_fasta[rows][base-2] in ["c", "C"]:
                            temp[base-span1] = 'y'
                        else:
                            temp[base-span1] = 'z'  

                    elif strand == '-' and tfbss_fasta[rows][base] in ["G", "g"]:
                        #CG condition | C base
                        if tfbss_fasta[rows][base-1] in ["c", "C"] :
                            
                            temp[base-span1] = 'x'
                        #CHG condition | C base
                        elif tfbss_fasta[rows][base-2] in ["c", "C"]:
                            temp[base-span1] = 'y'
                        #CHH condition | C base    
                        else:
                            temp[base-span1] = 'z'
                    elif strand == '-' and tfbss_fasta[rows][base] in ["c", "C"]:
                        #CG condition | G base
                        if tfbss_fasta[rows][base+1] in ["g","G"] :
                            
                            temp[base-span1] = 'X'
                        elif tfbss_fasta[rows][base+2] in ["g","G"]:
                            temp[base-span1] = 'Y'
                        else:
                            temp[base-span1] = 'Z'   

    return tfbss_fa.seqfn


def methylread_counter(TFBSFile, WGBSFile):
    """ Calc  methyl-info in each TFBS \n
    

    #### Args:
      - TFBSFile (_.bed file_) : tfbs of TF in celltype
      - WGBSFile (_.bed file_) : celltype with WGBS 
    #### Returns:
        - ctx_dict (_DataFrame_) : methyl-condition when given C/G base in each TFBS (x: CG, y: CHG, z: CHH; upper case: forward strand, lower case: reverse strand)
        - cread_dict(_DataFrame_)  : read count as cytosine/guanine of C/G (methylated read) in each TFBS
        - tread_dict(_DataFrame_) : read count as thymine/adenine of T/A (un-methylated read) in each TFBS
 
    """
    print('\nmethylread_counter starting... ')
    count = 1
    TFBSs = TFBSFile
    WGBSFile = [ filepath for file in  WGBSFile for filepath in file]
    

    for file in WGBSFile:
        WGBS = pybedtools.BedTool(file)
        start = time.time()
        Methylation_of_TFBSs= WGBS.intersect(TFBSs, wa = True, wb = True) 
        end = time.time()
        print('intersect finish,cost',(end-start)//60,'min')
        start = time.time()
        readMethylofTFBSsFile= open(Methylation_of_TFBSs.fn).readlines()
        print(str(count) + ' file finish......')
    

        for line in readMethylofTFBSsFile:
            Methylation_information_for_one_row = line.split()
            tfbsstartpos = Methylation_information_for_one_row[12]
            tfbsendpos = Methylation_information_for_one_row[13]
            chr = Methylation_information_for_one_row[0]
            base_MethylationInfo_pos = Methylation_information_for_one_row[1]
            capped_read = int(Methylation_information_for_one_row[9])
            methyl_read = capped_read * int(Methylation_information_for_one_row[10])//100
            strand = Methylation_information_for_one_row[-1]
    
            
            key = '>'+str(chr)+':'+str(tfbsstartpos)+'-'+str(tfbsendpos)+strand
            if strand == '+':
                #tread
                temp = total[key][1]
                temp[int(base_MethylationInfo_pos)-int(tfbsstartpos)] = temp[int(base_MethylationInfo_pos)-int(tfbsstartpos)] + capped_read- methyl_read
    
                #cread
                temp1 = total[key][2]
                temp1[int(base_MethylationInfo_pos)-int(tfbsstartpos)] = temp1[int(base_MethylationInfo_pos)-int(tfbsstartpos)] + methyl_read
            else:
                temp = total[key][1]
                key = -(int(base_MethylationInfo_pos)-int(tfbsstartpos))-1
                temp[key] = temp[key] + capped_read- methyl_read
    
                #cread
                temp1 = total[key][2]
                temp1[key] = temp1[key] + methyl_read
        end = time.time()

        count += 1

    ctx = []
    tread = []
    cread= []
    kKey = list(total.keys())
    print(kKey[0:6])
    for k in total:
        try:

            ctx.append(total[k][0])
            tread.append(total[k][1])
            cread.append(total[k][2])
        except:
            pass
    
    ctx_dict = pd.DataFrame(ctx)
    ctx_dict.to_csv("MYC_H1_whole_genome_binding_sites_ctx.csv", sep = '\t', index = False) 
    cread_dict = pd.DataFrame(cread)
    cread_dict.to_csv("MYC_H1_whole_genome_binding_sites_cread.csv", sep = '\t', index = False) 

    tread_dict = pd.DataFrame(tread)
    tread_dict.to_csv("MYC_H1_whole_genome_binding_sites_tread.csv", sep = '\t', index = False) 

    return ctx_dict,cread_dict, tread_dict









