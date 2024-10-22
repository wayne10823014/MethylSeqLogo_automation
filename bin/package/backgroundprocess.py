import os
import pandas as pd
import pybedtools 
import time
from multiprocessing import Pool
import pandas as pd
from functools import reduce
import multiprocessing as mp
from package.datapreprocess import  modify_bed

import numpy as np
import itertools

dir_path = os.path.dirname(os.path.realpath(__file__))
path = dir_path + "/../../Output1/"
global cpu_count
cpu_count = mp.cpu_count()

def isfile(path,flag=1):
    pd.set_option('display.float_format', lambda x: '%.4f' % x)
    if (flag):
        if os.path.isfile(path):
            pass
        else:
            print (path + " does not exist. Job cancelled.")
            print (path + "檔案不存在，工作終止。")
            print (path + "はありませんでした。ジョブがキャンセルされました。")
    if os.path.isfile(path):
        data = pd.read_table(path, sep= "\t", header= 0, index_col= 0)
        return data
    else:
        print (path + " does not exist. Job cancelled.")
        print (path + "檔案不存在，工作終止。")
        print (path + "はありませんでした。ジョブがキャンセルされました。")
        
    
    

def read_bgprob_table(species, celltype, region,TF=None):
    """
    Read background probability table \n
    #### bg_table content
    ['A', 'C', 'G', 'T', 'CpG', 'CHG', 'CHH', 'mCG', 'mCHG', 'mCHH', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    #### return 
    - bgpps : _DataFrame_  -> 
                    A :0.2955
                    AA:0.0979
    - mCG : _int_ 
    - mCHG : _int_
    - mCHH : _int_
    """
    if region == 'whole_genome':
        species_path = dir_path + "/../../Background_probability/"+ region + "/" + species + '_' + region + '_probability.txt'
        celltype_path = dir_path + "/../../Background_probability/"+ region + "/" + species + '_' + celltype + '_' + region  +'_methyl_probability.txt'
    else:
        species_path = dir_path + "/../../Background_probability/neighbor/"  + species + '_' + TF + '_' + celltype + '_' + region + '_probability.txt'
        # species_path = dir_path + "/../../Background_probability/"+ region + "/"  + species + '_' + celltype + '_' + region + '_probability.txt'
        celltype_path = dir_path + "/../../Background_probability/neighbor/" + species + '_' + TF + '_' + celltype + '_' + region  +'_methyl_probability.txt'
    
    # path = 'MethylSeqLogo_automation/Background_probability/Background_probability/human_HepG2_probability.txt' 
    species_path = isfile(species_path)
    celltype_path = isfile(celltype_path)
   
    region = 'neighbor' if region.isdigit()  else region
    
    data= pd.DataFrame(species_path)
    probmatrix= data[region].astype('float64')

    bgpps= probmatrix[['A', 'C', 'G', 'T', 'CpG', 'CHG', 'CHH', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']]
    
    # print ("background probabilities: ")
    # print (bgpps)
    # print ("\n")

    data= pd.DataFrame(celltype_path)
    probmatrix= data[region].astype('float64')

    bg_mCG= probmatrix['mCG'].astype('float64') 
    bg_mCHG= probmatrix['mCHG'].astype('float64') 
    bg_mCHH= probmatrix['mCHH'].astype('float64') 

    print ("Background methylation probabilities: ")
    print (probmatrix[['mCG', 'mCHG', 'mCHH']])
    print ("\n")

    return bgpps, bg_mCG, bg_mCHG, bg_mCHH


def calc(seq):
    
    if seq.startswith('>') or (seq.startswith('N') and len(set(seq))==2):
            return 0

    dimer_nucleobase_dict = {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0, 'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0, 'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0, 'TA': 0, 'TC': 0, 'TG': 0, 'TT': 0}
    monome_nucleobase_dict =  {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    cpg, chg, chh= 0, 0, 0

    
    for base in  range(0, len(seq), 1):
        #calculate ['A', 'C', 'G', 'T', 'AA', 'AC' ... 'TG', 'TT']
        if seq[base] in ["c", "C", "g", "G", "a", "A", "T", "t"]:
            monome_nucleobase_dict[seq[base].upper()] += 1
            if base  < len(seq)-1 and seq[base+1] in ["c", "C", "g", "G", "a", "A", "T", "t"]:
                dimer = seq[base] + seq[base+1]
                dimer_nucleobase_dict[dimer.upper()] += 1
                
        # calculate condition when given C/G
        if  base >= 2 and base <= len(seq)-3:
            
            if seq[base] in ["c", "C"]:
            
                #CG condition | C base +
                if seq[base+1] in ["g","G"] :
                    cpg += 1
                    
                #CHG condition | C base +
                elif seq[base+2] in ["g","G"]:
                    chg += 1
                    
                #CHH condition | C base +   
                else:
                    chh += 1
                    
            elif seq[base] in ["G", "g"]:
                
                #CG condition | G base +
                if seq[base-1] in ["c", "C"] :
                    cpg += 1
                    
                #CHG condition | G base +
                elif seq[base-2] in ["c", "C"]:
                    chg += 1
                    
                #CHH condition | G base +
                else:
                    chh += 1 
           


    return [cpg, chg, chh, monome_nucleobase_dict, dimer_nucleobase_dict] 
    # return [cpg, chg, chh, monome_nucleobase_dict, dimer_nucleobase_dict]                    


def read_file(file, species, celltype, region , dir,TF):    
    print('~~read_file~~') 
    filestime = time.time()
    cpg, chg, chh =[0, 0, 0] 

    monome, dimer = np.zeros(4), np.zeros(16)     

    # stime = time.time()
    # with open(file, 'r') as file:
    #     file = file.readlines()
    #     etime = time.time()            
    #     print('read file finish cost ', (etime-stime)/60,' min \n' )

    with Pool(cpu_count) as pool:
            # res = pool.map(calc,[file[seq] for seq in range(len(file))])
            res = pool.map(calc,file)
                
    for i in res:
        if isinstance(i, int):
            continue

        cpg += i[0]
        chg += i[1]
        chh += i[2]

        monome += np.array(list(i[3].values()))
        dimer += np.array(list(i[4].values()))

    sum = cpg + chh + chg
    
    cpg_chg_chh = ['{:.4f}'.format(cpg/sum), '{:.4f}'.format(chg/sum), '{:.4f}'.format(chh/sum)]
    # print('CpG', cpg,round(cpg /sum,4))
    # print('CHG',chg,round(chg /sum,4))
    # print('CHH',chh,round(chh /sum,4))  

    # print('monome',monome.tolist())
    # print('dimer',dimer.tolist())
    sum =  reduce(lambda x,y : x + y, monome.tolist())
    monome_prob = list(map(lambda x : round( x / sum, 4), monome.tolist()))
    # print(monome_prob)
    for monome in  monome_prob:
       
        cpg_chg_chh.append('{:.4f}'.format(monome))


    dimer_ = dimer.tolist()
    m = []
    index = 0
    for i in range(0, len(dimer_ ), 4):
        temp = dimer_ [i : i + 4]
        sum = reduce(lambda x, y : x + y, temp)
        sum = monome_prob[index] / sum
        for i in range(4):
            m.append(sum)
        index += 1
    m = np.array(m)    
    dimer_prob = np.multiply(dimer, m)
    etime = time.time()  
    # print(dimer_prob)

    for dimer in  dimer_prob:
        cpg_chg_chh.append('{:.4f}'.format(dimer))

    print('calc finish..', (etime-filestime)/ 60,' min' ) 
    df = pd.DataFrame(cpg_chg_chh, columns=[dir],index=[ 'CpG', 'CHG', 'CHH', 'A', 'C', 'G', 'T', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])
    # print(df)
    path = dir_path + '/../../Background_probability/'+ dir + '/'+ species + '_' + TF + '_' + celltype +'_'+region +'_probability.txt'
    df.to_csv(path, sep='\t')
    

def  read_big_file(big_fileA, big_fileB):

    while True:
        data_chunk_from_fileA = big_fileA.readline()
        data_chunk_from_fileB = big_fileB.readline()
        if data_chunk_from_fileA == '' or data_chunk_from_fileB == '':
            break
        
        yield data_chunk_from_fileA.replace('\n',''),data_chunk_from_fileB.replace('\n','')


def calc_mlevel(fileA, fileB):
    print('~~calc_mlevel~~')
    total, count, mismatch = 0, 0, 0
    length = 900000000
    s_time = time.time()
    with Pool(cpu_count-2) as pool:
        with open(fileA, 'r') as big_fileA, open(fileB, 'r') as big_fileB:
            for wgbs_from_fileA, wgbs_from_fileB in itertools.islice(read_big_file(big_fileA, big_fileB), length):
                wgbs_from_fileA, wgbs_from_fileB = wgbs_from_fileA.split('\t'), wgbs_from_fileB.split('\t')
                try:
                    if int(wgbs_from_fileA[-2]) < 4 and int(wgbs_from_fileB[-2]) < 4 :
                        continue
                    elif (wgbs_from_fileA[-2]) == '0':
                        total += int(wgbs_from_fileB[-1])/100 
                        count += 1
                    elif int(wgbs_from_fileB[-2]) == '0':   
                        total += int(wgbs_from_fileA[-1])/100
                        count += 1
                    else:
                        read_A = 4 if int(wgbs_from_fileA[-2]) >= 4 else int(wgbs_from_fileA[-2])
                        read_B = 4 if int(wgbs_from_fileB[-2]) >= 4 else int(wgbs_from_fileB[-2])
                        total += ((int(wgbs_from_fileA[-1])*read_A+ int(wgbs_from_fileB[-1])*read_B)/100) / (read_A+read_B)
                        count += 1
                except:
                    pass
        
            pool.close()
            pool.join()
    
    
    end =  time.time()
    print('花費時間：', (end - s_time)/60,'min') 
    print('m_prob:',total /count)
    print('mismatch',mismatch,'\n')
    return '{:.4f}'.format(total /count)



def culc_mlevel(tfbs_bed, wgbs_file, speciesname, celltype, region, dir, TF):
    print('~~culc_mlevel neighbor~~')
    TFBSs = tfbs_bed
    name = ['CG', 'CHG', 'CHH'] 
    total = []
    count = 0
    for file in wgbs_file:
        start = time.time()
        wgbs = pybedtools.BedTool(file[0])
        wgbs_tfbs_A = wgbs.intersect(TFBSs, wa = True)
        wgbs_tfbs_A.sort()
        wgbs_tfbs_A = pd.DataFrame(wgbs_tfbs_A)
        path_name_A = dir_path +  "/../temp/" + '/' + celltype + '_WGBS_'  + name[count] + "_" + dir + "_region_1.bed"
        wgbs_tfbs_A.to_csv(path_name_A , sep = '\t', index = False, header = False) 
        end = time.time()
        print(name[count] )
        print('wgbs_tfbs_A intersect finish,cost',(end-start)//60,'min')
    
        start = time.time()   
        wgbs = pybedtools.BedTool(file[1])  
        wgbs_tfbs_B = wgbs.intersect(TFBSs, wa = True) 
        wgbs_tfbs_B.sort()
        wgbs_tfbs_B = pd.DataFrame(wgbs_tfbs_B)
        path_name_B = dir_path +  "/../temp/" + '/' + celltype + '_WGBS_'  + name[count] + "_" + dir+"_region_2.bed"
        wgbs_tfbs_B.to_csv(path_name_B, sep = '\t', index = False, header = False) 
        end = time.time()
        print('wgbs_tfbs_B intersect finish,cost',(end-start)//60,'min')
    

        result = calc_mlevel(path_name_A, path_name_B)
        total.append(result)
        count += 1
        


       



    df = pd.DataFrame(total, columns=[dir],index=['mCG','mCHG', 'mCHH' ])
    print('～～存檔中~~')
    path = dir_path + '/../../Background_probability/' + dir + '/'+ speciesname + '_' + TF + '_'+ celltype +'_' + region +'_methyl_probability.txt'
    df.to_csv(path, sep='\t', float_format = '%.4f')    
    print('～～存好了~~')



def promoter(tfbs_bed, species, wgbs_file, celltype, region):
    tfbs_bed = modify_bed(tfbs_bed,50)
    fasta = pybedtools.BedTool(dir_path + '/../genome/'+species+'.fa')
    tfbss_fa  = tfbs_bed[1].sequence(fi = fasta,s=True)
    tfbss_fasta = open(tfbss_fa.seqfn).readlines()
    read_file(tfbss_fasta, species, celltype, region, 'promoter')
    culc_mlevel(tfbs_bed[1], wgbs_file, species, celltype, 'promoter')

    return read_bgprob_table(species, celltype, 'promoter')


def neighbor(tfbs_bed, species, wgbs_file, celltype, region, TF):
    
    fasta = pybedtools.BedTool(dir_path + '/../genome/'+species+'.fa')
    tfbss_fa  = tfbs_bed.sequence(fi = fasta,s=True)
    tfbss_fasta = open(tfbss_fa.seqfn).readlines()
    read_file(tfbss_fasta, species, celltype, region, 'neighbor',TF)
    culc_mlevel(tfbs_bed, wgbs_file, species, celltype, region, 'neighbor',TF)

    return read_bgprob_table(species, celltype, region,TF)

     
