import argparse
import time
from multiprocessing import Pool
import multiprocessing as mp
from functools import reduce
import numpy as np
import os
import itertools
import pandas as pd

'計算humun的 cpg, chg, chh, A, C, G, T,AA, AC..TT的比例'
'有用Multiproccess加速從 20 -> 10 min'
'只計算whole genome'

def get_parser():
    parser = argparse.ArgumentParser('background_probability')
    parser.add_argument('-c', '--celltype', choices= ['H1-hESC', 'HepG2', 'K562', 'GM12878', 'HeLa-S3', 'naive-hESC', 'primed-hESC', 'mESC', 'mMPGC_E13p5', 'mMPGC_E16p5', 'mFPGC_E13p5', 'mFPGC_E16p5', 'Leaf', 'Inflorescence', 'Root', 'shoot'], default= 'HepG2', help= 'What type of cells is the WGBS_file from? ', metavar= 'cell types')    
    parser.add_argument('-s', '--speciesname', help= 'species of the input data for MethylSeqLogo', metavar= 'species', default= 'human')
    parser.add_argument('-fa', '--fasta', help= 'fasta of species', metavar= 'species', default= '/genome/human.fa')

    
    wgbs = parser.add_argument_group('arg with methylation calculation')
    wgbs.add_argument('-CG','--mcg', nargs='+', type=str, help= 'WGBS at CpG')
    wgbs.add_argument('-CHG','--mchg',  nargs='+', type=str, help= 'WGBS at CHG')
    wgbs.add_argument('-CHH','--mchh',  nargs='+', type=str, help= 'WGBS at CHH')
    
    return parser

def  read_big_file(big_fileA, big_fileB):

    while True:
        data_chunk_from_fileA = big_fileA.readline()
        data_chunk_from_fileB = big_fileB.readline()
        if data_chunk_from_fileA == '' or data_chunk_from_fileB == '':
            break
        
        yield data_chunk_from_fileA.replace('\n',''),data_chunk_from_fileB.replace('\n','')


def calc_mlevel(fileA, fileB):

    total, count, mismatch = 0, 0, 0
    length = 900000000
    # s_time = time.time()
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
    
    
    # end =  time.time()
    # print('花費時間：', (end - s_time)/60,'min') 
    # print('m_prob:',total /count)
    print('mismatch',mismatch)
    return '{:.4f}'.format(total /count)


        

def to_dict(lst, val):
    res_dict = {lst[i]: val[i] for i in range(0, len(lst), 1)}
    return res_dict 

def calc(seq):
    
    if seq.startswith('>') or (seq.startswith('N') and len(set(seq))==2):
            return 0

    
    # monome_nucleobase_list = ['A', 'C', 'G', 'T', ]
    # dimer_nucleobase_list = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    
    # 不cal to_dict會加速嗎?
    # {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0, 'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0, 'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0, 'TA': 0, 'TC': 0, 'TG': 0, 'TT': 0}
    # {'A': 0, 'C': 0, 'G': 0, 'T': 0}
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


def read_file(file):    
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
    df = pd.DataFrame(cpg_chg_chh, columns=['whole_genome'],index=[ 'CpG', 'CHG', 'CHH', 'A', 'C', 'G', 'T', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])
    # print(df)
    path = dir_path + '/../Background_probability/whole_genome/'+ speciesname +'_whole_genome_probability.txt'
    df.to_csv(path, sep='\t')



def  main():
    
    global cpu_count, dir_path, speciesname, wgbs_file, celltype
    cpu_count = mp.cpu_count()
    dir_path = os.path.dirname(os.path.realpath(__file__))

    parser = get_parser()
    args = parser.parse_args()
  
    species_fa= dir_path + args.fasta
    celltype = args.celltype
    speciesname  =  args.speciesname
    wgbs_file = [args.mcg, args.mchg, args.mchh]

    print('\n')
    print('File:',species_fa)
    print('The computational will take approximately 10 minutes.\n')
    stime = time.time()  
    with open(species_fa, 'r') as file:
        file = file.readlines()
        etime = time.time()            
        print('read file finish cost ', (etime-stime)/60,' min \n' )
        read_file(file)
        print('\n')

    if wgbs_file[0] != None:
        print('start calculating mlevel...')


        start = time.time()
        CHG_result = calc_mlevel(wgbs_file[1][0], wgbs_file[1][1])
        end = time.time()
        print('CHG file finish cost ', (end-start)/60,' min \n' ) 

        start = time.time()
        CHH_result = calc_mlevel(wgbs_file[2][0], wgbs_file[2][1])
        end = time.time()
        print('CHH file finish cost ', (end-start)/60,' min \n' )   

        start = time.time()
        CG_result = calc_mlevel(wgbs_file[0][0], wgbs_file[0][1])
        end = time.time()
        print('CG file finish cost ', (end-start)/60,' min \n' )   



        df = pd.DataFrame([CG_result, CHG_result, CHH_result], columns=['whole_genome'],index=['mCG','mCHG', 'mCHH' ])
        path = dir_path + '/../Background_probability/whole_genome/'+ speciesname + '_'+ celltype +'_whole_genome_methyl_probability.txt'
        df.to_csv(path, sep='\t', float_format = '%.4f')



if __name__ == '__main__':
    main()