import argparse

def get_parser():
    parser = argparse.ArgumentParser('MethylSeqLogo')
    inputfile = parser.add_argument_group('required file input')
    inputfile.add_argument("-J", "--jaspar", required= True)
    inputfile.add_argument("-R", "--remap", required= True)
    inputfile.add_argument('-M','--methylationinfo', nargs='+', type=str, help= 'WGBS at CpG')
    inputfile.add_argument('-CG','--mcg', required= True, nargs='+', type=str, help= 'WGBS at CpG')
    inputfile.add_argument('-CHG','--mchg', required= True, nargs='+', type=str, help= 'WGBS at CHG')
    inputfile.add_argument('-CHH','--mchh', required= True, nargs='+', type=str, help= 'WGBS at CHH')
    

    tfbs = parser.add_argument_group('arg with TFBSs calculate')
    tfbs.add_argument('-TF', '--transcriptionfactor', help= 'transcription factor', metavar= 'TF', required= False, default= 'ZBTB33')
    tfbs.add_argument('-S', '--species', choices= ['human', 'mouse', 'arabidopsis', 'maize'], default= 'human', help= 'species of the input data', metavar= 'species')
    tfbs.add_argument('-C', '--celltypes', choices= ['H1-hESC', 'HepG2', 'K562', 'GM12878', 'HeLa-S3', 'naive-hESC', 'primed-hESC', 'mESC', 'mMPGC_E13p5', 'mMPGC_E16p5', 'mFPGC_E13p5', 'mFPGC_E16p5', 'Leaf', 'Inflorescence', 'Root', 'shoot'], default= 'HepG2', help= 'type of cell', metavar= 'cell types')    
    tfbs.add_argument('-T', '--threshold', help= 'threshold to display methylation probability', default= 0.4, ) #metavar= 'display threshold'
    tfbs.add_argument('-ROI', '--regionofinterest', default= 'whole_genome', help= 'genomic region of interest', ) #metavar= 'genomic region'
    # tfbs.add_argument('--bg')
    
    draw = parser.add_argument_group('arg with MethylSeqLogo display')
    draw.add_argument("-bot", "--beginningoftfbs", type= int, default= 1, required = False,
                      help= 'start position of TFBS to plot logo')
    draw.add_argument('-plen', '--plotlen', help= 'length of TFBS to plot logo', type= int, required = False)
    draw.add_argument('-m', '--mode', choices= ['Methyl', 'Regular'], default= 'Methyl', help= 'include methylation or not', ) #metavar= 'mode of plotting'
    draw.add_argument('-l', '--logotype', choices= ['Shannon', 'Kullback-Liebler', 'riverlake', 'all'], default= 'Kullback-Liebler', help= 'logo type', ) # metavar= 'type of logo'
    

    return parser
