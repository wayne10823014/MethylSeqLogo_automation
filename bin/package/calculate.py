import pandas as pd
import numpy as np
import math
import os
from functools import reduce
from operator import itemgetter
from collections import defaultdict

dinucs= ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
nucleotides = ['A', 'C', 'G', 'T']

dir_path = os.path.dirname(os.path.realpath(__file__))
path = dir_path + "/../../Output1/"


def calc_ppm(seqdata, TF, species, celltype, region, pseudocount=1):
	"""
	Calculate  nucleotides ('A', 'C', 'G', 'T') for each position(one column) in the motifs\n

	#### Return

	- ppm (_dataFrame_)
	"""

	pd.options.display.float_format = '{:,.4f}'.format
	ppm = []
	for pos_index in seqdata:
		freqs = []
		
		for base in nucleotides:
			freqs.append(
				float(seqdata[pos_index].tolist().count(base)) + pseudocount)
		nucleotides_sum = sum(freqs)
		pps = map(lambda X: X/nucleotides_sum, freqs)
		ppm.append(pps)

	ppm = pd.DataFrame(ppm, columns=nucleotides)

	ppmname = TF + '_' + species + '_' + celltype + '_' + region + '_' + 'ppm.txt'
	ppm.to_csv(path + ppmname, sep='\t', index=False,
			   header=True, float_format='%.4f')

	# print("Position probability matrix: ")
	# print(ppm)
	print('\n')

	return ppm


def calc_methyllevel(sitelist, creaddata, treaddata, index, J_bg_methyl_prob):
	# J_bg_methyl_prob: p(mCG, mCHG, mCHH|bg)

	mlevel = 0.0
	if len(sitelist) > 0:
		JiC_ps = []
		for site in sitelist:
			sumofcread_C_p = creaddata[index][site]
			sumoftread_C_p = treaddata[index][site]
			if (sumofcread_C_p + sumoftread_C_p) >= 4.0:  # depth >= 4.0
				JiC_ps.append((sumofcread_C_p) /
							  (sumofcread_C_p + sumoftread_C_p))
			else:
				JiC_ps.append(J_bg_methyl_prob)
		mlevel += round((((np.sum(JiC_ps)) + J_bg_methyl_prob) /
					   (len(JiC_ps)+1)), 4)
	else:
		mlevel += round(J_bg_methyl_prob, 4)

	return mlevel


# -----------------------------------------------------------
def calc_methylprob(ctxdata, creaddata, treaddata, bg_mCG, bg_mCHG, bg_mCHH, plotlen):
	"""
	Calculate conditional methylation probabilities: CG, CHG, CHH
	"""
	JiCs, PiCs, Cmethyls, Gmethyls, PiCs_, PiGs_, C_ratio, G_ratio= [], [], [], [], [], [], [], []	

	methylcond_given_c = {
						  'CG': ['X', bg_mCG],
						  'CHG': ['Y', bg_mCHG], 
						  'CHH': ['Z', bg_mCHH]
	}
	
	# for 迴圈跑的是motif的每個位置(序列長度)
	for pos_of_tfbs in list(ctxdata.columns.values):
	   
		sub_JiC, sub_PiC, sub_PiC_, sub_PiG_ = [], [], [], []

		total_C= float(ctxdata[pos_of_tfbs].isin(['X', 'Y', 'Z']).sum())
		total_G= float(ctxdata[pos_of_tfbs].isin(['x', 'y', 'z']).sum())

		# for 迴圈跑的是(CG, CHG, CHH)
		for key in methylcond_given_c:

			methylcond = methylcond_given_c[key][0]
			# mCG, mCHG, mCHH probability from forward(+) strand on per pos with each TFBSs
			sites = ctxdata.index[ctxdata[pos_of_tfbs].isin([ methylcond ])].tolist()
			forward_mlevel = calc_methyllevel(sites, creaddata, treaddata, pos_of_tfbs, methylcond_given_c[key][1])

			# mCG, mCHG, mCHH probability from reverse(-) strand on per pos with each TFBSs
			sites = ctxdata.index[ctxdata[pos_of_tfbs].isin([methylcond.lower()])].tolist()
			reverse_mlevel = calc_methyllevel(sites, creaddata, treaddata, pos_of_tfbs, methylcond_given_c[key][1])
			sub_JiC.append((forward_mlevel, reverse_mlevel))

			# Probability of C contexts, denominator is the number of total sites(計算每個位置的CG, CHG, CHH占比)
			condition_forward = ctxdata[pos_of_tfbs].isin([methylcond]).sum()
			condition_reverse = ctxdata[pos_of_tfbs].isin([methylcond.lower()]).sum()
			
			forward_prob = float(condition_forward) / len(ctxdata.index)
			reverse_prob = float(condition_reverse) / len(ctxdata.index)
			sub_PiC.append((forward_prob, reverse_prob))
			
			#p(CG, CHG, CHH|C/G), denominator is the number of total C or G
			if total_C!= 0.0:
				sub_PiC_.append(float(condition_forward) / total_C )
			else:
				sub_PiC_.append(0.0)
			
			if total_G!= 0.0:
				sub_PiG_.append(float(condition_reverse) / total_G)
			else:
				sub_PiG_.append(0.0)


		JiCs.append(sub_JiC) #每個位置3種情況(CG, CHG, CHH)的甲基化程度

		PiCs.append(sub_PiC) #每個位置3種情況(CG, CHG, CHH)各佔的比例

		Cmethyls.append(tuple([sub_JiC[i][0] for i in range(3)]))  #每個位置正股3種情況(CG, CHG, CHH)的甲基化程度
		Gmethyls.append(tuple([sub_JiC[i][1] for i in range(3)]))   #每個位置反股3種情況(CG, CHG, CHH)的甲基化程度
		C_ratio.append(tuple([sub_PiC[i][0] for i in range(3)])) 
		G_ratio.append(tuple([sub_PiC[i][1] for i in range(3)]))
		

		PiCs_.append(tuple(sub_PiC_))#每個位置已知正股 C 鹼基下3種情況(CG, CHG, CHH)各佔的比例
		PiGs_.append(tuple(sub_PiG_))#每個位置已知反股 G 鹼基下3種情況(CG, CHG, CHH)各佔的比例

	#　JiCs :　[[(0.0663, 0.0976), (0.0247, 0.0041), (0.007, 0.0033)],[(0.1158, 0.07), (0.0062, 0.0051), (0.0038, 0.0051)]]
	#　output : [[0.0663, 0.0976, 0.0247, 0.0041, 0.007, 0.0033], [0.1158, 0.07, 0.0062, 0.0051, 0.0038, 0.0051]]
	Methyls= [list(sum(i, ())) for i in JiCs] 
	Methyls= pd.DataFrame(Methyls, columns= ['mCpG_p', 'mCpG_m', 'mCHG_p', 'mCHG_m', 'mCHH_p', 'mCHH_m'])

	Methyls.insert(loc= 0, column= 'position', value= list(range(1, plotlen + 1)))	

	# print ("Methylation probabilities: ")
	# print (Methyls)
	# print ('\n')

	Freqs= [list(sum(i, ())) for i in PiCs]
	Freqs= pd.DataFrame(Freqs, columns= ['CpG_p', 'CpG_m', 'CHG_p', 'CHG_m', 'CHH_p', 'CHH_m'])
	Freqs.insert(loc= 0, column= 'position', value= list(range(1, plotlen + 1)))

	# print ("Context probabilities: ")
	# print (Freqs)
	# print ('\n')

	# Freqs_ = [list(sum(i, ())) for i in PiCs]
	Freqs_C= pd.DataFrame(PiCs_, columns= ['CpG_p', 'CHG_p', 'CHH_p'])
	Freqs_C.insert(loc= 0, column= 'position', value= list(range(1, plotlen + 1)))
	Freqs_G = pd.DataFrame(PiGs_, columns= ['CpG_m', 'CHG_m', 'CHH_m'])
	Freqs_C_only = pd.concat([Freqs_C, Freqs_G], axis= 1)

	# print ("Context probabilities (C only): ")
	# print (Freqs_C_only)
	# print ('\n')

	return C_ratio, G_ratio, Cmethyls, Gmethyls, Freqs_C_only




def calc_methylation_entropy(C_ratio, G_ratio, Cmethyls, Gmethyls, J_bCG, J_bCHG, J_bCHH, logotype):

	bg_prob = [J_bCG, J_bCHG, J_bCHH]
	Cents= []
	if logotype == 'Kullback-Liebler':
		for i in range(len(Cmethyls)):
			Cents_ = []	
			for cond in range(3):
				bg = bg_prob[cond]	
				px_p ,px_m = Cmethyls[i][cond], Gmethyls[i][cond]
  
				entropy_p = 0.0 if np.isnan(px_p) else C_ratio[i][cond] * ((px_p * math.log( px_p / bg, 2)) + (1 - px_p) * math.log((1 - px_p)/(1 - bg), 2))
				Cents_.append(entropy_p)
				
				entropy_m =  0.0 if np.isnan(px_m) else G_ratio[i][cond]  * ((px_m  * math.log( px_m  / bg, 2)) + (1 - px_m ) * math.log((1 -px_m )/(1 - bg), 2))
				Cents_.append(entropy_m)
				
			Cents.append(tuple(Cents_))
	else:
		for i in range(len(Cmethyls)):
			Cents_ = []	
			for cond in range(3):
				bg = bg_prob[cond]	
			    
				px_p ,px_m = Cmethyls[i][cond], Gmethyls[i][cond]

				entropy_p = 0.0 if np.isnan(px_p) else C_ratio[i][cond] * (1 - (px_p * math.log( px_p , 2)) + (1 - px_p) * math.log((1 - px_p), 2))
				Cents_.append(entropy_p)

				entropy_m = 0.0 if np.isnan(px_m) else  G_ratio[i][cond] * (1 - (px_m * math.log( px_m , 2)) + (1 - px_m) * math.log((1 -px_m), 2))
				Cents_.append(entropy_m)
				
			Cents.append(tuple(Cents_))
	Cents= pd.DataFrame(Cents, columns= ['CpG_p', 'CpG_m', 'CHG_p', 'CHG_m', 'CHH_p', 'CHH_m'])
	Cents['Methylation']= Cents[Cents.columns].sum(axis= 1)
	return Cents



def calc_totalEntropy(fourletterppm, bgpps, Cents, logotype, plotlen, TF, species, celltype, region):
	"""
	Calculate column entropy
	"""

	# entropys from bases, mode including Shannon and KL
	pwmentropys = []
	if logotype == 'Kullback-Liebler':
		print("Kullback-Liebler distance: ")
		for i in list(range(len(fourletterppm))):
			pos_ent = 0
			pps = fourletterppm.iloc[i]
			for (pp1, pp2) in zip(pps, bgpps):
				ent = (pp1 * math.log(pp1/pp2, 2))
				pos_ent += ent
			pwmentropys.append(pos_ent)
	else:
		print("Shannon entropy: ")
		for i in list(range(len(fourletterppm))):
			pos_ent = 0
			pps = fourletterppm.iloc[i]
			for pp1 in pps:
				ent = -(pp1 * math.log(pp1, 2))
				pos_ent += ent
			pwmentropys.append(2 - pos_ent)

	# make table of entropys
	# Mentropys = [sum(i) for i in Cents]
	entropys = Cents.assign(Base=pd.Series(pwmentropys))
	# entropys = pd.DataFrame.from_records(zip(Mentropys, pwmentropys), columns = ['Methylation', 'Base'])
	entropys['Total'] = entropys.loc[:, entropys.columns.isin(
		['Methylation', 'Base'])].sum(axis=1)
	entropys.insert(loc=0, column='position', value=range(1, plotlen + 1))
	print(entropys)
	print("\n")

	if logotype == 'Kullback-Liebler':
		icname = TF + '_' + species + '_' + celltype + \
			'_' + region + '_' + 'relative_entropy.txt'
		entropys.to_csv(path + icname, sep='\t', index=False,
						header=True, float_format='%.2f')

	else:
		icname = TF + '_' + species + '_' + celltype + \
			'_' + region + '_' + 'Shannon_entropy.txt'
		entropys.to_csv(path + icname, sep='\t', index=False,
						header=True, float_format='%.2f')

	return entropys


def calc_perBaseEntropy(entropys, mfourletterppm, mode, TF, species, celltype, region):
	"""
	Calculate per base nucleotide height (4-letter)
	"""
    # entropys 是每個位置的總和，這裡算位置中acgt的資訊量
	four_base_heights = []

	for i in list(range(len(mfourletterppm))):
		base_ents = []
		for j in nucleotides:
			base_ent = None
			if mode == 'Methyl':
				base_ent = round(entropys['Base'][i] * mfourletterppm[j][i], 5)
			base_ents.append((j, base_ent))
		base_ents = sorted(base_ents, key=itemgetter(1, 0))
		four_base_heights.append(base_ents)

	for i, four_base_height in enumerate(four_base_heights):
		four_base_height = sorted(four_base_height, key=itemgetter(1, 0))

	filename = TF + '_' + species + '_' + celltype + \
		'_' + region + '_per_base_information_content.txt'
	try:
		with open(path + filename, "w") as outfile:
			outfile.write("position" + "\t" + "A" + "\t" +
						  "C" + "\t" + "G" + "\t" + "T" + "\n")
			for i, four_base_height in enumerate(four_base_heights):
				four_base_height = sorted(four_base_height, key=itemgetter(0))
				outfile.write(str(i + 1) + "\t" + str("%.2f" % four_base_height[0][1]) + "\t" + str("%.2f" % four_base_height[1][1]) + "\t" + str(
					"%.2f" % four_base_height[2][1]) + "\t" + str("%.2f" % four_base_height[3][1]) + "\n")
	except IOError as e:
		print("Unable to open file: " + filename)
		print(filename + "檔案無法開啟。")
		print(filename + "ファイルをオプン出来ませんでした。")

	# print (four_base_heights)

	return four_base_heights


def to4basedippm(seqdata, plotlen):
	"""
	Build 4-letter dinucleotide probability matrix
	"""
	probs = []
	for i in list(range(plotlen - 1)):
		prob = []
		row = zip(seqdata.iloc[:, i], seqdata.iloc[:, i+1])
		dinucseq = []
		for j, k in enumerate(row):
			string = ''.join(k)
			dinucseq.append(string)
		for dinuc in dinucs:
			pro = float(dinucseq.count(dinuc)+1) / (len(seqdata.index)+16)
			prob.append(pro)  # length = 16
		probs.append(prob)

	dippm = pd.DataFrame(probs, columns=dinucs)
	dippm = dippm.T.astype('float64')

	# print("Dinucleotide probability matrix: ")
	# print(dippm)

	return dippm


def twomerBg(bgpps, dippm, plotlen):

	dinucs = []
	for i in nucleotides:
		for j in nucleotides:
			dinucs.append(i + j)

	dimerenrichments = {}
	for i in nucleotides:
		for j in nucleotides:
			dinuc = i+j
			entscore = math.log(bgpps.loc[dinuc], 2) - \
				math.log(bgpps.loc[i]*bgpps.loc[j], 2)
			dimerenrichments[dinuc] = entscore

	dientropis =  defaultdict(list) 
	for i in list(range(plotlen - 1)):
		for dinuc in dinucs:
			d1 = dinuc[0]
			d2 = dinuc[1]
			dientropis[dinuc].append(
				dippm.loc[dinuc][i] * (math.log(bgpps.loc[d1]*bgpps.loc[d2], 2) - math.log(bgpps.loc[dinuc], 2)))
	print("\n")
	# print("Dimer enrichment scores: " + "\n")
	series = pd.Series(dimerenrichments)
	# print(series)
	# print("\n")

	print("Depleted dimer binding scores: " + "\n")
	series2 = pd.DataFrame.from_dict(dientropis, orient='columns')
	series2['Total'] = series2.sum(axis=1)
	series2.insert(loc=0, column='position', value=[
				   str(i) + '-' + str(i+1) for i in list(range(1, plotlen))])
	print(series2)
	print("\n")

	bg_dientropis = {}
	for dinuc in dinucs:
		d1 = dinuc[0]
		d2 = dinuc[1]
		bg_dientropis[dinuc] = (
			(math.log(bgpps.loc[d1]*bgpps.loc[d2], 2) - math.log(bgpps.loc[dinuc], 2)))

	# print (bg_dientropis)

	return series2['CG'].tolist(), bg_dientropis['CG'], bg_dientropis['CC']
