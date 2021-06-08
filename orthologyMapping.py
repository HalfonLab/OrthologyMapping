#!/usr/bin/env python
################################################################
# Orthology mapping- 
#
# (c) Hasiba Asma July 2020
# Commented: June 2021
# use after post processing of SCRMshaw to get ortholog/paralog wrt DMEL
# input files can be downloaded from https://ezmeta.unige.ch/i5k/
#
################################################################


import os
import sys
import argparse
import pprint
import re

#This dictionary is created to save id mapping file of Spec_X. The format of this dictionary is something like this: [OFAS00001] = 'OFAS2:0001'
def idMap_dict(nameOfDict,file):
	with open(file,'r') as fi:
		rows = (line.split(' ') for line in fi)
		nameOfDict = {row[1].strip('\n'):row[0] for row in rows}
	return(nameOfDict)

#This dictionary is created to save best reciprocal hits (DMEL_SpecX or SpecX_DMEL format). The format of this dictionary is something like this: [7227:00001] = 'OFAS2:0001' 
def brh_dict(nameOfBrhDict,fileBrh,flipping):
	with open(fileBrh,'r') as fb:
		rowsB=(line.split('\t') for line in fb)
		if (flipping == 'true'):
			#print('true')
			nameOfBrhDict= {rowb[0]:rowb[1] for rowb in rowsB}
		else:
			#print('False')
			nameOfBrhDict= {rowb[1]:rowb[0] for rowb in rowsB}
	return(nameOfBrhDict)

#This dictionary is created to save id mapping file of DMEL. The format of this dictionary is something like this: [7227:0001] = 'FBgn0264125'
def idMap2_dict(nameOfDict2,file2):
	with open(file2,'r') as fi:
		rows = (line.split(' ') for line in fi)
		nameOfDict2 = {row[0]:row[1].strip() for row in rows}
	return(nameOfDict2)

#This function takes in the list of paralogs (via *.97.pc.txt file) and returns the LIST of list of paralogs e.g. ['OFAS2:004b90', 'OFAS2:004b2b']
def paralogs(nameoflist,filepc):
	listOflist=[]
	with open(filepc,'r') as fp:
		rows=(re.split(',|\t',line.strip()) for line in fp)
		for row in rows: listOflist.append(row)
		#for row in rows:
			#listOflist.append(row)
	return(listOflist)


def find(value,matrix):
	for list in matrix:
		if value in list:
			return [matrix.index(list),list.index(value)]
	return -1
	
def conv_dict(nameOfCDict,fileC):
	with open(fileC,'r') as fi:
		rows = (line.split('\t') for line in fi)
		nameOfCDict = {row[1].strip():row[0] for row in rows}
	return(nameOfCDict)
	

def geneSymbol_dict(nameofSDict,fileS):
	with open(fileS,'r') as fi:
		rows = (line.split('\t') for line in fi if not (line.startswith('#') or line.strip()==''))
		nameofSDict = {row[0]:row[2] for row in rows}
	return(nameofSDict)

def orthologs_dict(nameOfODict_orth,nameOfDict_para,fileO):
	nameOfODict_orth={}
	nameOfDict_para={}
	
	with open(fileO,'r') as fi:
		rows = (line.split('\t') for line in fi)
		#nameOfCDict = {row[0]:row[1]+row[2] for row in rows}
		for row in rows:
			if row[2]!='NULL':
				nameOfDict_para[row[0]]=row[2]
			if row[4]!='NULL\n':
				nameOfODict_orth[row[0]]=row[4].strip()
			

	return(nameOfODict_orth,nameOfDict_para)
#-----------------------------MAIN FUNCTION-------------------
def main():
	#global d1
	parser=argparse.ArgumentParser()
	parser.add_argument('-np1','--namesp1',help='name of species 1',required=True)
	parser.add_argument('-sp1id','--sp1idmap',help='species 1 id map text file',required=True)
	parser.add_argument('-sp1pc','--sp1pc97',help='species 1 paralogs 97 pc file',required=True)
	parser.add_argument('-brh','--brhsp1sp2',help='best reciprocal file sp1 to sp2',required=True)
	parser.add_argument('-sp2id','--sp2idmap',help='species 2 id map text file',required=True)
	parser.add_argument('-sp2pc','--sp2pc97',help='species 2 paralogs 97 pc file',required=True)
	parser.add_argument('-geneSet','--geneSet',help='gene set from gff same as gene set used for orthology mapping',required=True)
	parser.add_argument('-symb','--symbolGene',help='Gene Symbol',required=True)
	parser.add_argument('-conv','--conversion',help='if conversion req or not',default=False)
	parser.add_argument('-setConv','--setConvGene',help='conversion file')
	parser.add_argument('-flip','--flipped',help='if in brh DMEL:specie1 set it to False',default='True')
	parser.add_argument('-sep','--separator',help='if separator is not colon, provide its value',default=':')
	parser.add_argument('-so','--scrmshawOutput',help='SO')
	args = parser.parse_args()
	#absolute path
	namesp1=args.namesp1
	sp1idmap=os.path.abspath(args.sp1idmap)
	sp1pc97=os.path.abspath(args.sp1pc97)
	brhsp1sp2=os.path.abspath(args.brhsp1sp2)
	sp2idmap=os.path.abspath(args.sp2idmap)
	sp2pc97=os.path.abspath(args.sp2pc97)
	geneSet=os.path.abspath(args.geneSet)
	symbolGene=os.path.abspath(args.symbolGene)
	conversion=args.conversion
	flipped=(args.flipped).lower()
	separator=args.separator
	scrmshawOutput=(args.scrmshawOutput)
	scrmshawOutputPath=os.path.abspath(scrmshawOutput)

	# Depending on naming convention of genes used for SCRMshaw and orthoDB, you might need to convert them to same convention first
	#and if there is need of conversion (or conversion=true), user need to provide the respective file needed for conversion and that
	#file is used to create a dictionary to save them
	if (conversion=="TRUE" or conversion=="T" or conversion=="true"):
		convGeneSet=os.path.abspath(args.setConvGene)
		dict_conv=conv_dict('conv',convGeneSet)

	#creating a dictionary of SpecX using idmap file (format: )
	dict_sp1id=idMap_dict('sp1id',sp1idmap)
	#pprint.pprint(dict_sp1id)
	#print('Flipped',flipped)
	#creating a dictionary to store best reciprocal hits (format: )
	dict_brh=brh_dict('brh',brhsp1sp2,flipped)
	# creating a dictionary of DMEL using idmap file (format: )
	dict_sp2id=idMap2_dict('sp2id',sp2idmap)
	# creating a dictionary to save symbol of genes (format: )
	dict_symb=geneSymbol_dict('symbol',symbolGene)
	
	#pprint.pprint(dict_symb)
	#creating list of paralog lists for specX and DMEL
	pc1=paralogs('pc1',sp1pc97)
	pc2=paralogs('pc2',sp2pc97)

	#creating two intermediate files that will store information regarding orthologs and paralogs
	f=open(namesp1+'_temp.txt','w')
	f.write('GeneName_TC'+'\t'+'GeneName_OGid'+'\t'+'NoParalogOrListOfParalogs'+'\t'+'ifParalogHasOrtholog'+'\t'+'Ortholog'+'\n')
	f1=open(namesp1+'_final.txt','w')
	f1.write('GeneName'+'\t'+'ParalogsThatHaveOrthologs'+'\t'+'GeneSymbolParalogs'+'\t'+'Orthologs'+'\t'+'GeneSymbolOrthologs'+'\n')

	#opening geneSet file that has 1 gene per line (all of genes extracted from its gff) and finding out if there are any orthologs or/and paralogs present wrt DMEL, if present write into two temp files created before.. to use later on to edit SCRMshaw prediction file to add respective ortho/para data in it
	with open(geneSet,'r') as gS, open(namesp1+'_temp.txt','a') as g1,open(namesp1+'_final.txt','a') as f1:
		for line in gS:
			geneName=line.rstrip()
			g1.write(geneName)
			#f1.write(geneName+'\t')
			# if there is a step of conversion of gene Naming involved, then making sure to use the right one using the dictioanary created previously
			if (conversion=="TRUE" or conversion=="T" or conversion=="true"):
				#print('require c')
				if geneName.split(separator)[1] in dict_conv:
					f1.write(dict_conv[geneName.split(separator)[1]]+'\t')
				else:
					f1.write(geneName+'\t')
			# if not, then go ahead with line one of gene name
			else:
				f1.write(geneName+'\t')

			#if gene is present in the idMap dictionary of Spec X
			if geneName in dict_sp1id:
			   #print('present')
			   	# write the id of this gene in temp file
				g1.write('\t'+dict_sp1id[geneName])
			   #using the id of geneName check to see if its present in paralog file (using pc list created before)
				if any(dict_sp1id[geneName] in sublist for sublist in pc1):
					#g1.write('\t'+'Has paralog')
					paralogListIndex=find(dict_sp1id[geneName],pc1)
					
					#listing all paralogs
					g1.write('\t'+str(pc1[paralogListIndex[0]]))
					#g1.write('\t')
					#print(str(pc1[paralogListIndex[0]]))
					paraOrtho=0
					for i in range(0,len(pc1[paralogListIndex[0]])):
						
						#print(i)
						#print(pc1[paralogListIndex[0]][i])
						if i==0:
							g1.write('\t')
						#g1.write('\t')
						#check if paralogs have otholog / or present in brh?
						if pc1[paralogListIndex[0]][i] in dict_brh:
							paraOrtho+=1
							#name of paralog that has ortholog
							#g1.write(pc1[paralogListIndex[0]][i])
							if dict_brh[pc1[paralogListIndex[0]][i]] in dict_sp2id:
								#print("This paralog has ortholg")
								
								#final
								#g1.write(dict_sp2id[dict_brh[pc1[paralogListIndex[0]][i]]]+',')
								g1.write(pc1[paralogListIndex[0]][i]+'--'+dict_brh[pc1[paralogListIndex[0]][i]]+'__'+dict_sp2id[dict_brh[pc1[paralogListIndex[0]][i]]])
								f1.write(dict_sp2id[dict_brh[pc1[paralogListIndex[0]][i]]]+',')
						#else:
							# has paralogs but none of them has ortholog
						#	g1.write('paralog No ortholog')
							
								#check its symbol
								if dict_sp2id[dict_brh[pc1[paralogListIndex[0]][i]]] in dict_symb:
									g1.write('/'+dict_symb[dict_sp2id[dict_brh[pc1[paralogListIndex[0]][i]]]]+',')
									f1.write('\t'+dict_symb[dict_sp2id[dict_brh[pc1[paralogListIndex[0]][i]]]]+',')
								else:
									g1.write('/'+'NULL'+',')
									f1.write('\t'+'NULL'+',')
							
							#else:		
					if (paraOrtho==0):
						g1.write('No ortho of para')
						f1.write('NULL'+'\t'+'NULL')
					#g1.write('\t')
					#g1.write('\t'+str(pc1[find(dict_sp1id[geneName],pc1)[0]]))
					
					
					
				else:
					g1.write('\t'+'NO paralog'+'\t'+'No PARA No ortho of para')
					f1.write('NULL'+'\t'+'NULL')
				
				#check if it has orthologs----
				if dict_sp1id[geneName] in dict_brh:
					#g1.write('Has ortholog')
					if dict_brh[dict_sp1id[geneName]] in dict_sp2id: 
						#g1.write(geneName+'\t'+dict_sp1id[geneName]+'\t'+'Has parlogs:'+str('sublist')+'\t'+dict_brh[dict_sp1id[geneName]]+'\t'+dict_sp2id[dict_brh[dict_sp1id[geneName]]]+'\n')
						g1.write('\t'+dict_brh[dict_sp1id[geneName]]+'__'+dict_sp2id[dict_brh[dict_sp1id[geneName]]])				
						f1.write('\t'+dict_sp2id[dict_brh[dict_sp1id[geneName]]])
						
						# Gene symbol:
						if dict_sp2id[dict_brh[dict_sp1id[geneName]]] in dict_symb:
							g1.write('/'+dict_symb[dict_sp2id[dict_brh[dict_sp1id[geneName]]]])
							f1.write('\t'+dict_symb[dict_sp2id[dict_brh[dict_sp1id[geneName]]]])
						else:
							g1.write('/'+'NULL')
							f1.write('\t'+'NULL')
					
					else:
						print('p in o')
					
				else:
					g1.write('\t'+'NO ortholog')
					f1.write('\t'+'NULL'+'\t'+'NULL')
				
			else:
				g1.write('\t'+'NULL'+'\t'+'NULL'+'\t'+'NULL'+'\t'+'NULL'+'\t'+'NULL')
				f1.write('NULL'+'\t'+'NULL'+'\t'+'NULL'+'\t'+'NULL')
			   #for sublist in pc1:
				   #for item in sublist:
			g1.write('\n')
			f1.write('\n')
	f.close()
	f1.close
	
	fof=os.path.abspath(namesp1+'_final.txt')
	dict_orthologs,dict_paralogs=orthologs_dict('orthologs','paralogs',fof)
	print(len(dict_orthologs))
	pprint.pprint(dict_orthologs)
	print(len(dict_paralogs))
	pprint.pprint(dict_paralogs)	
	
	orthologOutput='SO_'+scrmshawOutput
	with open(scrmshawOutputPath,'r') as so,open(orthologOutput,'w') as fo:
		for line in so:
			cols=line.split('\t')
			
			#Check col 5 and col 10
			#To reduce calculation
			if cols[5]==cols[10]:
				#if one gene 
				if cols[5].find(',') == -1:
					#In theory shouldnt enter here ever
					#check if col 5 has ortholog or paralog or both			
					if (cols[5] in dict_orthologs) and (cols[5] not in dict_paralogs):
						ortho_para1 = dict_orthologs[cols[5]]+'(o)'
					if (cols[5] in dict_paralogs) and (cols[5] not in dict_orthologs):
						ortho_para1 = dict_paralogs[cols[5]]+'(p)'
					if (cols[5] not in dict_paralogs) and (cols[5] not in dict_orthologs):
						ortho_para1 = '-'
					if (cols[5] in dict_paralogs) and (cols[5] in dict_orthologs):
						#check they both are similar or not
						if dict_orthologs[cols[5]] != dict_paralogs[cols[5]].strip(','):
							ortho_para1 = dict_orthologs[cols[5]]+'(o),'+ dict_paralogs[cols[5]]+'(p)'
						else:
							ortho_para1 = dict_orthologs[cols[5]]+'(op)'
			
					fo.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+ortho_para1+'\t'+cols[7]+'\t'+cols[8]+'\t'+cols[9]+'\t'+cols[10]+'\t'+ortho_para1+'\t'+cols[12]+'\t'+cols[13]+'\t'+cols[14]+'\t'+cols[15]+'\t'+cols[16]+'\t'+cols[17])

				#if multiple genes
				
				if cols[5].find(',') != -1:
					ortho_para1=''
					genes1=cols[5].split(',')	
					for i in range(0,len(genes1)):
				
						#check if col 5 1st gene has ortholog or paralog or both			
						if (genes1[i] in dict_orthologs) and (genes1[i] not in dict_paralogs):
							ortho_para1 = ortho_para1+ dict_orthologs[genes1[i]]+'(o),'
						if (genes1[i] in dict_paralogs) and (genes1[i] not in dict_orthologs):
							ortho_para1 = ortho_para1+ dict_paralogs[genes1[i]]+'(p),'
						if (genes1[i] not in dict_paralogs) and (genes1[i] not in dict_orthologs):
							ortho_para1 = ortho_para1+ '-'
						if (genes1[i] in dict_paralogs) and (genes1[i] in dict_orthologs):
							#check they both are similar or not
							if dict_orthologs[genes1[i]] != dict_paralogs[genes1[i]].strip(','):
								ortho_para1 = ortho_para1+ dict_orthologs[genes1[i]]+'(o),'+ dict_paralogs[genes1[i]]+'(p),'
							else:
								ortho_para1 = ortho_para1+ dict_orthologs[genes1[i]]+'(op),'
					
					fo.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+ortho_para1+'\t'+cols[7]+'\t'+cols[8]+'\t'+cols[9]+'\t'+cols[10]+'\t'+ortho_para1+'\t'+cols[12]+'\t'+cols[13]+'\t'+cols[14]+'\t'+cols[15]+'\t'+cols[16]+'\t'+cols[17])
			
			
			
			elif cols[5] != cols[10]:
			
				#CASE A
				#Neither col 5 nor col 10 has more than one gene on it---
				if cols[5].find(',') == -1 and cols[10].find(',') == -1:
					#check if col 5 has ortholog or paralog or both			
					if (cols[5] in dict_orthologs) and (cols[5] not in dict_paralogs):
						ortho_para1 = dict_orthologs[cols[5]]+'(o)'
					if (cols[5] in dict_paralogs) and (cols[5] not in dict_orthologs):
						ortho_para1 = dict_paralogs[cols[5]]+'(p)'
					if (cols[5] not in dict_paralogs) and (cols[5] not in dict_orthologs):
						ortho_para1 = '-'
					if (cols[5] in dict_paralogs) and (cols[5] in dict_orthologs):
						#check they both are similar or not
						if dict_orthologs[cols[5]] != dict_paralogs[cols[5]].strip(','):
							ortho_para1 = dict_orthologs[cols[5]]+'(o),'+ dict_paralogs[cols[5]]+'(p)'
						else:
							ortho_para1 = dict_orthologs[cols[5]]+'(op)'
			
					#fo.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+ortho_para+'\t'+cols[7]+'\t'+cols[8]+'\t'+cols[9]+'\t'+cols[10]+'\t'+cols[11]+'\t'+cols[12]+'\t'+cols[13]+'\t'+cols[14]+'\t'+cols[15]+'\t'+cols[16]+'\t'+cols[17])

					#check if col 10 has ortholog or paralog or both	
					if (cols[10] in dict_orthologs) and (cols[10] not in dict_paralogs):
						ortho_para2 = dict_orthologs[cols[10]]+'(o)'
					if (cols[10] in dict_paralogs) and (cols[10] not in dict_orthologs):
						ortho_para2 = dict_paralogs[cols[10]]+'(p)'
					if (cols[10] not in dict_paralogs) and (cols[10] not in dict_orthologs):
						ortho_para2 = '-'
					if (cols[10] in dict_paralogs) and (cols[10] in dict_orthologs):
						#check they both are similar or not
						if dict_orthologs[cols[10]] != dict_paralogs[cols[10]].strip(','):
							ortho_para2 = dict_orthologs[cols[10]]+'(o),'+ dict_paralogs[cols[10]]+'(p)'
						else:
							ortho_para2 = dict_orthologs[cols[10]]+'(op)'
			
					fo.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+ortho_para1+'\t'+cols[7]+'\t'+cols[8]+'\t'+cols[9]+'\t'+cols[10]+'\t'+ortho_para2+'\t'+cols[12]+'\t'+cols[13]+'\t'+cols[14]+'\t'+cols[15]+'\t'+cols[16]+'\t'+cols[17])


				#CASE B
				#Column 5 has more than 1 gene but col 10 has 1 gene only
				elif cols[5].find(',') != -1 and cols[10].find(',') == -1:
					ortho_para1=''
					genes1=cols[5].split(',')	
					for i in range(0,len(genes1)):
				
						#check if col 5 1st gene has ortholog or paralog or both			
						if (genes1[i] in dict_orthologs) and (genes1[i] not in dict_paralogs):
							ortho_para1 = ortho_para1+ dict_orthologs[genes1[i]]+'(o),'
						if (genes1[i] in dict_paralogs) and (genes1[i] not in dict_orthologs):
							ortho_para1 = ortho_para1+ dict_paralogs[genes1[i]]+'(p),'
						if (genes1[i] not in dict_paralogs) and (genes1[i] not in dict_orthologs):
							ortho_para1 = ortho_para1+ '-'
						if (genes1[i] in dict_paralogs) and (genes1[i] in dict_orthologs):
							#check they both are similar or not
							if dict_orthologs[genes1[i]] != dict_paralogs[genes1[i]].strip(','):
								ortho_para1 = ortho_para1+ dict_orthologs[genes1[i]]+'(o),'+ dict_paralogs[genes1[i]]+'(p),'
							else:
								ortho_para1 = ortho_para1+ dict_orthologs[genes1[i]]+'(op),'
			
					
					#check if col 10 has ortholog or paralog or both	
					if (cols[10] in dict_orthologs) and (cols[10] not in dict_paralogs):
						ortho_para2 = dict_orthologs[cols[10]]+'(o)'
					if (cols[10] in dict_paralogs) and (cols[10] not in dict_orthologs):
						ortho_para2 = dict_paralogs[cols[10]]+'(p)'
					if (cols[10] not in dict_paralogs) and (cols[10] not in dict_orthologs):
						ortho_para2 = '-'
					if (cols[10] in dict_paralogs) and (cols[10] in dict_orthologs):
						#check they both are similar or not
						if dict_orthologs[cols[10]] != dict_paralogs[cols[10]].strip(','):
							ortho_para2 = dict_orthologs[cols[10]]+'(o),'+ dict_paralogs[cols[10]]+'(p)'
						else:
							ortho_para2 = dict_orthologs[cols[10]]+'(op)'
			
					fo.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+ortho_para1+'\t'+cols[7]+'\t'+cols[8]+'\t'+cols[9]+'\t'+cols[10]+'\t'+ortho_para2+'\t'+cols[12]+'\t'+cols[13]+'\t'+cols[14]+'\t'+cols[15]+'\t'+cols[16]+'\t'+cols[17])


				#CASE C
	# 			#Column 5 has 1 gene but col 10 has more than 1 gene 
	#
				elif cols[5].find(',') == -1 and cols[10].find(',') != -1:
			
					if (cols[5] in dict_orthologs) and (cols[5] not in dict_paralogs):
						ortho_para1 = dict_orthologs[cols[5]]+'(o)'
					if (cols[5] in dict_paralogs) and (cols[5] not in dict_orthologs):
						ortho_para1 = dict_paralogs[cols[5]]+'(p)'
					if (cols[5] not in dict_paralogs) and (cols[5] not in dict_orthologs):
						ortho_para1 = '-'
					if (cols[5] in dict_paralogs) and (cols[5] in dict_orthologs):
						#check they both are similar or not
						if dict_orthologs[cols[5]] != dict_paralogs[cols[5]].strip(','):
							ortho_para1 = dict_orthologs[cols[5]]+'(o),'+ dict_paralogs[cols[5]]+'(p)'
						else:
							ortho_para1 = dict_orthologs[cols[5]]+'(op)'
			
			
				
					ortho_para2=''
					genes2=cols[10].split(',')	
					for i in range(0,len(genes2)):
				
						#check if col 5 1st gene has ortholog or paralog or both			
						if (genes2[i] in dict_orthologs) and (genes2[i] not in dict_paralogs):
							ortho_para2 = ortho_para2+ dict_orthologs[genes2[i]]+'(o),'
						if (genes2[i] in dict_paralogs) and (genes2[i] not in dict_orthologs):
							ortho_para2 = ortho_para2+ dict_paralogs[genes2[i]]+'(p),'
						if (genes2[i] not in dict_paralogs) and (genes2[i] not in dict_orthologs):
							ortho_para2 = ortho_para2+ '-'
						if (genes2[i] in dict_paralogs) and (genes2[i] in dict_orthologs):
							#check they both are similar or not
							if dict_orthologs[genes2[i]] != dict_paralogs[genes2[i]].strip(','):
								ortho_para2 = ortho_para2+ dict_orthologs[genes2[i]]+'(o),'+ dict_paralogs[genes2[i]]+'(p),'
							else:
								ortho_para2 = ortho_para2+ dict_orthologs[genes2[i]]+'(op),'

					fo.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+ortho_para1+'\t'+cols[7]+'\t'+cols[8]+'\t'+cols[9]+'\t'+cols[10]+'\t'+ortho_para2+'\t'+cols[12]+'\t'+cols[13]+'\t'+cols[14]+'\t'+cols[15]+'\t'+cols[16]+'\t'+cols[17])



				#CASE D
				#Column 5 and col 10 both have more than 1 gene  but different
				elif cols[5].find(',') != -1 and cols[10].find(',') != -1:
					ortho_para1=''
					genes1=cols[5].split(',')	
					for i in range(0,len(genes1)):
				
						#check if col 5 1st gene has ortholog or paralog or both			
						if (genes1[i] in dict_orthologs) and (genes1[i] not in dict_paralogs):
							ortho_para1 = ortho_para1+ dict_orthologs[genes1[i]]+'(o),'
						if (genes1[i] in dict_paralogs) and (genes1[i] not in dict_orthologs):
							ortho_para1 = ortho_para1+ dict_paralogs[genes1[i]]+'(p),'
						if (genes1[i] not in dict_paralogs) and (genes1[i] not in dict_orthologs):
							ortho_para1 = ortho_para1+ '-'
						if (genes1[i] in dict_paralogs) and (genes1[i] in dict_orthologs):
							#check they both are similar or not
							if dict_orthologs[genes1[i]] != dict_paralogs[genes1[i]].strip(','):
								ortho_para1 = ortho_para1+ dict_orthologs[genes1[i]]+'(o),'+ dict_paralogs[genes1[i]]+'(p),'
							else:
								ortho_para1 = ortho_para1+ dict_orthologs[genes1[i]]+'(op),'
							
				
					ortho_para2=''
					genes2=cols[10].split(',')	
					for i in range(0,len(genes2)):
				
						#check if col 5 1st gene has ortholog or paralog or both			
						if (genes2[i] in dict_orthologs) and (genes2[i] not in dict_paralogs):
							ortho_para2 = ortho_para2+ dict_orthologs[genes2[i]]+'(o),'
						if (genes2[i] in dict_paralogs) and (genes2[i] not in dict_orthologs):
							ortho_para2 = ortho_para2+ dict_paralogs[genes2[i]]+'(p),'
						if (genes2[i] not in dict_paralogs) and (genes2[i] not in dict_orthologs):
							ortho_para2 = ortho_para2+ '-'
						if (genes2[i] in dict_paralogs) and (genes2[i] in dict_orthologs):
							#check they both are similar or not
							if dict_orthologs[genes2[i]] != dict_paralogs[genes2[i]].strip(','):
								ortho_para2 = ortho_para2+ dict_orthologs[genes2[i]]+'(o),'+ dict_paralogs[genes2[i]]+'(p),'
							else:
								ortho_para2 = ortho_para2+ dict_orthologs[genes2[i]]+'(op),'

					fo.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+ortho_para1+'\t'+cols[7]+'\t'+cols[8]+'\t'+cols[9]+'\t'+cols[10]+'\t'+ortho_para2+'\t'+cols[12]+'\t'+cols[13]+'\t'+cols[14]+'\t'+cols[15]+'\t'+cols[16]+'\t'+cols[17])

#./orthologyMapping.py -np1 tcas -sp1id TCAST.idmap.txt -sp1pc TCAST.97pc.txt -brh TCAST_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet OSG2geneSet.txt -symb fb_synonym_fb_2020_02.tsv -conv true -setConv OGS3toOGS2_conversion.txt -so scrmshawOutput_peaksCalled_antennal_lobe_imm_1388_peaks.bed -sep ':' -flip TRUE
#../orthologyMapping.py -np1 agamb -sp1id AGAMB.idmap.txt -sp1pc AGAMB.97pc.txt -brh AGAMB_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet AGAM_4.9_genes -symb fb_synonym_fb_2020_02.tsv -conv false -so scrmshawOutput_peaksCalled_adult_circulatory_imm_MedianPointAmplitudeCurve_633_peaks.bed 
#./orthologyMapping.py -np1 apis -sp1id AMELL.idmap.txt -sp1pc AMELL.97pc.txt -brh AMELL_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet OSG2geneSet.txt -symb fb_synonym_fb_2020_02.tsv -conv true -setConv OSG3toOSG2_conversion.txt -so scrmshawOutput_peaksCalled_adult_circulatory_imm_MedianPointAmplitudeCurve_575_peaks.bed -sep '|' -flip false

main()
