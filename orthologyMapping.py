#!/usr/bin/env python
################################################################
# Orthology mapping-
#
# (c) Hasiba Asma July 2020
# Commented: June 2021
# this script can be used to get ortholog/paralog of species X wrt DMEL for each gene (list provided by user)
# or after post processing of SCRMshaw to get ortholog/paralog wrt DMEL
# input files can be downloaded from https://ezmeta.unige.ch/i5k/
# files needed for X and DMEL species include DMELA.idmap.txt, X.idmap.txt, DMELA.97.pc.txt, X.97.pc.txt, X_DMELA.brh, listOfgenes.txt(1 gene in a line), fb_synonym_fb_2020_02.tsv (fly symbol file)
# e.g if you have files for AGAM from ^ above link you can run the following command to get list of orthologs and paralogs
# ../orthologs.py -np1 agamb -sp1id AGAMB.idmap.txt -sp1pc AGAMB.97pc.txt -brh AGAMB_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet AGAM_4.9_genes -symb fb_synonym_fb_2020_02.tsv -conv false
# and in addition if you want to convert SCRMshaw output to file which has DMEL ortho/paralogs use -so flag and provide scrmshaw output
#updated and commented June 2021
################################################################


import os
import sys
import argparse
import pprint
import re
import csv


# This dictionary is created to save id mapping file of Spec_X. The format of this dictionary is something like this: [OFAS00001] = 'OFAS2:0001'
def idMap_dict(nameOfDict, file):
    with open(file, 'r') as fi:
        rows = (line.split(' ') for line in fi)
        nameOfDict = {row[1].strip('\n'): row[0] for row in rows}
    return (nameOfDict)


# This dictionary is created to save best reciprocal hits (DMEL_SpecX or SpecX_DMEL format). The format of this dictionary is something like this: [7227:00001] = 'OFAS2:0001'
def brh_dict(nameOfBrhDict, fileBrh, flipping):
    with open(fileBrh, 'r') as fb:
        rowsB = (line.split('\t') for line in fb)
        if (flipping == 'true'):
            # print('true')
            nameOfBrhDict = {rowb[0]: rowb[1] for rowb in rowsB}
        else:
            # print('False')
            nameOfBrhDict = {rowb[1]: rowb[0] for rowb in rowsB}
    return (nameOfBrhDict)


# This dictionary is created to save id mapping file of DMEL. The format of this dictionary is something like this: [7227:0001] = 'FBgn0264125'
def idMap2_dict(nameOfDict2, file2):
    with open(file2, 'r') as fi:
        rows = (line.split(' ') for line in fi)
        nameOfDict2 = {row[0]: row[1].strip() for row in rows}
    return (nameOfDict2)


# This function takes in the list of paralogs (via *.97.pc.txt file) and returns the LIST of list of paralogs e.g. ['OFAS2:004b90', 'OFAS2:004b2b']
def paralogs(nameoflist, filepc):
    listOflist = []
    with open(filepc, 'r') as fp:
        rows = (re.split(',|\t', line.strip()) for line in fp)
        for row in rows: listOflist.append(row)
    # for row in rows:
    # listOflist.append(row)
    return (listOflist)


# this function is used to find a value and return the index of that value
def find(value, matrix):
    for list in matrix:
        if value in list:
            return [matrix.index(list), list.index(value)]
    return -1


# this function is used to create a dictionary, when there is need to convert gene names using the file provided by user convGeneSet
def conv_dict(nameOfCDict, fileC):
    with open(fileC, 'r') as fi:
        rows = (line.split('\t') for line in fi)
        nameOfCDict = {row[1].strip(): row[0] for row in rows}
    return (nameOfCDict)


# this function creates a dictionary to store FBgn gene symbol of drosophilla
def geneSymbol_dict(nameofSDict, fileS):
    with open(fileS, 'r') as fi:
        rows = (line.split('\t') for line in fi if not (line.startswith('#') or line.strip() == ''))
        nameofSDict = {row[0]: row[2] for row in rows}
    return (nameofSDict)


# this function creates two dictionaries (ortholog and paralog) using x_final file created before..
def orthologs_dict(nameOfODict_orth, nameOfDict_para, fileO):
    nameOfODict_orth = {}
    nameOfDict_para = {}

    #modifying this by taking into account of dmel paralogs
    with open(fileO, 'r') as fi:
        rows = (line.split('\t') for line in fi)
        # nameOfCDict = {row[0]:row[1]+row[2] for row in rows}
        for row in rows:
            if row[2] != 'NULL' and row[2] != '-':
                if row[2]!= 'NoSymbolFound':
                    nameOfODict_orth[row[0]] = row[2]
                else:
                    nameOfODict_orth[row[0]] = row[1]

            if row[4] != 'NULL' and row[4] != '-':
                pOpsall = row[4].lstrip('[').rstrip(']\n').replace("'", "")
                pOps=pOpsall.split(',')
                if row[2] not in pOps:
                    nameOfDict_para[row[0]]=row[2]+pOpsall
                else:
                    nameOfDict_para[row[0]]=pOpsall

            elif (row[4] =='-' or row[4] == 'NULL') and row[1]=='NoDirectOrtholog':

                if (row[6] != 'NULL' and row[6] != '-') and (row[8] =='NULL' or row[8]=='-\n') :
                    if row[6]!= 'NoSymbolFound':
                        nameOfDict_para[row[0]] = row[6].strip()
                    else:
                        nameOfDict_para[row[0]] = row[5].strip()
                elif (row[6] != 'NULL' and row[6] != '-') and (row[8] !='NULL' and row[8]!='-\n'):
                    ppsall= row[8].lstrip('[').rstrip(']\n').replace("'","")
                    pps=ppsall.split(',')
                    if row[6].strip() not in pps:
                        nameOfDict_para[row[0]] = row[6].strip() + ppsall
                    else:
                        nameOfDict_para[row[0]] = ppsall

    return (nameOfODict_orth, nameOfDict_para)


# this function creates files using dictionaries
def file_from_dict(nameOfFile, nameOfDict):
    nameOfFileCSV = nameOfFile + '.csv'
    w = csv.writer(open(nameOfFileCSV, "w"))
    for key, val in nameOfDict.items():
        w.writerow([key, val])


# -----------------------------MAIN FUNCTION-------------------
def main():
    # global d1
    parser = argparse.ArgumentParser()
    parser.add_argument('-np1', '--namesp1', help='name of species 1', required=True)
    parser.add_argument('-sp1id', '--sp1idmap', help='species 1 id map text file', required=True)
    parser.add_argument('-sp1pc', '--sp1pc97', help='species 1 paralogs 97 pc file', required=True)
    parser.add_argument('-brh', '--brhsp1sp2', help='best reciprocal file sp1 to sp2', required=True)
    parser.add_argument('-sp2id', '--sp2idmap', help='species 2 id map text file', required=True)
    parser.add_argument('-sp2pc', '--sp2pc97', help='species 2 paralogs 97 pc file', required=True)
    parser.add_argument('-geneSet', '--geneSet', help='gene set from gff same as gene set used for orthology mapping',
                        required=True)
    parser.add_argument('-symb', '--symbolGene', help='Gene Symbol', required=True)
    parser.add_argument('-conv', '--conversion', help='if conversion req or not', default=False)
    parser.add_argument('-setConv', '--setConvGene', help='conversion file')
    parser.add_argument('-flip', '--flipped', help='if in brh DMEL:specie1 set it to False', default='True')
    parser.add_argument('-sep', '--separator', help='if separator is not colon, provide its value', default=':')
    parser.add_argument('-so', '--scrmshawOutput', help='if you want to get ortho/paralogs in scrmshaw output, provide scrmhsaw output file here, else leave this parameter',default='NoSCRM')
    args = parser.parse_args()
    # absolute path
    namesp1 = args.namesp1
    sp1idmap = os.path.abspath(args.sp1idmap)
    sp1pc97 = os.path.abspath(args.sp1pc97)
    brhsp1sp2 = os.path.abspath(args.brhsp1sp2)
    sp2idmap = os.path.abspath(args.sp2idmap)
    sp2pc97 = os.path.abspath(args.sp2pc97)
    geneSet = os.path.abspath(args.geneSet)
    symbolGene = os.path.abspath(args.symbolGene)
    conversion = args.conversion
    flipped = (args.flipped).lower()
    separator = args.separator
    scrmshawOutput = (args.scrmshawOutput)
    if (scrmshawOutput != 'NoSCRM'):
        scrmshawOutputPath = os.path.abspath(scrmshawOutput)

    # Depending on naming convention of genes used for SCRMshaw and orthoDB, you might need to convert them to same convention first
    # and if there is need of conversion (or conversion=true), user need to provide the respective file needed for conversion and that
    # file is used to create a dictionary to save them
    if (conversion == "TRUE" or conversion == "T" or conversion == "true"):
        convGeneSet = os.path.abspath(args.setConvGene)
        dict_conv = conv_dict('conv', convGeneSet)

    # creating a dictionary of SpecX using idmap file (format: )
    dict_sp1id = idMap_dict('sp1id', sp1idmap)
    # pprint.pprint(dict_sp1id)
    # print('Flipped',flipped)
    # creating a dictionary to store best reciprocal hits (format: )
    dict_brh = brh_dict('brh', brhsp1sp2, flipped)
    # creating a dictionary of DMEL using idmap file (format: )
    dict_sp2id = idMap2_dict('sp2id', sp2idmap)
    # creating a dictionary to save symbol of genes (format: )
    dict_symb = geneSymbol_dict('symbol', symbolGene)

    # pprint.pprint(dict_symb)
    # creating list of paralog lists for specX and DMEL
    pc1 = paralogs('pc1', sp1pc97)
    pc2 = paralogs('pc2', sp2pc97)

    # creating two intermediate files that will store information regarding orthologs and paralogs
    g1 = open(namesp1 + '_temp.txt', 'w')
    g1.write(
        'GeneName_TC' + '\t' + 'GeneName_OGid' + '\t' + 'Ortholog_id_symb' + '\t' +'Dmel_paralogs'+'\t'+'Dmel_paralogs_ids'+'\t'+'Dmel_paralogs_Symbols'+'\t'+ 'listOfParalogsIfAny' + '\t' + 'ifParalogHasOrtholog_id_symb' + '\t' +'Dmel_Pparalogs'+'\t'+'Dmel_Pparalogs_ids'+'\t'+'Dmel_Pparalogs_Symbols'+'\n')
    f1a = open(namesp1 + '_final.txt', 'w')
    f1a.write(
        'GeneName' + '\t' + 'Orthologs' + '\t' + 'GeneSymbolOrthologs' + '\t' +'Dmel_paralogs'+'\t'+'GeneSymbolDmel_paralogs'+'\t' +'ParalogsThatHaveOrthologs' + '\t' + 'GeneSymbolParalogs'+ '\t' +'Dmel_Pparalogs'+'\t'+'GeneSymbolDmel_Pparalogs'+'\n')

    # opening geneSet file that has 1 gene per line (all of genes extracted from its gff) and
    # finding out if there are any orthologs or/and paralogs present wrt DMEL, if present write into two temp files created before..
    # to use later on to edit SCRMshaw prediction file to add respective ortho/para data in it
    with open(geneSet, 'r') as gS, open(namesp1 + '_temp.txt', 'a') as g1, open(namesp1 + '_final.txt', 'a') as f1a:
        for line in gS:
            geneName = line.rstrip('\n')
            # g1.write(geneName)
            # f1.write(geneName+'\t')
            # if there is a step of conversion of gene Naming involved, then making sure to use the right one using the
            # dictioanary created previously
            if (conversion == "TRUE" or conversion == "T" or conversion == "true"):
                # print('require c')
                if geneName.split(separator)[1] in dict_conv:
                    # f1.write(dict_conv[geneName.split(separator)[1]]+'\t')
                    geneName2 = dict_conv[geneName.split(separator)[1]]
                else:
                    # f1.write(geneName+'\t')
                    geneName2 = geneName
            # if not, then go ahead with line one of gene name
            else:
                # f1.write(geneName+'\t')
                geneName2 = geneName

            # if gene is present in the idMap dictionary of Spec X
            if geneName in dict_sp1id:
                # print('present')
                # writing genename in temp file using idmap file
                # g1.write('\t'+dict_sp1id[geneName])

                # STEP 1: Check ORTHOLOGS
                # First of all check if it has any Ortholog or not..
                # check if it has orthologs-forget about paralogs---
                if dict_sp1id[geneName] in dict_brh:
                    # g1.write('Has ortholog')
                    if dict_brh[dict_sp1id[geneName]] in dict_sp2id:

                        #initializing all variables used for dmels paralogs
                        dmelsParalogsList=''
                        dmelsParalogsListFBgn=[]
                        dmelsParalogsListFBgnSymb=[]


                        # print(geneName2)
                        # print(dict_sp1id[geneName2])
                        # print(dict_brh[dict_sp1id[geneName]]+'__'+dict_sp2id[dict_brh[dict_sp1id[geneName]]])
                        # print("---Has Ortholog--") #probably no need to look at paralogs then, even if they are there, they wont be in brh file right?
                        # g1.write(geneName+'\t'+dict_sp1id[geneName]+'\t'+'Has parlogs:'+str('sublist')+'\t'+dict_brh[dict_sp1id[geneName]]+'\t'+dict_sp2id[dict_brh[dict_sp1id[geneName]]]+'\n')
                        g1.write(
                            geneName2 + '\t' + dict_sp1id[geneName] + '\t' + dict_brh[dict_sp1id[geneName]] + '__' +
                            dict_sp2id[dict_brh[dict_sp1id[geneName]]])
                        #	f1.write('\t'+dict_sp2id[dict_brh[dict_sp1id[geneName]]])
                        f1a.write(geneName2 + '\t' + dict_sp2id[dict_brh[dict_sp1id[geneName]]])
                        noOrigOrtholog = 'False'

                        # Gene symbol:
                        if dict_sp2id[dict_brh[dict_sp1id[geneName]]] in dict_symb:
                            g1.write('/' + dict_symb[dict_sp2id[dict_brh[dict_sp1id[geneName]]]])
                            # f1.write('\t'+dict_symb[dict_sp2id[dict_brh[dict_sp1id[geneName]]]])
                            f1a.write('\t' + dict_symb[dict_sp2id[dict_brh[dict_sp1id[geneName]]]])
                        else:
                            g1.write('/' + 'NoSymb')
                            # f1.write('\t'+'NULL')
                            f1a.write('\t' + 'NoSymbolFound')


                        #---FORGOT DMEL paralogs
                        # check if dmel's id (whose ortholog is found) has any paralogs?
                        #if dict_brh[dict_sp1id[geneName]] in pc2:
                        if any(dict_brh[dict_sp1id[geneName]] in sublist2 for sublist2 in pc2):
                            #has paralogs , printing those here
                            # find where is it located on the list
                            paralogListIndex2 = find(dict_brh[dict_sp1id[geneName]] , pc2)
                         #   print('list of dmel paralogs: ' + str(pc2[paralogListIndex2[0]]))
                            dmelsParalogsList=pc2[paralogListIndex2[0]]

                            # go through each of these paralogs and list their symbols
                            for i in range(0, len(pc2[paralogListIndex2[0]])):

                                # check each para for its symbol
                                #first check if it has Fbgn id
                                if pc2[paralogListIndex2[0]][i] in dict_sp2id:
                                    #print("Fbgn id for this paralog, ",pc2[paralogListIndex2[0]][i]," is: ",dict_sp2id[pc2[paralogListIndex2[0]][i]])
                                    dmelsParalogsListFBgn.append(dict_sp2id[pc2[paralogListIndex2[0]][i]])
                                    #check for symbols
                                    if dict_sp2id[pc2[paralogListIndex2[0]][i]] in dict_symb:
                                       # print("and the symbol is: ", dict_symb[dict_sp2id[pc2[paralogListIndex2[0]][i]]])
                                        dmelsParalogsListFBgnSymb.append(dict_symb[dict_sp2id[pc2[paralogListIndex2[0]][i]]])
                                    else:
                                     #   print("but doesnt have symbol")
                                        dmelsParalogsListFBgnSymb.append('-')
                                else:
                                   # print("no FbgnID found")
                                    dmelsParalogsListFBgn.append('-')
                        if dmelsParalogsList!='':
                           # print(dmelsParalogsList)
                            #print(dmelsParalogsListFBgn)
                            #print(dmelsParalogsListFBgnSymb)
                            g1.write('\t'+str(dmelsParalogsList)+'\t'+str(dmelsParalogsListFBgn)+'\t'+str(dmelsParalogsListFBgnSymb))
                            f1a.write('\t'+str(dmelsParalogsListFBgn)+'\t'+str(dmelsParalogsListFBgnSymb))
                        else:
                            g1.write('\t'+'NoParalogs'+'\t'+'-'+'\t'+'-')
                            f1a.write('\t'+'NoParalogs' + '\t' + '-')
                        #end-DMELs paralogs for direct orthologs


                    else:  # id is not found in dictionary of ids--so no ortholog can be found
                        # print('id is not found in dictionary of ids')
                        g1.write(geneName2 + '\t' + dict_sp1id[geneName2] + '\t' + 'NULL' + '\t' + 'NULL'+'\t' + 'NULL' + '\t' + 'NULL'+ '\t' + 'NULL'+'\t' + 'NULL' + '\t' + 'NULL')
                        f1a.write(geneName2 + '\t' + 'NULL' + '\t' + 'NULL'+'\t' + 'NULL' + '\t' + 'NULL'+'\t'+'NULL' + '\t' + 'NULL')

                else:  # if not present in BRH file
                    # print(dict_sp1id[geneName])
                    # print("---NO Ortholog--")
                    g1.write(geneName2 + '\t' + dict_sp1id[geneName] + '\t' + 'NO ortholog'+'\t'+'-'+'\t'+'-'+'\t'+'-')
                    f1a.write(geneName2 + '\t' + 'NoDirectOrtholog' + '\t' + 'NULL'+'\t'+'-'+'\t'+'-')
                    noOrigOrtholog = 'True'

                # STEP 2: Check PARALOGS.. (if no ortholog found)
                # SECOND STEP --check if there is any paralog present, if so: whether that paralog has any ortholog [although it won't if gene itself has already an ortholog]
                # using the id of geneName, check to see if its present in paralog file (using pc list created before)

                # if there wasnt any ortholog then go ahead and check its paralogs
                if noOrigOrtholog == 'True':
                    # try start
                    if any(dict_sp1id[geneName] in sublist for sublist in pc1):
                        # find where is it located on the list
                        paralogListIndex = find(dict_sp1id[geneName], pc1)

                        # if this index [1] item is 0 : that means this does have paralogs, but none of paralogs have any ortholog
                        # if  this index [1] item is 1: that means the paralog at position 0 might have ortholog might not > check if it does
                        if paralogListIndex[1] == 0:
                            #	print("No Ortholog of its following paralogs",str(pc1[paralogListIndex[0]]))
                            f1a.write('\t' + 'HasParalogButNoChanceOfOrthologOfParalogs' + '\t' + '-'+'\t'+'-'+'\t'+'-')
                            g1.write('\t' + str(pc1[paralogListIndex[0]]) + '\t' + 'NULL'+'\t'+'-'+'\t'+'-'+'\t'+'-')
                        #paralog can only have ortholog if its not the first one in 97 pc file
                        elif paralogListIndex[1] != 0:
                            # print("Its paralog may have an ortholog",str(pc1[paralogListIndex[0]]))

                            if pc1[paralogListIndex[0]][0] in dict_brh:  # its paralog does have  ortholog
                                #	print("Yes its paralog DOES have an ortholog- see below")
                                if dict_brh[pc1[paralogListIndex[0]][0]] in dict_sp2id:

                                    # initializing all variables used for dmels paralogs
                                    dmelsParalogsParalogsList = ''
                                    dmelsParalogsParalogsListFBgn = []
                                    dmelsParalogsParalogsListFBgnSymb = []

                                    #	print(pc1[paralogListIndex[0]][0]+'--'+dict_brh[pc1[paralogListIndex[0]][0]]+'__'+dict_sp2id[dict_brh[pc1[paralogListIndex[0]][0]]])
                                    g1.write('\t' + str(pc1[paralogListIndex[0]]) + '\t' + pc1[paralogListIndex[0]][
                                        0] + '--' + dict_brh[pc1[paralogListIndex[0]][0]] + '__' + dict_sp2id[
                                                 dict_brh[pc1[paralogListIndex[0]][0]]])
                                    f1a.write('\t' + dict_sp2id[dict_brh[pc1[paralogListIndex[0]][0]]] + '\t')
                                    # check its symbol
                                    if dict_sp2id[dict_brh[pc1[paralogListIndex[0]][0]]] in dict_symb:
                                        g1.write('/' + dict_symb[dict_sp2id[dict_brh[pc1[paralogListIndex[0]][0]]]])
                                        f1a.write(dict_symb[dict_sp2id[dict_brh[pc1[paralogListIndex[0]][0]]]])
                                    else:
                                        g1.write('/' + 'NoSymb')
                                        f1a.write('NoSymbolFound')


                                    # ---FORGOT DMEL paralogs
                                    # check if dmel's id (whose ortholog is found) has any paralogs?
                                    # if dict_brh[dict_sp1id[geneName]] in pc2:
                                    if any(dict_brh[pc1[paralogListIndex[0]][0]] in sublist2 for sublist2 in pc2):
                                        # has paralogs , printing those here
                                        # find where is it located on the list
                                        paralogListIndex3 = find(dict_brh[pc1[paralogListIndex[0]][0]], pc2)
                                    #    print('P list of dmel paralogs: ' + str(pc2[paralogListIndex3[0]]))
                                        dmelsParalogsParalogsList = pc2[paralogListIndex3[0]]
                                        # go through each of these paralogs and list their symbols
                                        for i in range(0, len(pc2[paralogListIndex3[0]])):

                                            # check each para for its symbol
                                            # first check if it has Fbgn id
                                            if pc2[paralogListIndex3[0]][i] in dict_sp2id:
                                               # print("Fbgn id for this paralog, ", pc2[paralogListIndex3[0]][i],
                                                #      " is: ", dict_sp2id[pc2[paralogListIndex3[0]][i]])
                                                dmelsParalogsParalogsListFBgn.append(dict_sp2id[pc2[paralogListIndex3[0]][i]])
                                                # check for symbols
                                                if dict_sp2id[pc2[paralogListIndex3[0]][i]] in dict_symb:
                                                   # print("and the symbol is: ",
                                                       #   dict_symb[dict_sp2id[pc2[paralogListIndex3[0]][i]]])
                                                    dmelsParalogsParalogsListFBgnSymb.append(
                                                        dict_symb[dict_sp2id[pc2[paralogListIndex3[0]][i]]])
                                                else:
                                                    #print("but doesnt have symbol")
                                                    dmelsParalogsParalogsListFBgnSymb.append('-')
                                            else:
                                                # print("no FbgnID found")
                                                dmelsParalogsParalogsListFBgn.append('-')
                                    if dmelsParalogsParalogsList != '':
                                        # print(dmelsParalogsList)
                                        # print(dmelsParalogsListFBgn)
                                        # print(dmelsParalogsListFBgnSymb)
                                        g1.write('\t' + str(dmelsParalogsParalogsList) + '\t' + str(
                                            dmelsParalogsParalogsListFBgn) + '\t' + str(dmelsParalogsParalogsListFBgnSymb))
                                        f1a.write('\t' + str(dmelsParalogsParalogsListFBgn) + '\t' + str(
                                            dmelsParalogsParalogsListFBgnSymb))
                                    else:
                                        g1.write('\t' + 'NoParalogsParalogs' + '\t' + '-' + '\t' + '-')
                                        f1a.write('\t' + 'NoParalogsParalogs' + '\t' + '-')
                                    # end-DMELs paralogs for paralogs orthologs

                                    # end-DMELs









                            else:  # its paralog also doesnt have any ortholog
                                # print("no its paralog also doesnt have any ortholog")
                                g1.write('\t' + str(pc1[paralogListIndex[0]]) + '\t' + '-'+'\t'+'-'+'\t'+'-'+'\t'+'-')
                                f1a.write('\t' + 'HasParalogButNoOrthologOfParalogs' + '\t' + '-'+'\t'+'-'+'\t'+'-')

                    else:  # not found in 97 pc file
                        # print(" doesnt have paralogs")
                        f1a.write('\t' + 'NoParalogsFound' + '\t' + '-'+'\t'+'-'+'\t'+'-')
                        g1.write('\t' + 'NoParalogs' + '\t' + '-'+'\t'+'-'+'\t'+'-'+'\t'+'-')

                else:  # ortholog is present, no need to look at paralog
                    #	print(" doesnt need paralogs")
                    f1a.write('\t' + 'NoNeedOfParalogs' + '\t' + '-'+'\t'+'-'+'\t'+'-')
                    g1.write('\t' + '-' + '\t' + '-'+'\t'+'-'+'\t'+'-'+'\t'+'-')

            else:  # gene id not found in idmap--nothing else could be looked up
                g1.write(geneName2 + '\t' + '-' + '\t' + '-' + '\t' + '-' + '\t' + '-'+ '\t' + '-' + '\t' + '-'+'\t'+'-'+'\t'+'-'+'\t'+'-'+'\t'+'-')
                f1a.write(geneName + '\t' + 'geneIDnotMapped' + '\t' + '-' + '\t' + '-' + '\t' + '-'+ '\t' + '-' + '\t' + '-'+'\t'+'-'+'\t'+'-')
            # for sublist in pc1:
            # for item in sublist:
            g1.write('\n')
            f1a.write('\n')
    g1.close()
    f1a.close()

    fof = os.path.abspath(namesp1 + '_final.txt')
    dict_orthologs, dict_paralogs = orthologs_dict('orthologs', 'paralogs', fof)
    print("number of orthologs found:" + str(len(dict_orthologs)))
    # pprint.pprint(dict_orthologs)
    print("number of paralogs found:" + str(len(dict_paralogs)))
    # pprint.pprint(dict_paralogs)

    # creating csv files for ortholog and paralog list
    file_from_dict(namesp1 + '_orthologList', dict_orthologs)
    file_from_dict(namesp1 + '_paralogList', dict_paralogs)

    orthologOutput = 'SO_' + scrmshawOutput

    #need to fix it at some point, exclude anything useless-- op stuff (June 22)
    if (scrmshawOutput != 'NoSCRM'):
        with open(scrmshawOutputPath, 'r') as so, open(orthologOutput, 'w') as fo:
            for line in so:
                cols = line.split('\t')

                # Check col 5 and col 10
                # To reduce calculation
                if cols[5] == cols[10]:
                    # if one gene
                    if cols[5].find(',') == -1:
                        # In theory shouldnt enter here ever
                        # check if col 5 has ortholog or paralog or both
                        if (cols[5] in dict_orthologs) and (cols[5] not in dict_paralogs):
                            ortho_para1 = dict_orthologs[cols[5]] + '_o_'
                        if (cols[5] in dict_paralogs) and (cols[5] not in dict_orthologs):
                            ortho_para1 = str(dict_paralogs[cols[5]]) + '_p_'
                        if (cols[5] not in dict_paralogs) and (cols[5] not in dict_orthologs):
                            ortho_para1 = '-'
                        if (cols[5] in dict_paralogs) and (cols[5] in dict_orthologs):
                            # check they both are similar or not
                            if dict_orthologs[cols[5]] != dict_paralogs[cols[5]].strip(','):
                                ortho_para1 = dict_orthologs[cols[5]] + '_o_,' + str(dict_paralogs[cols[5]]) + '_p_'
                            else:
                                ortho_para1 = dict_orthologs[cols[5]] + '_op_'

                        fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
                            5] + '\t' + ortho_para1 + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
                                     10] + '\t' + ortho_para1 + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' +
                                 cols[15] + '\t' + cols[16] + '\t' + cols[17])

                    # if multiple genes

                    if cols[5].find(',') != -1:
                        ortho_para1 = ''
                        genes1 = cols[5].split(',')
                        for i in range(0, len(genes1)):

                            # check if col 5 1st gene has ortholog or paralog or both
                            if (genes1[i] in dict_orthologs) and (genes1[i] not in dict_paralogs):
                                ortho_para1 = ortho_para1 + dict_orthologs[genes1[i]] + '_o_,'
                            if (genes1[i] in dict_paralogs) and (genes1[i] not in dict_orthologs):
                                ortho_para1 = ortho_para1 + str(dict_paralogs[genes1[i]]) + '_p_,'
                            if (genes1[i] not in dict_paralogs) and (genes1[i] not in dict_orthologs):
                                ortho_para1 = ortho_para1 + '-'
                            if (genes1[i] in dict_paralogs) and (genes1[i] in dict_orthologs):
                                # check they both are similar or not
                                if dict_orthologs[genes1[i]] != dict_paralogs[genes1[i]].strip(','):
                                    ortho_para1 = ortho_para1 + dict_orthologs[genes1[i]] + '_o_,' + dict_paralogs[
                                        genes1[i]] + '_p_,'
                                else:
                                    ortho_para1 = ortho_para1 + dict_orthologs[genes1[i]] + '_op_,'

                        fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
                            5] + '\t' + ortho_para1 + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
                                     10] + '\t' + ortho_para1 + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' +
                                 cols[15] + '\t' + cols[16] + '\t' + cols[17])



                elif cols[5] != cols[10]:

                    # CASE A
                    # Neither col 5 nor col 10 has more than one gene on it---
                    if cols[5].find(',') == -1 and cols[10].find(',') == -1:
                        # check if col 5 has ortholog or paralog or both
                        if (cols[5] in dict_orthologs) and (cols[5] not in dict_paralogs):
                            ortho_para1 = dict_orthologs[cols[5]] + '_o_'
                        if (cols[5] in dict_paralogs) and (cols[5] not in dict_orthologs):
                            ortho_para1 = str(dict_paralogs[cols[5]]) + '_p_'
                        if (cols[5] not in dict_paralogs) and (cols[5] not in dict_orthologs):
                            ortho_para1 = '-'
                        if (cols[5] in dict_paralogs) and (cols[5] in dict_orthologs):
                            # check they both are similar or not
                            if dict_orthologs[cols[5]] != dict_paralogs[cols[5]].strip(','):
                                ortho_para1 = dict_orthologs[cols[5]] + '_o_,' + dict_paralogs[cols[5]] + '_p_'
                            else:
                                ortho_para1 = dict_orthologs[cols[5]] + '_op_'

                        # fo.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+ortho_para+'\t'+cols[7]+'\t'+cols[8]+'\t'+cols[9]+'\t'+cols[10]+'\t'+cols[11]+'\t'+cols[12]+'\t'+cols[13]+'\t'+cols[14]+'\t'+cols[15]+'\t'+cols[16]+'\t'+cols[17])

                        # check if col 10 has ortholog or paralog or both
                        if (cols[10] in dict_orthologs) and (cols[10] not in dict_paralogs):
                            ortho_para2 = dict_orthologs[cols[10]] + '_o_'
                        if (cols[10] in dict_paralogs) and (cols[10] not in dict_orthologs):
                            ortho_para2 = str(dict_paralogs[cols[10]]) + '_p_'
                        if (cols[10] not in dict_paralogs) and (cols[10] not in dict_orthologs):
                            ortho_para2 = '-'
                        if (cols[10] in dict_paralogs) and (cols[10] in dict_orthologs):
                            # check they both are similar or not
                            if dict_orthologs[cols[10]] != dict_paralogs[cols[10]].strip(','):
                                ortho_para2 = dict_orthologs[cols[10]] + '_o_,' + dict_paralogs[cols[10]] + '_p_'
                            else:
                                ortho_para2 = dict_orthologs[cols[10]] + '_op_'

                        fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
                            5] + '\t' + ortho_para1 + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
                                     10] + '\t' + ortho_para2 + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' +
                                 cols[15] + '\t' + cols[16] + '\t' + cols[17])


                    # CASE B
                    # Column 5 has more than 1 gene but col 10 has 1 gene only
                    elif cols[5].find(',') != -1 and cols[10].find(',') == -1:
                        ortho_para1 = ''
                        genes1 = cols[5].split(',')
                        for i in range(0, len(genes1)):

                            # check if col 5 1st gene has ortholog or paralog or both
                            if (genes1[i] in dict_orthologs) and (genes1[i] not in dict_paralogs):
                                ortho_para1 = ortho_para1 + dict_orthologs[genes1[i]] + '_o_,'
                            if (genes1[i] in dict_paralogs) and (genes1[i] not in dict_orthologs):
                                ortho_para1 = ortho_para1 + str(dict_paralogs[genes1[i]]) + '_p_,'
                            if (genes1[i] not in dict_paralogs) and (genes1[i] not in dict_orthologs):
                                ortho_para1 = ortho_para1 + '-'
                            if (genes1[i] in dict_paralogs) and (genes1[i] in dict_orthologs):
                                # check they both are similar or not
                                if dict_orthologs[genes1[i]] != dict_paralogs[genes1[i]].strip(','):
                                    ortho_para1 = ortho_para1 + dict_orthologs[genes1[i]] + '_o_,' + dict_paralogs[
                                        genes1[i]] + '_p_,'
                                else:
                                    ortho_para1 = ortho_para1 + dict_orthologs[genes1[i]] + '_op_,'

                        # check if col 10 has ortholog or paralog or both
                        if (cols[10] in dict_orthologs) and (cols[10] not in dict_paralogs):
                            ortho_para2 = dict_orthologs[cols[10]] + '_o_'
                        if (cols[10] in dict_paralogs) and (cols[10] not in dict_orthologs):
                            ortho_para2 = str(dict_paralogs[cols[10]]) + '_p_'
                        if (cols[10] not in dict_paralogs) and (cols[10] not in dict_orthologs):
                            ortho_para2 = '-'
                        if (cols[10] in dict_paralogs) and (cols[10] in dict_orthologs):
                            # check they both are similar or not
                            if dict_orthologs[cols[10]] != dict_paralogs[cols[10]].strip(','):
                                ortho_para2 = dict_orthologs[cols[10]] + '_o_,' + dict_paralogs[cols[10]] + '_p_'
                            else:
                                ortho_para2 = dict_orthologs[cols[10]] + '_op_'

                        fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
                            5] + '\t' + ortho_para1 + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
                                     10] + '\t' + ortho_para2 + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' +
                                 cols[15] + '\t' + cols[16] + '\t' + cols[17])


                    # CASE C
                    # 			#Column 5 has 1 gene but col 10 has more than 1 gene
                    #
                    elif cols[5].find(',') == -1 and cols[10].find(',') != -1:

                        if (cols[5] in dict_orthologs) and (cols[5] not in dict_paralogs):
                            ortho_para1 = dict_orthologs[cols[5]] + '_o_'
                        if (cols[5] in dict_paralogs) and (cols[5] not in dict_orthologs):
                            ortho_para1 = str(dict_paralogs[cols[5]]) + '_p_'
                        if (cols[5] not in dict_paralogs) and (cols[5] not in dict_orthologs):
                            ortho_para1 = '-'
                        if (cols[5] in dict_paralogs) and (cols[5] in dict_orthologs):
                            # check they both are similar or not
                            if dict_orthologs[cols[5]] != dict_paralogs[cols[5]].strip(','):
                                ortho_para1 = dict_orthologs[cols[5]] + '_o_,' + dict_paralogs[cols[5]] + '_p_'
                            else:
                                ortho_para1 = dict_orthologs[cols[5]] + '_op_'

                        ortho_para2 = ''
                        genes2 = cols[10].split(',')
                        for i in range(0, len(genes2)):

                            # check if col 5 1st gene has ortholog or paralog or both
                            if (genes2[i] in dict_orthologs) and (genes2[i] not in dict_paralogs):
                                ortho_para2 = ortho_para2 + dict_orthologs[genes2[i]] + '_o_,'
                            if (genes2[i] in dict_paralogs) and (genes2[i] not in dict_orthologs):
                                ortho_para2 = ortho_para2 + str(dict_paralogs[genes2[i]]) + '_p_,'
                            if (genes2[i] not in dict_paralogs) and (genes2[i] not in dict_orthologs):
                                ortho_para2 = ortho_para2 + '-'
                            if (genes2[i] in dict_paralogs) and (genes2[i] in dict_orthologs):
                                # check they both are similar or not
                                if dict_orthologs[genes2[i]] != dict_paralogs[genes2[i]].strip(','):
                                    ortho_para2 = ortho_para2 + dict_orthologs[genes2[i]] + '_o_,' + dict_paralogs[
                                        genes2[i]] + '_p_,'
                                else:
                                    ortho_para2 = ortho_para2 + dict_orthologs[genes2[i]] + '_op_,'

                        fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
                            5] + '\t' + ortho_para1 + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
                                     10] + '\t' + ortho_para2 + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' +
                                 cols[15] + '\t' + cols[16] + '\t' + cols[17])



                    # CASE D
                    # Column 5 and col 10 both have more than 1 gene  but different
                    elif cols[5].find(',') != -1 and cols[10].find(',') != -1:
                        ortho_para1 = ''
                        genes1 = cols[5].split(',')
                        for i in range(0, len(genes1)):

                            # check if col 5 1st gene has ortholog or paralog or both
                            if (genes1[i] in dict_orthologs) and (genes1[i] not in dict_paralogs):
                                ortho_para1 = ortho_para1 + dict_orthologs[genes1[i]] + '_o_,'
                            if (genes1[i] in dict_paralogs) and (genes1[i] not in dict_orthologs):
                                ortho_para1 = ortho_para1 + str(dict_paralogs[genes1[i]]) + '_p_,'
                            if (genes1[i] not in dict_paralogs) and (genes1[i] not in dict_orthologs):
                                ortho_para1 = ortho_para1 + '-'
                            if (genes1[i] in dict_paralogs) and (genes1[i] in dict_orthologs):
                                # check they both are similar or not
                                if dict_orthologs[genes1[i]] != dict_paralogs[genes1[i]].strip(','):
                                    ortho_para1 = ortho_para1 + dict_orthologs[genes1[i]] + '_o_,' + dict_paralogs[
                                        genes1[i]] + '_p_,'
                                else:
                                    ortho_para1 = ortho_para1 + dict_orthologs[genes1[i]] + '_op_,'

                        ortho_para2 = ''
                        genes2 = cols[10].split(',')
                        for i in range(0, len(genes2)):

                            # check if col 5 1st gene has ortholog or paralog or both
                            if (genes2[i] in dict_orthologs) and (genes2[i] not in dict_paralogs):
                                ortho_para2 = ortho_para2 + dict_orthologs[genes2[i]] + '_o_,'
                            if (genes2[i] in dict_paralogs) and (genes2[i] not in dict_orthologs):
                                ortho_para2 = ortho_para2 + str(dict_paralogs[genes2[i]]) + '_p_,'
                            if (genes2[i] not in dict_paralogs) and (genes2[i] not in dict_orthologs):
                                ortho_para2 = ortho_para2 + '-'
                            if (genes2[i] in dict_paralogs) and (genes2[i] in dict_orthologs):
                                # check they both are similar or not
                                if dict_orthologs[genes2[i]] != dict_paralogs[genes2[i]].strip(','):
                                    ortho_para2 = ortho_para2 + dict_orthologs[genes2[i]] + '_o_,' + dict_paralogs[
                                        genes2[i]] + '_p_,'
                                else:
                                    ortho_para2 = ortho_para2 + dict_orthologs[genes2[i]] + '_op_,'

                        fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
                            5] + '\t' + ortho_para1 + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
                                     10] + '\t' + ortho_para2 + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' +
                                 cols[15] + '\t' + cols[16] + '\t' + cols[17])


# ./orthologyMapping.py -np1 tcas -sp1id TCAST.idmap.txt -sp1pc TCAST.97pc.txt -brh TCAST_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet OSG2geneSet.txt -symb fb_synonym_fb_2020_02.tsv -conv true -setConv OGS3toOGS2_conversion.txt -so scrmshawOutput_peaksCalled_antennal_lobe_imm_1388_peaks.bed -sep ':' -flip TRUE
# ../orthologyMapping.py -np1 agamb -sp1id AGAMB.idmap.txt -sp1pc AGAMB.97pc.txt -brh AGAMB_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet AGAM_4.9_genes -symb fb_synonym_fb_2020_02.tsv -conv false -so scrmshawOutput_peaksCalled_adult_circulatory_imm_MedianPointAmplitudeCurve_633_peaks.bed
# ./orthologyMapping.py -np1 apis -sp1id AMELL.idmap.txt -sp1pc AMELL.97pc.txt -brh AMELL_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet OSG2geneSet.txt -symb fb_synonym_fb_2020_02.tsv -conv true -setConv OSG3toOSG2_conversion.txt -so scrmshawOutput_peaksCalled_adult_circulatory_imm_MedianPointAmplitudeCurve_575_peaks.bed -sep '|' -flip false

main()
