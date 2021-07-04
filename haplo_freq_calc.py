from numpy.lib.function_base import append
import pyslim
import tskit
import numpy as np
import pandas as pd
from pandas import DataFrame
import getopt, sys
import re

path = None

#Commandline options

full_cmd_arguments = sys.argv
argument_list = full_cmd_arguments[1:]
short_options = "hp:"
long_options = ["help", "path="]
try:
    arguments, values = getopt.getopt(argument_list, short_options, long_options)
except getopt.error as err:
    print (str(err))
    sys.exit(2)
for current_argument, current_value in arguments:
    if current_argument in ("-h", "--help"):
        print("-p/--path specifies the path for the .trees-files (add slash or backslash at the end)\n-n/--nHaplogroups defines the number of Haplogroups\n-t/--time defines the time in how many generations before admixture haplogroups should be defined")
    elif current_argument in ("-p", "--path"):
        path = str(current_value)


def haploFreqCalc(input):

    #Importing and simplifying tree sequence

    treeseq_raw = pyslim.load(input)
    all_nodes = treeseq_raw.samples()
    keep_nodes_F = []
    keep_nodes_M = []

    for u in all_nodes:
        if treeseq_raw.mutation_at(u,0) != -1 and treeseq_raw.individual(treeseq_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_FEMALE:
            keep_nodes_F.append(u)
        if treeseq_raw.mutation_at(u,100) != -1 and treeseq_raw.individual(treeseq_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_MALE:
            keep_nodes_M.append(u)
        else:
            pass

    treeseq_mtDNA = treeseq_raw.simplify(keep_nodes_F)
    treeseq_YChrom = treeseq_raw.simplify(keep_nodes_M)

    if(treeseq_mtDNA.num_trees != 1): raise ValueError("more than one tree!")
    if(treeseq_YChrom.num_trees != 1): raise ValueError("more than one tree!")

    tree_mtDNA = treeseq_mtDNA.first()
    tree_YChrom = treeseq_YChrom.first()

    #Creating dicts for subpopulations

    malePopulations = {}
    for i in range(treeseq_YChrom.num_populations):
        malePopulations[f"Subpopulation {i+1}"] = [ind for ind in tree_YChrom.leaves(tree_YChrom.root) if tree_YChrom.population(ind) == i]

    femalePopulations = {}
    for i in range(treeseq_mtDNA.num_populations):
        femalePopulations[f"Subpopulation {i+1}"] = [ind for ind in tree_mtDNA.leaves(tree_mtDNA.root) if tree_mtDNA.population(ind) == i]

    #Define Haplogroups

    def mutateEdge(tree: tskit.Tree, node, dataSource: str):
        edgeMutationDict = dict()
        if dataSource == "YChrom":
            L=900000
            m=3.01*10**-8
        elif dataSource == "mtDNA":
            L=16500
            m=6.85*10**-7
        else:
            raise Exception("invalid dataSource argument!")
        return 1-((1-m)**(L*tree.branch_length(node))) > np.random.uniform(0,1)


    def defineHaplogroup(tree: tskit.Tree, dataSource, time):
        haplogroups = {}
        for u in tree.nodes(root=tree.root, order="timedesc"):
            if u == tree.root or tree.time(tree.parent(u)) <= time:
                pass
            else:
                if mutateEdge(tree, u, dataSource):
                    for i in haplogroups.values():
                        for v in tree.leaves(u):
                            i.discard(v)
                    haplogroups[f"H_{u}"] = set([u for u in tree.leaves(u)])
        haplogroups = {key:val for key, val in haplogroups.items() if len(val) != 0}
        return haplogroups


    parameters = re.findall("\d+", input)
    impParameters = [int(parameters[-1]),int(parameters[-3])+int(parameters[-1])]
    time = np.mean(impParameters)

    haplogroups_mtDNA = defineHaplogroup(tree_mtDNA, "mtDNA", time)
    haplogroups_YChrom = defineHaplogroup(tree_YChrom, "YChrom", time)


    #Calculate the Haplogroup Frequencies for every Population

    def calculateHFrequency(dataSource,subpop,*haplogroup):
        subpopHFrequencies = {}
        for j,i in enumerate(haplogroup):
            hgroupCounter = 0
            for m in subpop:
                if m in i:
                    hgroupCounter += 1
            subpopHFrequencies[f'{dataSource}_H{j+1}'] = hgroupCounter/len(subpop)
        return subpopHFrequencies

    yChromDict = {}
    yChromDict["subpop1"] = calculateHFrequency('Y', malePopulations['Subpopulation 1'], *haplogroups_YChrom.values())
    yChromDict["subpop2"] = calculateHFrequency('Y', malePopulations['Subpopulation 2'], *haplogroups_YChrom.values())
    yChromDict["subpop3"] = calculateHFrequency('Y', malePopulations['Subpopulation 3'], *haplogroups_YChrom.values())
    ychromDf = DataFrame(yChromDict)

    mtDnaDict = {}
    mtDnaDict["subpop1"] = calculateHFrequency('M', femalePopulations['Subpopulation 1'], *haplogroups_mtDNA.values())
    mtDnaDict["subpop2"] = calculateHFrequency('M', femalePopulations['Subpopulation 2'], *haplogroups_mtDNA.values())
    mtDnaDict["subpop3"] = calculateHFrequency('M', femalePopulations['Subpopulation 3'], *haplogroups_mtDNA.values())
    mtDnaDf = DataFrame(mtDnaDict)

    Df_untransposed = pd.concat([ychromDf, mtDnaDf])
    Df = Df_untransposed.transpose()

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(Df)

    output = open(outputFile, "w")
    output.write(Df.to_csv())
    output.close()

inputFile = path
outputFile = inputFile.replace(".trees", ".csv")

haploFreqCalc(inputFile)