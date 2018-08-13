#!/usr/bin/python

"""
Created on Fri Mar 17 12:31:23 2017
@author: francinecamacho
"""

"""This script takes in a BLAST tabular output file, a fasta file used to BLAST against itself to parse out the
 duplicated clusters, proteins or genes from the fasta file. Assumption is that the BLAST tabular file is run with 
 --perc_identity parameter in BLAST. This is necessary for BLAST to calculate the query coverage on hits that passed
 with the desired percent identity. The output is the de-replicated fasta file, and a text file with the unique cluster, 
 proteins or gene names."""

from Bio import SeqIO
import os
import pandas as pd

"Cluster match class hold, the subject name of blast match"
class Cluster(object):
    def __init__(self, name, length):
        self.name = name
        self.length = length
     
     
"""
Function to make a panda data frame to iterate our BLAST tabular file output. """
def makeDataFrame(PATH, tabularHeader, coverageCutOff, perc_identity):
    names = tabularHeader.split(" ")
    dataframe = pd.read_csv(PATH, sep="\t", names=names, header=None)
    # filter dataframe to just have hits that meet qcovs threshold and pident threshold and qseqid are not identical to sseqid
    filteredDF = dataframe[(dataframe.sseqid != dataframe.qseqid) & (dataframe.qcovs >= coverageCutOff) &
                           (dataframe.pident >= perc_identity)]
    # retrieve list of unique queries
    uniqueQID = dataframe['qseqid'].unique()
    print("Number of inputted BGCs:", len(uniqueQID))
    return uniqueQID, filteredDF


"""
Function to create dictionary of BGCs with their respective length and other BGCs that match
at pident and coverage of 95/98 percent or greater
resultDict['V1.CD1.0__scaffold_6__140632__unknownMolecule__arylpolyene__ANTISMASH__107405_140632'] ==
{'matches': [class(O2.UC41.2__scaffold_123__38976__unknownMolecule__arylpolyene__ANTISMASH__1616_38976),
class(V1.UC5.3__scaffold_364__61116__unknownMolecule__arylpolyene__ANTISMASH__0_44388)], 'length': [33227]}"""
def createResultsDict(df):
    resultDict = {}

    # iterate through pandas df structure to get metadata of blast hit
    for index, row in df.iterrows():
        queryName = df.at[index, 'qseqid']
        subjectName = df.at[index, 'sseqid']
        qlength = df.at[index, 'qlen']
        subjectLen = df.at[index, 'slen']

        matchesDict = {}

        if queryName not in resultDict:
            # create matches class (subjects information)
            subjClusterClass = Cluster(subjectName, subjectLen)
            matchesDict.setdefault('qlength', []).append(qlength)
            matchesDict.setdefault('matches', []).append(subjClusterClass)
            resultDict[queryName] = matchesDict


        else:  # queryname is a key in the resultDict
            # if we have another hit at 98 pident or more
            # check if subjectid is already a match in the query['matches']
            clusterIndex = indexOfCluster(subjectName, resultDict[queryName]["matches"])
            if clusterIndex != -1:  # we already have a class with the subjName we don't need to do anything else
                continue
            else:
                subjClusterClass = Cluster(subjectName, subjectLen)
                resultDict[queryName]['matches'].append(subjClusterClass)

    return resultDict


# Function to determine is subjectid is already reported in queryid['matches'].
def indexOfCluster(name, clusterArray):
    i = 0

    for cluster in clusterArray:
        if cluster.name == name:
            return i
        i += 1

    return -1


# Function to find the longest match for a query in terms of sequences length
def findMaxSubject(matchesArray):
    maxSubject = matchesArray[0]
    for match in matchesArray:
        if maxSubject.length < match.length:
            maxSubject = match
    return maxSubject


# Function to determine uniqueness and return a list of unique proteins, genes
def findUniqueClusters(filteredDict, uniqueIDArray, outdir, outfile):
    os.chdir(outdir)
    notUniqueKeys = []  # stores not unique items

    for query in filteredDict:

        if query not in notUniqueKeys:
            # matches of query
            clusterMatches = filteredDict[query]['matches']
            # check to determine if there is a match or not per query
            if len(clusterMatches) != 0:
                # find the match that has the greatest sequence length out
                # of all
                maxMatch = findMaxSubject(clusterMatches)
                # iterate the matches to check the length
                for subject in clusterMatches:
                    # check if query < subject
                    if filteredDict[query]['qlength'][0] > subject.length:
                        if subject.name not in notUniqueKeys:
                            notUniqueKeys.append(subject.name)  # check if the key is in the notunique list
                    # if subject > query
                    elif filteredDict[query]['qlength'][0] < subject.length:
                        # check if the key is in the notunique list
                        if query not in notUniqueKeys:
                            notUniqueKeys.append(query)
                        # need to check that the subject is no the max Match and not already in the list
                        if subject.name != maxMatch.name and subject.name not in notUniqueKeys:
                            notUniqueKeys.append(subject.name)
                    else:  # they are equal
                        if subject.name not in notUniqueKeys:
                            notUniqueKeys.append(subject.name)
            else:
                continue
        else:
            continue

    # Find all unique keys
    uniqueKeys = set(uniqueIDArray) - set(notUniqueKeys)
    uniqueKeys_list = list(uniqueKeys)
    print("Number of unique BGCs:", len(uniqueKeys_list))
    if len(uniqueKeys_list) == 0:
        print("Error: Could not identify any duplicates from input file")
    else:
        uniqueKeys_series = pd.Series(uniqueKeys_list)
        uniqueKeys_DF = uniqueKeys_series.to_frame("bgcName")
        bgcNameFile = outfile + "_uniqueBGCNames.txt"
        uniqueKeys_DF.to_csv(bgcNameFile, index=False,
                             sep='\t')  # write dataframe to csv format (text file) of unique BGCs name

    return uniqueKeys_list  # return unique cluster list to map list to our BGC master list to make a fasta file


"""Function to create a fasta file using the list of uniwue BGCs and map those to the original master bgc fasta file """


def createFastaFile(uniqueBGCList, bgcMasterList, outdir, outfile):
    # uniqueBGC list
    os.chdir(outdir)
    fastafileName = outfile + "_uniqueBGCs.fa"
    with open(fastafileName, 'w') as uniqueFile:
        for seq_record in SeqIO.parse(bgcMasterList, 'fasta'):
            if seq_record.id in uniqueBGCList:
                seq_record.description = ""
                SeqIO.write(seq_record, uniqueFile, "fasta")
            else:
                continue
    uniqueFile.close()


def main(tabular_dir, tabular_file, outdir, outfile, perc_identity, coverage_cutoff, bgc_master_file,
         tabular_file_header):
    os.chdir(tabular_dir)
    tabularFilePath = os.path.join(tabular_dir, tabular_file)

    statinfo = os.stat(tabularFilePath)
    if statinfo.st_size != 0:  # if tabular file is not an empty continue else skip

        uniqueQIDArray, dfObject = makeDataFrame(tabularFilePath, tabular_file_header, coverage_cutoff, perc_identity)

        matchesDict = createResultsDict(dfObject)
        uniqueKeys_list = findUniqueClusters(matchesDict, uniqueQIDArray, outdir, outfile)
        createFastaFile(uniqueKeys_list, bgc_master_file, outdir, outfile)

    else:
        print("ERROR: Inputted tabular file is empty.")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--tabular_dir', type=str, required=True, help='directory where tabular file input is located')
    parser.add_argument('--tabular_file', type=str, required=True, help='blast tabular file results')
    parser.add_argument('--outdir', type=str,
                        help='directory to output fasta file unique BGCs and text file of unique BGCs name')
    parser.add_argument('--outfile', type=str, required=True, help='name of cohort or project')
    parser.add_argument('--perc_identity', type=int, required=False)
    parser.add_argument('--coverage_cutoff', type=int, required=True)
    parser.add_argument('--bgc_master_file', type=str, required=True,
                        help='fasta file used for BLAST and de-replication')
    parser.add_argument('--tabular_file_header', type=str, required=False,
                        default="sseqid qseqid slen qlen qcovs pident Evalue qstart qend")

    args = parser.parse_args()

    main(args.tabular_dir, args.tabular_file, args.outdir, args.outfile, args.perc_identity, args.coverage_cutoff,
         args.bgc_master_file, args.tabular_file_header)
