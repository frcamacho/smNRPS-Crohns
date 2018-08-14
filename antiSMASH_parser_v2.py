#!/bin/env python
"""
Created on Wed Jun  13 12:56:23 2016
@author: francinecamacho
"""
from Bio import SeqIO
import os 


def parseGenbankFile(gbkExtension, sampleID, gbkFilePath, outputFile, outputDir):
	os.chdir(outputDir)
	with open(outputFile, 'w') as parsedFile:
		
		for fileName in os.listdir(gbkFilePath):

			os.chdir(gbkFilePath)
			if fileName.endswith(gbkExtension): #antiSMASH output file 

				for seq_record in SeqIO.parse(fileName, 'genbank'): # iterate through features in gb file
					scaffold_id = seq_record.description
					features = seq_record.features

					#initialization flags 
					inCluster = False 
					clusterLength = None 
					clusterType = None 
					biosytheticStart = []
					biosytheticEnd = [] 
					featureCount = 0 
					for gb_feature in features: # iterating through the features in SeqIO structure 
						if gb_feature.type == 'cluster': # if the feature is a cluster 
							clusterType = ''.join(gb_feature.qualifiers.get('product'))
							inCluster = True  
							clusterLength = gb_feature.location.nofuzzy_end
#                       clusternum = gb_feature.qualifiers.get("note")[0].replace("Cluster number: ", "")

						if gb_feature.type == 'CDS' and inCluster:
							
							biosytheticStart.append(gb_feature.location.nofuzzy_start)
							biosytheticEnd.append(gb_feature.location.nofuzzy_end)  
						checkNextFeature = featureCount+1
						printRow = False 
						try:
								if features[checkNextFeature].type == "cluster" or features[checkNextFeature].type == "source":
									printRow = True
									
						except IndexError:
								continue 
						if printRow == True and inCluster == True:
							clusterStart = min(biosytheticStart)
							clusterEnd = max (biosytheticEnd)
							clusterSeq = str((seq_record.seq))[clusterStart:clusterEnd] 

							#clusterStart here add 1 for the originial cluster start since in the gbk file starts at 1. 
							row = ">"+ sampleID+'__'+scaffold_id+ '__'+ str(clusterLength) +'__'+ clusterType +'__'+'ANTISMASH'+'__'+ str(clusterStart+1)+'_'+ str(clusterEnd) + '\n'+ clusterSeq+'\n'
							os.chdir(outputDir)
							parsedFile.write(row)
							biosytheticStart = []
							biosytheticEnd = [] 
						featureCount+=1 
						
						

			else: 
				continue

		
	parsedFile.close()

def main(sampleID, gbk_path, outdir, outfile):

	gbkExtension = ".final.gbk"
	parseGenbankFile(gbkExtension, sampleID, gbk_path, outfile, outdir)


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--sampleID', type=str)
	parser.add_argument('--gbk_path', type=str)
	parser.add_argument('--outdir', type= str) 
	parser.add_argument('--outfile', type= str) 


	args = parser.parse_args()


	main(args.sampleID, args.gbk_path, args.outdir, args.outfile)
