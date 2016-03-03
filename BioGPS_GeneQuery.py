# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 09:57:39 2016
@author: ellisrj2
"""
import requests
import json 
import itertools as it
import csv 
import pandas

def BioGPS_gene_expression(gene_list, brain_region, species):
    #read list of genes
    gene_list = [gene.rstrip('\n') for gene in open(gene_list)] 
    genes = list(set(gene_list)) # removes duplicates
    
    gene_reads = []
    ids = []
    probes = []
    
    for gene in genes:
        # for each gene, query the database
        gene_query = requests.get('http://mygene.info/v2/query?q=' + gene + '&species=' + species) # Response object
        BGPSdict = json.loads(gene_query.content) # Turns it into a dictionary
        # if the response works (200)
        if len(BGPSdict['hits']) > 0: # if BioGPS found hits for the gene
            gene_reads.append(gene)
            hits = BGPSdict['hits']
            for i in range(len(hits)):
                #if one of the hits matches the queried gene
                if hits[i]['symbol'] == gene:
                    ids.append(hits[i]['_id']) # Accession number
                    id_query = requests.get('http://mygene.info/v2/gene/' + hits[i]['_id']) #Response object
                    id_dict = json.loads(id_query.content) # Turns it into a dictionary
                    if 'reporter' in id_dict:
                        if 'Mouse430_2' in id_dict['reporter']:
                            probes.append(id_dict['reporter']['Mouse430_2']) # Affymetrix probe
                            break
                        else:
                            probes.append('No probe from the Mouse430_2 data set')
                            break
                    else:
                        probes.append('No probes in BioGPS')
                        break
                # If no hits matching exactly by name, take first hit
                elif len(BGPSdict['hits']) > 0:
                    ids.append(hits[0]['_id'])
                    id_query = requests.get('http://mygene.info/v2/gene/' + hits[0]['_id']) # Response object
                    id_dict = json.loads(id_query.content) # Turns it into a dictionary        
                    if 'reporter' in id_dict:
                        if 'Mouse430_2' in id_dict['reporter']:
                            probes.append(id_dict['reporter']['Mouse430_2']) # Affymetrix probe
                            break
                        else:
                            probes.append('No probe from the Mouse430_2 data set')
                            break
                    else:
                        probes.append('No probes in BioGPS')
                        break
                else:
                    ids.append('No hits matched the gene')
                    probes.append('No hits matched the gene')
                    break
        else:
            gene_reads.append(gene)
            ids.append('No hits in BioGPS')
            probes.append('No hits in BioGPS')
    
    # Most genes have multiple probes which are made into a list. However, when a gene has only one probe, that one probe is stored as a string. This loop turns any single probes affiliated with a gene into a list of length 1 which makes it easier to index later in the script.         
    for i in range(len(probes)):
        if type(probes[i]) != list:
            probes[i] = [probes[i]]
    
    row_genes = []
    ids_genes = []
    
    # list the gene once for each probe so that on a spreadsheet, the genes and probes are correctly matched
    for i in range(len(gene_reads)):
        for j in range(len(probes[i])):
            row_genes.append(gene_reads[i])
            ids_genes.append(ids[i])
            
    allprobes = it.chain(*probes)
    allprobes = list(allprobes)
    
    # Read CSV
    df = pandas.read_csv('geneatlas_MOE430_20090327.raw.avg.csv')
    
    df_probes = list(df['Probe']) # read probe list
    df_brain_region = list(df[brain_region]) # read striatum expression values list
    
    expression_vals = []
    for probe in allprobes:
        if probe in df_probes:
            probe_index = df_probes.index(probe)
            expression_vals.append(df_brain_region[probe_index])
        else:
            expression_vals.append('probe not in list')
    
    with open(brain_region + '_expression.csv', 'wb') as thefile:
        writer = csv.writer(thefile)
        writer.writerow(["Gene", "ID", "Probe", brain_region + " Expression"])
        writer.writerows(it.izip_longest(row_genes, ids_genes, allprobes, expression_vals))