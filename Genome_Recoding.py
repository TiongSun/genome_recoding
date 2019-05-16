'''
@Author: TS 27/11/2017
 
This script can automatically:
    1. identify all targeted codons & overlapping genes;
    2. report any embedded genes;
    3. split only overlapping genes where recoding will introduce mutation(s);
    4. identify the orientation of upstream and downstream overlapping genes;
    5. add synthetic insertion of customized size (20bp by default) in all split genes except when upstream gene is forward and downstream gene is reverse;
    6. re-annotate positions of all features in the genbank file after split;
    7. annotate synthetic insertion;
    8. recode the genome and annoate all recoded codons;
    9. export a final genbank file as output;
    10. print the total processing time.  

Usage: Just fill in the input session and click run. That is it!    

'''
#==============================================================================
# Input
#==============================================================================
Input_GenBank_File = 'mds AP012306 2.gb'
Output_GenBank_File = 'recoded_genome.gb'

recoding_scheme = [{'target': 'TAG', 'recode': 'TAA'},
                   {'target': 'TCG', 'recode': 'AGC'},
                   {'target': 'TCA', 'recode': 'AGT'}]

synthetic_insertion_size= 20

#==============================================================================
# Import Packages
#==============================================================================
from Bio import SeqIO
from Bio.Seq import Seq
import itertools 
import time


#==============================================================================
# Define functions 
#==============================================================================
# get a list of codon properties: 0 = codon, 1 = codon start, 2 = codon end, 3 = codon position set; 4 = strand ('f' or 'r'), 5= new genbank annotation
def get_codon(genbank, target_codon):  
    codon_pos_list = []
    if type(target_codon) != list:
        target_codon = [target_codon] 
    for feature in genbank.features:
        if feature.type == "CDS":
            start=feature.location.start.position
            end=feature.location.end.position
            sense=feature.strand   
            if sense == 1:                                                  
                orf = genbank.seq[start:end]   
                codon = [orf[i:i+3] for i in range(0,len(orf),3)]              # read codon in triplet
                codon_pos = [{'codon': str(x),
                              'start': (i*3)+start,
                              'end': (i*3)+start+3,
                              'position': set(range((i*3)+start, (i*3)+start+3)), # set(range(start, end))
                              'strand': 'f'} for i, x in enumerate(codon) if x in target_codon] # extract codon info       
            else:
                orf = genbank.seq[start:end].reverse_complement()              # if antisense, do reverse transcription
                codon = [orf[i:i+3] for i in range(0,len(orf),3)]                          
                codon_pos = [{'codon': str(x),
                              'start': ((len(codon)-i)*3)+start,
                              'end': ((len(codon)-i)*3)+start-3,
                              'position': set(range(((len(codon)-i)*3)+start-3, ((len(codon)-i)*3)+start)), # set(range(start, end))
                              'strand': 'r'} for i, x in enumerate(codon) if x in target_codon] 
            codon_pos_list.append(codon_pos)    
    return list(itertools.chain.from_iterable(codon_pos_list))

# get CDS propertie: 0 = index; 1 = gene; 2= start; 3 = end; 4 = strand
def get_CDS(genbank):
    count = 0
    CDS_info = []
    for feature in genbank.features:
        if feature.type == "CDS":
            start = feature.location.start.position
            end = feature.location.end.position
            sense = feature.strand
            strand = ['f' if x == 1 else 'r' for x in [sense]][0]
            try: 
                gene = feature.qualifiers["gene"][0]
            except KeyError:
                gene = 'NA'
            if sense == 1:
                aa = genbank.seq[start:end].translate()
            else:
                aa = genbank.seq[start:end].reverse_complement().translate()              
            CDS_info.append({'index': count, 'gene': gene, 'start': start, 'end': end,
                             'position': set(range(start, end)), 'strand': strand, 'aa': aa})
            count += 1
    return CDS_info

def report_embedded_genes(CDS): # CDS = output from get_CDS
    for i in range(1,len(CDS)):
        if len(CDS[i].get('position') & CDS[i-1].get('position')) == len(CDS[i].get('position')):
            print('Warning: {} at position {} is an embedded gene. '.format(CDS[i].get('gene'),CDS[i].get('start')) + 
                      'Please consider how this gene should be split')
            
# get overlapping info: 0 = index; 1 = overlapping codon; 2= upstream gene; 3 = downstream gene; 4 = start position; 5 = end position; 6 = orientation
def get_overlap(codon_list): # codon_list = output from get_codon_info
    overlap_list = []
    for i in range(1,len(CDS)):  
        upstream = CDS[i-1]
        downstream = CDS[i]
        upstream_set = set(range(upstream.get('start'), upstream.get('end')+1))
        downstream_set = set(range(downstream.get('start'), downstream.get('end')+1))
        overlap_region = (upstream_set & downstream_set)
        if bool(overlap_region) == True:
            overlap_codon = [codon_list[i].get('start') for i in range(len(codon_list)) if bool(codon_list[i].get('position') & overlap_region)==True]
            if bool(overlap_codon) == True:
                overlap_list.append({'index': i,
                                'codon': overlap_codon,
                                'upstream gene': upstream.get('gene'),
                                'downstream gene': downstream.get('gene'),
                                'start': min(list(overlap_region)),
                                'end': max(list(overlap_region)),
                                'orientation': upstream.get('strand') + downstream.get('strand')})
    return overlap_list
                
def recoding(seq, codon_list): 
    mutable_seq = seq.tomutable()
    for i in codon_list:
        start = i.get('start')
        end = i.get('end')
        codon = i.get('codon')
        if i.get('strand') == 'f':
            mutable_seq[start: end] = [a.get('recode') for a in recoding_scheme if a.get('target') == codon][0]
        else:
            mutable_seq[end: start] = [str(Seq(a.get('recode')).reverse_complement()) for a in recoding_scheme if a.get('target') == codon][0]
    return mutable_seq.toseq()


 
if __name__ == '__main__':
    start_time = time.time()
    #==============================================================================
    # Get essential information of the genbank file
    #==============================================================================
    mds = SeqIO.read(Input_GenBank_File, "genbank")                                # Read genbank files    
    target_codon = [i.get('target') for i in recoding_scheme]                      # get a list of target codon
    codon_list = get_codon(mds, target_codon)                                      # get info on the target codons
    CDS = get_CDS(mds)                                                             # get info on all CDS in the genbank file
    report_embedded_genes(CDS)                                                     # check for embedded gene(s)
    overlap = get_overlap(codon_list)                                              # get info on all the overlapping CDS in the genbank
    
    #==============================================================================
    # Identify overlaps that will introduce mutation after recoding
    #==============================================================================
    test_recoding = SeqIO.read(Input_GenBank_File, 'genbank')                      # create a genbank file to test whether recoding will introduce mutation at overlap regions
    test_recoding.seq = recoding(test_recoding.seq, codon_list)# recoding
    CDS_test_recoding = get_CDS(test_recoding) 
    
    overlap_mutated = []
    for i in range(len(overlap)):
        index = overlap[i].get('index')
        if CDS_test_recoding[index-1].get('aa') != CDS[index-1].get('aa') or CDS_test_recoding[index].get('aa') != CDS[index].get('aa'):
                overlap_mutated.append(index)
            
    #==============================================================================
    # Split overlapping genes
    #==============================================================================
    from Bio.SeqFeature import FeatureLocation
    from Bio.SeqFeature import SeqFeature
    
    split = mds.seq
    for i in range(len(overlap)-1, -1,-1):                                         # read overlap in reverse direction
        if overlap[i].get('index') in overlap_mutated:                             # overlap need to be fix or else it will introduce missense mutation
            start = overlap[i].get('start') 
            end = overlap[i].get('end') 
            orientation = overlap[i].get('orientation')
            length = end - start 
            if orientation == 'fr':                                                # no synthetic insertion
                split = (split[:end] + split[start:])
            else:
                split = (split[:end] + split[start - synthetic_insertion_size: start] + split[start:])           
            for feature in mds.features:                                           # After split, all subsequent genbank annotated has to be adjusted
                if feature.location.start.position >= start:
                    if orientation == 'fr':                                        # no synthetic insertion
                        new_start = int(feature.location.start.position + length)
                        new_end = int(feature.location.end.position + length)
                    else:
                        new_start = int(feature.location.start.position + length + synthetic_insertion_size)
                        new_end = int(feature.location.end.position + length + synthetic_insertion_size)
                    feature.location = FeatureLocation(new_start, new_end, feature.location.strand)
            if orientation != 'fr':                                                # annotate synthetic insertion
                mds.features.append(SeqFeature(FeatureLocation(end, end + synthetic_insertion_size), 
                                           qualifiers = {"note": 'Synthetic Insertion'}, type = "misc_feature"))
    
    mds.seq = split           
    
    #==============================================================================
    # Recode and annotate the genome
    #==============================================================================
    codon_list_after_split = get_codon(mds, target_codon)  
    mds.seq = recoding(mds.seq, codon_list_after_split)                            # recoding
        
    for i in codon_list_after_split:
        recoded_codon = [recoding_scheme[a].get('recode') for a in range(len(recoding_scheme)) if i.get('codon') == recoding_scheme[a].get('target')][0]  # get recoded codon from corresponding target codon
        annotation = (i.get('codon') + ' to ' + recoded_codon)
        if i.get('strand') == 'f':
            mds.features.append(SeqFeature(FeatureLocation(i.get('start'), i.get('end')), 
                                           qualifiers = {"note": annotation}, type = "misc_feature"))
        else: 
            mds.features.append(SeqFeature(FeatureLocation(i.get('end'), i.get('start')), 
                                           qualifiers = {"note": annotation}, type = "misc_feature"))
        
    #==============================================================================
    # Export recoded genbank file and calculate processing time
    #==============================================================================                                                # convert mutable sequence back into non-mutable SeqObject                                                      # replace the DNA sequence in the original genbank file with the new recoded DNA sequence
    SeqIO.write(mds, Output_GenBank_File, "genbank")                         
    
    end_time = time.time() - start_time
    print('Total processing time = {}'.format(time.strftime('%H:%M:%S', time.gmtime(end_time))))
    




