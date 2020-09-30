#!/usr/bin/env python3
'''Binders and Tools for Clustering of Bacterial Sequences'''
def blastprimer(query, db='/home/boom/amp/primers/unambiguous_primers_2020.02'):
    import os
    cwd = os.getcwd()
    os.system('blastn -db %s -query %s -task "blastn-short" -out %s -outfmt 6'%(db,query,query+'.tsv'))

def blast16(query, db='/home/boom/Documents/holoweb/static/fasta/dbs/2020.02.C'):
    import os
    cwd = os.getcwd()
    os.system('blastn -db %s -query %s -out %s -outfmt 6'%(db,query,query+'.tsv'))

def parseblast(file):
    '''tsv is expected'''
    import pandas as pd
    df = pd.read_csv(file,sep = '\t',header = None)
    df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    df.sort_values(by='pident').reset_index()
    return df

def addsptofasta(fasta):
    '''take a multifasta, break it appart, blast each against 16s, add the species/pctid to the name, output new multifasta'''
    import os
    import shutil
    import holotools.biop as biop
    import sys
    cwd = os.getcwd()
    spfasta = open('spblast_'+fasta,'w')
    try:
        os.mkdir('tmp')
    except:
        print('tmp directory already exists, will overwrite')
        shutil.rmtree('tmp')
        os.mkdir('tmp')
    d = biop.fdict(fasta)
    n = len(d)
    c = 0
    for k,v in d.items():
        sys.stdout.write('\r')
        c+=1
        j = (c + 1) / n
        sys.stdout.write("[%-20s] %d%%" % ('~'*int(20*j), 100*j))
        sys.stdout.flush()
        out = open(os.getcwd()+'/tmp/'+k,'w')
        out.write('>%s\n%s\n'%(str(k),str(v)))
        out.close()
        try:
            blast16(cwd+'/tmp/'+k)
            df = parseblast(os.getcwd()+'/tmp/'+k+'.tsv')
            sp = df.sseqid.iloc[0]
            pid = df.pident.iloc[0]
            spfasta.write('>%s\n%s\n'%(k+'_'+sp+'_'+str(pid),str(v)))
        except:
            print('There was an issue with %s'%k)
    spfasta.close()
    shutil.rmtree('tmp')
