#!/usr/bin/env python3
'''Truncate Sequences to primer locations'''
def truncate2primer(file, primerlist= '/home/boom/amp/primers/unambiguous_primers_2020.02.fna', primerblastdb = '/home/boom/amp/primers/unambiguous_primers_2020.02', startprimers = ['515F_mod'], endprimers = ['806R_Mod','806F_Mod'], min_pct = 90, plusleft = 0, plusright = 0, min_len_pct = 90,sides = 'both',del=True):
    from holotools import biop
    import os
    import shutil
    import pandas as pd

    cwd = os.getcwd()
    # temp file set up for data dumps and manipulation
    try:
        os.mkdir('tmp')
    except:
        print('tmp directory already exists, will overwrite')
        shutil.rmtree('tmp')
        os.mkdir('tmp')
    # read in query file
    d = biop.fdict(file)
    p = biop.fdict(primerlist)
    fdf = pd.DataFrame()
    na = open('nopfound.tsv','w')
    na.write('query\tissue\tdetails\n')
    trunc = {}
    # make tmp files of each seq blast and parse
    for k,v in d.items():
        out = open('%s/tmp/%s.fna'%(cwd,k),'w')
        out.write('>%s\n%s\n'%(k,v))
        out.close()
        os.system('blastn -db %s -query %s -task "blastn-short" -out %s -outfmt 6'%(primerblastdb,cwd+'/tmp/'+k+'.fna',cwd+'/tmp/'+k+'.tsv'))
        try:
            df = pd.read_csv(cwd+'/tmp/'+k+'.tsv', sep = '\t',header = None)
            df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
            df = df.loc[df.pident>=min_pct]
            # end primers
            s = df[df.sseqid.isin(endprimers)]
            for i in endprimers:
                plen = len(p[i])
                ep = s.loc[s.sseqid==i]
                ep = ep.loc[ep.length>=plen*(min_len_pct/100)]
                if len(ep)<1:
                    na.write(k+'\tprimer not found\t%s\n'%i)
                else:
                    mmin = min(ep.qstart)
                fdf = pd.concat([fdf,ep])
            # start primers
            s = df[df.sseqid.isin(startprimers)]
            for i in startprimers:
                plen = len(p[i])
                ep = s.loc[s.sseqid==i]
                ep = ep.loc[ep.length>=plen*(min_len_pct/100)]
                if len(ep)<1:
                    na.write(k+'\tprimer not found\t%s\n'%i)
                if side = 'right':
                    mmax = max(ep.qend)
                fdf = pd.concat([fdf,ep])
        except:
            na.write(k+'\tprimer not found\tall\n')
    if del == True:
        shutil.rmtree('tmp')
    if side == 'right':
        trunc[k]=v[:mmax]
    elif side == 'left':
        trunc[k]=v[mmin:]
    if side == 'both':
        trunc[k]=v[mmin:mmax]

    na.close()
    return fdf, trunc

# truncate2primer('test16.fna')
