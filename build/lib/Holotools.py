#!/usr/bin/env python3
'''Binders and Tools for Clustering of Bacterial Sequences'''
from ht_biop import *
from ht_trim import *


d = fdict('test_16.fna')
td = {}
for k,v in d.items():
    print(k)
    td[k]=sliding_window(v)
print(td)
