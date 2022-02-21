#!/usr/bin/env python3

import glob
import random
import multiprocessing
import numpy
import scipy.stats
import scipy.spatial.distance

def _parsefasta(fname):
  data = {}
  with open(fname, 'r') as handle:
    line = handle.readline()
    while line:
      if line[0] == '>':
        id = line[1:].strip()
        se = ''
        line = handle.readline()
        while line and line[0] != '>':
          se += line.strip()
          line = handle.readline()
        data[id] = se
  return data

def _shannon(seqs1, seqs2):
  aacodes = '-ARNDCEQGHILKMFPSTWYV'
  dists = []
  # iterate over alignment columns
  for i in range(len(seqs1[0])):
    col1 = [1.0]*len(aacodes)
    col2 = [1.0]*len(aacodes)
    for seq in seqs1:
      idx = aacodes.find(seq[i])
      if idx > -1:
        col1[idx] += 1.0
    for seq in seqs2:
      idx = aacodes.find(seq[i])
      if idx > -1:
        col2[idx] += 1.0
    dist = scipy.spatial.distance.jensenshannon(col1,col2)
    dists.append(dist)
  return dists

def _bootstrap(seqs, n):
  bootstraps = []
  nboots = 100
  for i in range(nboots):
    # shuffle seqs
    random.shuffle(seqs)
    # split into ingroup and outgroup
    ing  = seqs[:n]
    outg = seqs[n+1:]
    # calculate distances
    bootstraps.append(_shannon(ing,outg))
  return bootstraps

def _calcp(observed, distn, index):
  null = []
  for x in distn:
    null.append(x[index])
  if numpy.std(null) < 1.0e-10 or observed < 1.0e-10:
    return 1.0
  kde = scipy.stats.gaussian_kde(null)
  pvalue = 1.0 - kde.integrate_box_1d(0.0,observed)
  #null.sort()
  #i = 0
  #while i < len(null) and observed > null[i]:
  #  i += 1
  #pvalue = 1.0 - (i/len(null))
  #if pvalue < 1.0 / len(null):
  #  pvalue = 1.0 / len(null)
  return pvalue


DPATT = 'D_'
EPATT = 'E_'

def run(fname):
  # read sequence data
  seqdict = _parsefasta(fname)
  for patt in [DPATT, EPATT]:
    # divide data into matching/nonmatching
    ingrp  = []
    outgrp = []
    allgrp = []
    for k in seqdict.keys():
      if k.find(DPATT) == 0 or k.find(EPATT) == 0:
        allgrp.append(seqdict[k])
        if k.find(patt) == 0:
          ingrp.append(seqdict[k])
        else:
          outgrp.append(seqdict[k])
    # calculate shannon distances
    dists = _shannon(ingrp, outgrp)
    bootdists = _bootstrap(allgrp, len(ingrp))
    # print info to output file
    outfname = patt + fname.split('.fasta')[0] + '.extantde_jsdists.csv'
    with open(outfname, 'w') as outf:
      outf.write('col,distance,pvalue\n')
      for i in range(len(dists)):
        pvalue = _calcp(dists[i],bootdists,i)
        outf.write('{},{:.4f},{:.4e}\n'.format(i+1,dists[i],pvalue))


fnamelist = list(glob.glob('*.fasta'))
with multiprocessing.Pool(len(fnamelist)) as pool:
  pool.map(run, fnamelist, chunksize=1)

