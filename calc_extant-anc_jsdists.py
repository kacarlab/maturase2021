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

def _calcdistn(seqs):
  aacodes = '-ARNDCEQGHILKMFPSTWYV'
  distns = []
  # iterate over alignment columns
  for i in range(len(seqs[0])):
    col = [1.0]*len(aacodes)
    for seq in seqs:
      idx = aacodes.find(seq[i])
      if idx > -1:
        col[idx] += 1.0
    # normalize distribution
    tot  = sum(col)
    dist = [x/tot for x in col]
    # append column distribution
    distns.append(dist)
  # return all column distributions
  return distns


def _shannon(distns1, distns2):
  dists = []
  # iterate over alignment columns
  for i in range(len(distns1)):
    d1 = distns1[i]
    d2 = distns2[i]
    dist = scipy.spatial.distance.jensenshannon(d1,d2)
    dists.append(dist)
  return dists


def run(fname):
  # read sequence data
  seqdict = _parsefasta(fname)
  d_patt  = 'D_'
  e_patt  = 'E_'
  d_e_distns = [[],[]]
  # divide data into matching/nonmatching
  dgrp = []
  egrp = []
  for k in seqdict.keys():
    if k.find(d_patt) == 0:
      dgrp.append(seqdict[k])
    elif k.find(e_patt) == 0:
      egrp.append(seqdict[k])
    # calculate distributions
  d_e_distns = [_calcdistn(dgrp), _calcdistn(egrp)]

  # read ancestral distributions
  ancfname = fname.split('.')[0] + '.asr.probdists.csv'
  node_col_distns = {}
  with open(ancfname, 'r') as handle:
    handle.readline()
    for line in handle:
      linearr = line.strip().split(',')
      node = linearr[0]
      col  = int(linearr[1])
      gap  = float(linearr[2])
      if linearr[-1] == '*':
        aas = [ float(x)*(1.0-gap) for x in linearr[3:-1] ]
      else:
        aas = [ float(x)*(1.0-gap) for x in linearr[3:] ]
      cold = [gap]
      cold.extend(aas)
      if node not in node_col_distns.keys():
        node_col_distns[node] = {}
      node_col_distns[node][col] = cold

  outfname = fname.split('.fasta')[0] + '.shannondists.csv'
  with open(outfname, 'w') as outf:
    outf.write('node,col,D_dist,E_dist,D-E\n')
    # get column-wise shannon distances
    for node in node_col_distns.keys():
      # fill in column distributions for this node
      n_cols = list(node_col_distns[node].keys())
      n_cols.sort()
      n_distns = []
      for i in n_cols:
        n_distns.append(node_col_distns[node][i])
      # calculate shannon distances
      d_dists = _shannon(n_distns, d_e_distns[0])
      e_dists = _shannon(n_distns, d_e_distns[1])
      # print info to output file
      for j in range(len(d_dists)):
        outf.write('{},{},{:.4f},{:.4f},{:0.6f}\n'.format(node,
                                                  j+1,
                                                  d_dists[j],
                                                  e_dists[j],
                                                  d_dists[j]-e_dists[j]))


fnamelist = list(glob.glob('*.extanc.fasta'))
with multiprocessing.Pool(len(fnamelist)) as pool:
  pool.map(run, fnamelist, chunksize=1)
