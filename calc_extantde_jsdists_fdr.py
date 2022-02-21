#!/usr/bin/env python3

import glob
import statsmodels.stats.multitest

for fname in glob.iglob('*.extantde_jsdists.csv'):
  outfname = fname
  header = ''
  lines = []
  pvals = []
  with open(fname,'r') as handle:
    header = handle.readline().strip() + ',qvalue\n'
    for line in handle:
      pvals.append(float(line.split(',')[-1]))
      lines.append(line.strip())

  rejected,qvals = statsmodels.stats.multitest.fdrcorrection(pvals)

  with open(outfname, 'w') as outf:
    outf.write(header)
    i = 0
    for line in lines:
      outf.write(line)
      outf.write(',{:.4e}\n'.format(qvals[i]))
      i += 1
