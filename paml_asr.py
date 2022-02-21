#!/usr/bin/env python3

import sys
import os
import shutil
import re
from optparse import OptionParser, OptionGroup

VERSION = "asr v1.0"

# set default executables and model #
PAMLBS = "/path/to/paml" #Set path to PAML
CODEML = PAMLBS + "/bin/codeml"
BASEML = PAMLBS + "/bin/baseml"
SUBMOD_NAME = PAMLBS + "/dat/lg.dat"


### global variables ###
AMINO_ACIDS = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
GAPS = ["-", "X"]

### default options       ###
FIX_BLENS     = 1           # optimize branch lengths, starting from user-supplied tree
FIX_ALPHA     = 0           # optimize gamma shape parameter
INITIAL_ALPHA = 1.5         # starting gamma shape parameter


### BEG HELPER FUNCTIONS #######################################################
#

### normalize a probability distribution so it sums to 1.0                   ###
def normalize(pdist):
    thedist = []
    for p in pdist:
        if p < 0.0:
            thedist.append(0.0)
        else:
            thedist.append(p)
    total = sum(thedist)
    newdist = [x/total for x in thedist]
    return newdist

### converts an alignment from residues and gaps to W/S coding for paml      ###
### (binary=False) or to 0/1 coding for raxml (binary=True)                  ###
def convertAlignment(alnfname, gapalnfname, binary=False):
    GAP = "W"
    RES = "S"

    if binary:
        GAP = "0"
        RES = "1"

    fastaid = ">"
    realgap = "-"

    handle = open(alnfname, "r")
    outf   = open(gapalnfname, "w")
    for line in handle:
        if line[0] == fastaid:
            outf.write(line)
        else:
            newseq = ""
            for c in line.strip():
                if c == realgap:
                    newseq += GAP
                else:
                    newseq += RES
            outf.write("%s\n" % newseq)
    handle.close()
    outf.close()
    return

### combine a tree with branch lengths and a tree with node labels into      ###
### a single tree with both                                                  ###
def combineTrees(bltree, labeltree):
    fulltree = ""

    bltreearr    = bltree[:-1].split(")")
    labeltreearr = labeltree[:-1].split(")")

    for i in range(len(bltreearr)):
        fulltree += bltreearr[i]
        if i+1 < len(labeltreearr):
            fulltree += ")" + labeltreearr[i+1].split(",")[0]

    return fulltree + ";"

### parse ancestral state probability distribution from paml output          ###
def parseASRprobs(handle, warn=True):
    # N x M array, where N is the sequence length, and  #
    # M=20, one for each amino acid residue. This will  #
    # hold the probability distributions over each site #
    # in the sequence alignment                         #
    site_res_prob_dists = []

    # assume we just read "Prob distribution at node ..."  #
    # so we need to skip over a few lines int the rst file #
    handle.readline()
    handle.readline()
    handle.readline()
    line = handle.readline().strip()
    while line:
        aastrarr = line.split()[3:]
        aadist = []
        i = 0
        for s in aastrarr:
            aa = s[0]
            pp = float(s.split("(")[1][:-1])
            if warn and aa != AMINO_ACIDS[i]:
                sys.stderr.write("WARNING: amino-acid names don't appear to be in the right order!\n")
            aadist.append(pp)
            i += 1
        site_res_prob_dists.append(aadist)
        line = handle.readline().strip()

    return site_res_prob_dists

### parse ML ancestral sequences from paml output                            ###
def parseASRseqs(handle):
    mlseqs = {}

    # assume we just read "List of extant and ...", so #
    # we need to skip over a few lines in the rst file #
    handle.readline()
    handle.readline()
    handle.readline()
    line = handle.readline().strip()
    while line:
        linearr = line.split()
        if linearr[0] == "node":
            nodeid = linearr[1][1:]
            seq    = "".join(linearr[2:])
            mlseqs[nodeid] = seq
        line = handle.readline().strip()

    return mlseqs

### parse the entire paml rst file                                           ###
def parseRST(warn=True):
    bltree         = ""
    labeltree      = ""
    asr_prob_dists = {}
    ml_asrseqs = ""

    nodepattern = re.compile("Prob distribution at node (.+), by site")

    handle = open("rst", "r")
    for line in handle:
        line = line.strip()
        if line.find("Ancestral reconstruction by") == 0:
            handle.readline()
            bltree = handle.readline().strip().replace(" ", "")
        elif line == "tree with node labels for Rod Page's TreeView":
            labeltree = handle.readline().strip().replace(" ", "")
        elif line == "List of extant and reconstructed sequences":
            ml_asrseqs = parseASRseqs(handle)
        else:
            match = nodepattern.search(line)
            if match:
                nodeid    = match.group(1)
                probdists = parseASRprobs(handle,warn)
                asr_prob_dists[nodeid] = probdists
    handle.close()

    # convert separate node-labelled tree and branch-length tree #
    # into a single Newick tree that has both node labels and    #
    # branch lengths                                             #
    complete_tree = combineTrees(bltree, labeltree)

    return(complete_tree, ml_asrseqs, asr_prob_dists)

#
### END HELPER FUNCTIONS #######################################################



### BEG PAML ###################################################################
#

### amino-acid control file template for codeml                              ###
AACTL = """
      seqfile = %s
     treefile = %s
      outfile = %s

        noisy = 3     * 0,1,2,3,9: how much rubbish on the screen
      verbose = 2     * 1: detailed output, 0: concise output
      runmode = 0     * 0: user tree

      seqtype = 2     * 2: AAs
   aaRatefile = %s    * name of substitution model file
        model = 3     * models for AAs or codon-translated AAs:
                        * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                        * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
        Mgene = 0     * 0: same model for all genes

    fix_alpha = %d    * 0: estimate gamma shape parameter; 1: fix it
        alpha = %.2f  * initial alpha, 0->infinity (constant rate)
       Malpha = 0     * 0: same alpha for all genes
        ncatG = 4     * # of categories in the dG or AdG models of rates

        clock = 0     * 0: no molecular clock
        getSE = 0     * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1     * (1/0): rates (alpha>0) or ancestral states (alpha=0)

  Small_Diff = 0.5e-6 * optimization stopping criterion
    cleandata = 0     * 0: don't remove ambiguous sites
        ndata = 1     * we only have a single alignment
  fix_blength = %d    * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 1     * 0: simultaneous; 1: one branch at a time
"""

### presence-absence control file template for codeml                        ###
GAPCTL = """
      seqfile = %s
     treefile = %s
      outfile = %s

        noisy = 3     * 0,1,2,3,9: how much rubbish on the screen
      verbose = 2     * 1: detailed output, 0: concise output
      runmode = 0     * 0: user tree

        model = 2     * 0: JC69; 2: F81 (includes state frequencies)

        Mgene = 0     * 0: same model for all genes

    fix_alpha = 0     * 0: estimate gamma shape parameter; 1: fix it
        alpha = %.2f  * initial alpha, 0->infinity (constant rate)
       Malpha = 0     * 0: same alpha for all genes
        ncatG = 4     * # of categories in the dG or AdG models of rates

        clock = 0     * 0: no molecular clock
        nhomo = 0     * 0: homogeneous
       getSE =  0     * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1     * (1/0): rates (alpha>0) or ancestral states (alpha=0)

  Small_Diff = 0.5e-6 * optimization stopping criterion
    cleandata = 0     * 0: don't remove ambiguous sites
        ndata = 1     * we only have a single alignment
  fix_blength = 1     * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 1     * 0: simultaneous; 1: one branch at a time
"""

### execute ancestral sequence reconstruction using paml                     ###
def doASR(alnfname, phyfname, verbose, clean):
    # file names #
    pamlctlfname     = "paml.ctl"
    pamloutfname     = "paml.out"
    pamllogfname     = "paml.log"
    pamldeffilenames = ["lnf", "rates", "rst", "rst1", "rub", "2base.t", "in.basemlg"]
    gapalnfname      = alnfname + ".presabsaln"

    codeml = CODEML
    submod = SUBMOD_NAME

    if not os.path.exists(codeml):
        sys.stderr.write("ERROR: can't find codeml executable %s\n" % codeml)
        sys.exit(1)

    if not os.path.exists(submod):
        sys.stderr.write("ERROR: can't find paml's substitution model file %s\n" % submod)
        sys.exit(1)

    if verbose:
        sys.stdout.write("executing ancestral residue inference... ")
        sys.stdout.flush()

    # generate codeml control file #
    ctl = AACTL % (alnfname, phyfname, pamloutfname, submod, FIX_ALPHA, INITIAL_ALPHA, FIX_BLENS)
    handle = open(pamlctlfname, "w")
    handle.write(ctl)
    handle.close()

    # execute codeml to reconstruct ancestral amino-acid sequences #
    cmd = "%s %s > %s" % (codeml,pamlctlfname,pamllogfname)
    os.system(cmd)

    if verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write("parsing ancestral residue information... ")
        sys.stdout.flush()

    # parse codeml output #
    (tree, mlSeqs, orig_asrprobs) = parseRST()

    # normalize all the ASR probability distributions #
    ASRProbs = {}
    for (id,pdists) in orig_asrprobs.items():
        newpdists = []
        for pdist in pdists:
            newpdists.append(normalize(pdist))
        ASRProbs[id] = newpdists

    if not clean:
        cfbase = "codeml."
        if os.path.exists(pamlctlfname):
            shutil.move(pamlctlfname, cfbase+pamlctlfname)
        if os.path.exists(pamloutfname):
            shutil.move(pamloutfname, cfbase+pamloutfname)
        if os.path.exists(pamllogfname):
            shutil.move(pamllogfname, cfbase+pamllogfname)
        for f in pamldeffilenames:
            if os.path.exists(f):
                shutil.move(f, cfbase+f)

    if verbose:
        sys.stdout.write("done.\n")

    ## done with codeml. ancestral amino-acid sequences computed ##
    ## now we need to infer ancestral gap characters; this will  ##
    ## be a bit of a hack!                                       ##

    if verbose:
        sys.stdout.write("executing ancestral indel inference... ")
        sys.stdout.flush()

    # convert amino-acid alignment to presence-absence alignment #
    convertAlignment(alnfname, gapalnfname)

    # generate baseml control file #
    baseml = shutil.which("baseml")
    # note that branch lengths and gamma shape will probably be poor for #
    # the gap model, so we will always estimate them.                    #
    ctl = GAPCTL % (gapalnfname, phyfname, pamloutfname, INITIAL_ALPHA)
    handle = open(pamlctlfname, "w")
    handle.write(ctl)
    handle.close()

    # execute baseml to reconstruct ancestral presence-absence states #
    cmd = "%s %s > %s" % (baseml,pamlctlfname,pamllogfname)
    os.system(cmd)

    if verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write("parsing ancestral indel information... ")
        sys.stdout.flush()

    # parse baseml output #
    (gap_tree, baseml_mlseqs, baseml_asrprobs) = parseRST(warn=False)

    # convert baseml asrprobs to presence-absence asrprobs #
    gap_asrprobs = {}
    for (k,v) in baseml_asrprobs.items():
        id = k
        pdists = []
        for parr in v:
            # baseml reports amino-acids in order: T,C,A,G #
            # gap = T or A; no-gap = C or G                #
            # note that we're defining that gap goes first #
            gdist = normalize([parr[0]+parr[2],parr[1]+parr[3]])
            pdists.append(gdist)
        gap_asrprobs[id] = pdists

    if not clean:
        cfbase = "baseml."
        if os.path.exists(pamlctlfname):
            shutil.move(pamlctlfname, cfbase+pamlctlfname)
        if os.path.exists(pamloutfname):
            shutil.move(pamloutfname, cfbase+pamloutfname)
        if os.path.exists(pamllogfname):
            shutil.move(pamllogfname, cfbase+pamllogfname)
        for f in pamldeffilenames:
            if os.path.exists(f):
                shutil.move(f, cfbase+f)

    if verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write("integrating posterior probability distribtuions... ")
        sys.stdout.flush()

    ## done with baseml. ancestral presenece-absence sequences computed ##
    ## now we need to combine the presence-absence sequences with the   ##
    ## previously-computed amino-acid sequences.                        ##

    # calculate the full 'probability distributions', combining gap and   #
    # sequence distributions. Note that we are not integrating these into #
    # a combined distribuiton. We have separate distributions for GAPs    #
    # and for residues.                                                   #
    complete_probdists = {}
    for id in ASRProbs.keys():
        aaprobdists = ASRProbs[id]
        gpprobdists = gap_asrprobs[id]
        probdists = []
        for i in range(len(aaprobdists)):
            aadist  = aaprobdists[i]
            gpdist  = gpprobdists[i]
            newdist = []
            newdist.append(gpdist[0])
            for aa in aadist:
                newdist.append(aa)
            probdists.append(newdist)
        complete_probdists[id] = probdists

    if verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write("calculating maximum-likelihood ancestral sequences... ")
        sys.stdout.flush()

    # calculate ML ancestral sequences with gaps from the probability #
    # distributions, rather than relying on PAML's ML sequences       #
    # If P(gap) >= 0.5, then we set the ML state to gap. Otherwise,   #
    # we set the state to the amino-acid residue with the highest     #
    # posterior probability. In the case of 'ties', the first         #
    # amino-acid residue (in 'alphabetical' order) is returned as the #
    # ML state                                                        #
    complete_seqs = {}
    for (id,probdists) in complete_probdists.items():
        se = ""
        for pdist in probdists:
            if pdist[0] >= 0.5:
                se += GAPS[0]
            else:
                max_i = 0
                max_p = 0.0
                for i in range(1,len(pdist)):
                    if pdist[i] > max_p:
                        max_i = i
                        max_p = pdist[i]
                se += AMINO_ACIDS[max_i-1]
        complete_seqs[id] = se

    if verbose:
        sys.stdout.write("done.\n")

    ## clean up temporary files ##
    if clean:
        cmd = "rm -f %s %s %s %s" % (pamlctlfname, pamloutfname, pamllogfname, gapalnfname)
        for f in pamldeffilenames:
            cmd += " " + f
        os.system(cmd)

    return (tree, complete_seqs, complete_probdists)

#
### END PAML ###################################################################



### BEG MAIN ###################################################################
#

if __name__ == "__main__":
    ### set up command-line argument parsing ###
    optparser = OptionParser(usage="usage: %prog [options] ALIGNMENT.fasta PHYLOGENY.tre",
                             version=VERSION,
                             description="reconstruct ancestral protein sequences")

    optparser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                         help="print some runtime info to the screen [default: %default]")
    optparser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                         help="no news is (probably) good news")
    optparser.add_option("--clean", action="store_true", dest="clean",
                         help="after running the analysis, remove any intermediary or log files, keeping only the results files [default: %default]")
    optparser.add_option("--noclean", action="store_false", dest="clean",
                         help="don't remove intermediary or log files after analysis")

    # options for helper applications #
    group = OptionGroup(optparser, "Helper Options", "where to find helper programs")
    group.add_option("--codeml", action="store", type="string", dest="codeml",
                     help="set name of codeml execuable program to FILE [default: %default]",
                     metavar="FILE")
    group.add_option("--baseml", action="store", type="string", dest="baseml",
                     help="set name of baseml execuable program to FILE [default: %default]",
                     metavar="FILE")
    optparser.add_option_group(group)

    # options controlling the ancestral reconstruction #
    group = OptionGroup(optparser, "ASR Options", "how to perform the ancestral state reconstruction")
    group.add_option("--model", action="store", type="string", dest="pamlmodel",
                     help="set codeml's substitution model to FILE [default: %default]",
                     metavar="FILE")
    group.add_option("--fixblens", action="store_true", dest="fixblens",
                     help="don't let codeml optimize branch lengths [default: %default]")
    group.add_option("--fixalpha", action="store_true", dest="fixalpha",
                     help="don't let codeml optimize gamma shape parameter [default: %default]")
    group.add_option("--initialalpha", action="store", type="float", dest="initialalpha",
                     help="set codeml's initial gamma shape parameter to NUM [default: %default]",
                     metavar="NUM")

    optparser.add_option_group(group)


    optparser.set_defaults(verbose=True,
                           clean=True,
                           codeml=CODEML,
                           baseml=BASEML,
                           pamlmodel=SUBMOD_NAME,
                           fixblens=False,
                           fixalpha=False,
                           initialalpha=INITIAL_ALPHA)

    (options, args) = optparser.parse_args()

    if len(args) < 2:
        optparser.error("incorrect number of arguments: ALIGNMENT.fasta PHYLOGENY.tre required")

    for i in range(len(args)):
        if args[i].find("./") == 0:
            args[i] = args[i][2:]

    alnfname = args[0]
    phyfname = args[1]

    ### sort out the working directory                           ###

    # get full path to phylogeny                                   #
    # this is needed in case the phylogeny is in a different       #
    # directory from the alignment. we will work in the alignment  #
    # directory                                                    #
    phyfname = os.path.abspath(phyfname)

    # figure out where the alignment file is and move there        #
    # this is necessary, as paml writes particular temporary files #
    workingdir = os.path.abspath(os.path.dirname(alnfname))
    alnfname   = os.path.basename(alnfname)

    returndir = os.getcwd()
    os.chdir(workingdir)

    ### done sorting out the working directory                   ###


    ### sort out ASR parameters                                  ###

    CODEML      = options.codeml
    BASEML      = options.baseml

    SUBMOD_NAME   = options.pamlmodel
    INITIAL_ALPHA = options.initialalpha

    if options.fixblens:
        FIX_BLENS = 2
    else:
        FIX_BLENS = 1

    if options.fixalpha:
        FIX_ALPHA = 1
    else:
        FIX_ALPHA = 0

    ### done sorting out ASR parameters                          ###


    # get list of all-gap columns, and collapse alignment #
    # Note that PAML is okay with all-gap columns, but    #
    # having too many all-gap columns might               #
    # cause indel reconstruction to be off, so we'll      #
    # remove them in any case. We'll re-insert them later #

    # read original alignment #
    orig_alignment = {}
    orig_aln_len = 0
    handle = open(alnfname, "r")
    for line in handle:
        if line[0] == ">":
            id = line[1:].strip()
            orig_alignment[id] = ""
        else:
            orig_alignment[id] += line.strip().replace(" ", "")
            orig_aln_len = len(orig_alignment[id])
    handle.close()

    # calculate colums having only ambiguous characters #
    all_ambig_cols = []
    for i in range(orig_aln_len):
        all_gaps = True
        for orig_seq in orig_alignment.values():
            if orig_seq[i] != "-" and orig_seq[i] != "X" and orig_seq[i] != "x":
                all_gaps = False
                break
        if all_gaps:
            all_ambig_cols.append(i)

    # write temporary compact alignment #
    compact_aln_fname = "ASRTEMP_compactalignment_ASRTEMP.fasta"
    handle = open(compact_aln_fname, "w")
    for (id,seq) in orig_alignment.items():
        handle.write(">%s\n" % id)
        for i in range(orig_aln_len):
            if i not in all_ambig_cols:
                handle.write(seq[i])
        handle.write("\n")
    handle.close()


    ### execute ASR ###
    (tree, mlseqs, asrprobs) = doASR(compact_aln_fname, phyfname, options.verbose, options.clean)


    ### print output files ###

    basename = alnfname.split(".fasta")[0]
    treefname   = basename + ".asr.tre"
    mlseqsfname = basename + ".asr.mlseqs.fasta"
    probsfname  = basename + ".asr.probdists.csv"

    handle = open(treefname, "w")
    handle.write("%s\n" % tree)
    handle.close()

    if options.verbose:
        sys.stdout.write("wrote labelled tree to %s/%s\n" % (workingdir,treefname))

    handle = open(mlseqsfname, "w")
    for (id,se) in mlseqs.items():
        handle.write(">%s\n" % id)
        col_id = 0
        for orig_col_id in range(orig_aln_len):
            if orig_col_id in all_ambig_cols:
                handle.write("-")
            else:
                handle.write(se[col_id])
                col_id += 1
        handle.write("\n")
    handle.close()

    if options.verbose:
        sys.stdout.write("wrote aligned ancestral sequences to %s/%s\n" % (workingdir,mlseqsfname))

    handle = open(probsfname, "w")
    handle.write("Node,Position,%s" % GAPS[0])
    for res in AMINO_ACIDS:
        handle.write(",%s" % res)
    handle.write(",*-ambiguousCol")
    handle.write("\n")
    for (id,pdists) in asrprobs.items():
        col_id = 0
        for orig_col_id in range(orig_aln_len):
            handle.write("%s,%d" % (id, orig_col_id))
            if orig_col_id in all_ambig_cols:
                handle.write(",%f" % 1.0)
                handle.write((",%f" % 0.0) * len(AMINO_ACIDS))
                handle.write(",*")
                handle.write("\n")
            else:
                pdist = pdists[col_id]
                for v in pdist:
                    handle.write(",%f" % v)
                handle.write("\n")
                col_id += 1
    handle.close()

    if options.verbose:
        sys.stdout.write("wrote posterior probability distributions to %s/%s\n" % (workingdir,probsfname))


    ### clean up ###

    if options.clean:
        if os.path.exists(compact_aln_fname):
            os.remove(compact_aln_fname)


    # return back home before exiting #
    os.chdir(returndir)

    if options.verbose:
        sys.stdout.write("finished.\n")

    sys.exit(0)

#
### END MAIN ###################################################################
