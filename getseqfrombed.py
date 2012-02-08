#!/usr/bin/env python
"""
This script is used to extract sequences from bed file.

USAGE 
  getseqfrombed.py {OPTIONS} <input .bed file> <reference .fasta file> 

OPTIONS

  -b [error file]\tSpecify the positional error profile to be used. The file should include at least 100 lines, each containing a positive number. The number at line x is the weight that an error is occured at x% position of the read. If no positional error file specified, uniform weight is assumed.

  -r [error rate]\tSpecify the overall error rate, a number between 0 and 1. Default is 0 (no errors).

  -l [read length]\tSpecify the read length. Default is 75. 

  -f [seq]\tFill at the end of each read by the sequence seq, if the read is shorter than the read length. This option is useful when simulating poly-A tails in RNA-Seq reads.

NOTE

  1. The input .bed file is best to sort according to chromosome names. Use - to input from STDIN.
  2. Biopython and numpy package are required.
  
  3. When applying models, we assume that all sequences are in the same length. The length information is given by the -l parameter. If the sequence length is greater than read length, nucleotides outside the read length will not be simulated for error.

HISTORY

  Aug 25, 2011:
    1.Rename makebedseq.py to getseqfrombed.py
    2.Print results to stdout
"""

import sys;
import pydoc;
import os;
import random;
import bisect;
import math;
import numpy;
from Bio import SeqIO;
from Bio.SeqRecord import SeqRecord;

if len(sys.argv)<2:
  print>>sys.stderr, (pydoc.render_doc(sys.modules[__name__]));
  sys.exit();

# analyzing parameters
posweight=[];
errrate=0.00;
readlength=75;
forcelength=False;
filledseq='A';

for i in range(len(sys.argv)):
  if i<len(sys.argv)-1:
    if sys.argv[i]=='-b':
      bline=0;
      tbweight=0;
      print>>sys.stderr, ('Using pos bias file'+sys.argv[i+1]);
      for lines in open(sys.argv[i+1]):
        bline=bline+1;
        if bline>100:
          break;
        tbweight=float(lines.strip());
        posweight.append(tbweight);
      if len(posweight)!=100:
        print>>sys.stderr, ('Error: the bias file should include at least 100 lines.');
        sys.exit();
    if sys.argv[i]=='-r':
      errrate=float(sys.argv[i+1]);
      if errrate<0 or errrate>1:
        print>>sys.stderr, ('Error: the error rate should be between 0-1.');
        sys.exit();
      print>>sys.stderr,('Error rate: '+str(errrate));
    if sys.argv[i]=='-l':
      readlength=int(sys.argv[i+1]);
      print>>sys.stderr,('Read length:'+str(readlength));
    if sys.argv[i]=='-f':
      forcelength=True;
      filledseq=sys.argv[i+1];
      print>>sys.stderr,('Force same read length with filled :'+(filledseq));
    


# construct weight probability for read length, if possible
rlenweight=[];
if len(posweight)!=0:
  kweight=0;
  for i in xrange(readlength):
    nfrac=i*100.0/readlength;
    lower=int(math.floor(nfrac));
    higher=int(math.ceil(nfrac));
    if higher==lower: higher=lower+1;
    #print('higher:'+str(higher)+',lower:'+str(lower));
    if higher<100:
      val=posweight[lower]*(nfrac-lower)+posweight[higher]*(higher-nfrac);
    else:
      val=posweight[99];
    kweight+=val;
    rlenweight.append(kweight);

bedfile=sys.argv[-2];
reffile=sys.argv[-1];
#ofastafile=sys.argv[-1];

# build reference
seqref=SeqIO.index(reffile,'fasta');
refkeys=seqref.keys();

# read bed file, and ready for writing
if bedfile!="-":
  fid=open(bedfile);
else:
  fid=sys.stdin;
#ofid=open(ofastafile,'w');
ofid=sys.stdout

nlines=0;

prevchr='';
previndex='';

for lines in fid:
  # update line counter
  nlines=nlines+1;
  if nlines %10000==1:
    print>>sys.stderr,('Processing '+str(nlines)+' lines...');
  # parse lines
  bedfield=lines.strip().split('\t');
  if len(bedfield)!=12:
    print>>sys.stderr,('Error: incorrect number of fields at line %d (should be 12, observed %d)' % (nlines, len(bedfield)) );
    continue;
  # clustering
  fieldrange=[int(bedfield[1]),int(bedfield[2])];
  # parse all exons
  exonlen=[int(x) for x in bedfield[10][:-1].split(',')];
  exonstart=[int(x)+fieldrange[0] for x in bedfield[11][:-1].split(',')];
  if not bedfield[0] in refkeys:
    print>>sys.stderr,('Warning: '+bedfield[0]+ ' not in the reference. Ignore...' );
    continue;
  if bedfield[0]!=prevchr:
    print>>sys.stderr, ('Switching to %s ...' % bedfield[0]);
    prevchr=bedfield[0];
    previndex=seqref[bedfield[0]];
  # extract sequences
  thisseq=SeqRecord('');
  for i in range(len(exonlen)):
    thisseq+=previndex[exonstart[i]:(exonstart[i]+exonlen[i])];
  if forcelength:
    if sum(exonlen)<readlength:
      thisseq+=filledseq*(readlength-sum(exonlen));
  thisseq.id=bedfield[3];
  thisseq.description='';
  # reverse-complement the sequence if it is on the negative strand
  if bedfield[5]=='-':
    thisseq.seq=thisseq.seq.reverse_complement();
  # mutation
  nmut=numpy.random.poisson(errrate);
  if nmut>0:
    newseq=thisseq.seq;
    for n in xrange(nmut):
      if len(posweight)==0:
        # uniform distrib
        modifyposition=random.choice(xrange(len(newseq)));
      else:
        rchosen=random.random()*kweight;
        modifyposition=bisect.bisect_right(posweight,rchosen);
      # mutate the position
      if len(newseq)>modifyposition:
        topos=random.choice('ATGC');
        while topos==newseq[modifyposition]:
          topos=random.choice('ATGC');
        print >>sys.stderr,('MUTATION at position '+str(modifyposition)+','+newseq[modifyposition]+'->'+topos);
        newseq=newseq[:modifyposition]+topos+newseq[(modifyposition+1):];
    #print>>sys.stderr,('NMUTATION:'+str(nmut));
    #print>>sys.stderr, (str(thisseq.seq));
    #print>>sys.stderr,(newseq);
    #thisseq.seq=newseq;
  # write to record
  SeqIO.write(thisseq,ofid,'fasta');
  
        

# ofid.close();
if bedfile!="-":
  fid.close();
    



