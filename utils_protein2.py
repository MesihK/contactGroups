import sys
from protein import protein
from AAmap import AAmap
import os.path
from collections import defaultdict
import numpy as np
import commp as cp


def seqaafreq(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_protein2.py seqaafreq 1t3r.pdb all (output: 1t3r.pdb.all.saf)')
	pdbfile = arglist[0]
	chainid = arglist[1]

	p = protein(pdbfile, chain=chainid)
	seqaa = p.seq.upper()
	aafreq = defaultdict(lambda:0)
	for aa in seqaa: 
		aafreq[aa]+=1

	outstr = ' '.join([('%d' % aafreq[aa]) for aa in cp.aat01])
	outfile = '%s.%s.saf' % (pdbfile, chainid)
	with open(outfile ,'w') as fout:
		fout.write('%s\n' % outstr)
	cp._info('save to %s' % outfile)

# write cutoff contact by specific method
# output: 1. method residue pdb file
# 		  2. contact file	
#def contactbycutoff(arglist):
def writecontact(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_protein2.py contactbycutoff 1t3r.pdb chain sgc cutoff')
	pdbfile = arglist[0] #sys.argv[2]
	chainid = arglist[1] #sys.argv[3]
	method = arglist[2] #sys.argv[4]
	cutoff = float(arglist[3]) #float(sys.argv[5])

	if method not in ['sgc', 'tip', 'ca']:
		cp._err('invalid method: %s' % method)

	p = protein(pdbfile, chain=chainid)
	if method =='sgc':
		ralist = p.atomsbyscgmcenter()
	elif method == 'tip':
		ralist = p.atomsbytip('AAtips.py')
	elif method == 'ca':
		ralist = p.ca
	# continue ...
	cglist = p.contactbycutoff(ralist, cutoff)

	outrafile = '%s.%s.%s' % (pdbfile, chainid, method)
	with open(outrafile, 'w') as fp:
		for a in ralist:
			fp.write(a.writeAtom())
	cp._info('save residue pdb: %s' % outrafile)

	outcgfile = '%s.%s.%s.cg' % (pdbfile, chainid, method)
	with open(outcgfile, 'w') as fp:
		for a,b in cglist:
			fp.write('%d %s %d %s\n' % (a.resSeq, cp.aa2a[a.resName], b.resSeq, cp.aa2a[b.resName]))
	cp._info('save cg file: %s' % outcgfile)

def foo(arglist):
	cp._info(repr(arglist))

if __name__ == '__main__':
	cp.dispatch(__name__)
