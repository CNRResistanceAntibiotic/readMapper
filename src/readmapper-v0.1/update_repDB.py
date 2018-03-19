#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/06/2017

import argparse
from Bio import SeqIO
import os


def load_rep_db(inpfile):
    repList = []
    idList  = []
    seqList = []
    f = open(inpfile, 'r')
    for n,rec in enumerate(SeqIO.parse(f, 'fasta')):
        if rec.id in idList:
            rec.id = rec.id + '_%i' % n
        idList.append(rec.id)
        if str(rec.seq) not in seqList:
            repList.append(rec)
            seqList.append(str(rec.seq))
    f.close()
    print '\nLoading of %s done!' % inpfile
    print 'Number of records: %i' % len(repList)
    return repList


def write_ariba_db(repList, db_name, db_ver, outdir, db):

    txt = ''
    for rec in repList:
        ID   = rec.id
        mol_type = 0
        comment = rec.id
        search_type = 0
        txt = txt + '%s\t%s\t%s\t%s\t%s\t%s\n' % (ID, mol_type, search_type, '.', '.', comment)

    outPrefix = '%s_%s_%s' % (os.path.join(outdir, db_name), db_ver, db)
    SeqIO.write(repList, open(outPrefix + '.fa','w'), 'fasta')

    f = open(outPrefix + '.tsv','w')
    f.write(txt)
    f.close()
    return outPrefix


def prepareref(outPrefix):

    wkDir = os.path.dirname(outPrefix)
    if os.path.exists(outPrefix) == True:
        cmd = 'rm -Rf %s' % outPrefix
        os.system(cmd)
    txt = '#!/bin/bash\n'
    txt = txt + 'ariba prepareref -f %s.fa -m %s.tsv %s\n' % (outPrefix, outPrefix, outPrefix)

    exeFile = os.path.join(wkDir, 'updaterepDB.sh')
    f = open(exeFile,'w')
    f.write(txt)
    f.close()
    cmd = 'chmod a+x %s' % exeFile
    #print cmd
    os.system(cmd)
    cmd = '%s' % exeFile
    #print cmd
    #os.system(cmd)
    #cmd = 'rm %s' % exeFile
    #os.system(cmd)

def main(args):

    repFile  = os.path.abspath(args.repFile)
    if '_' in os.path.basename(repFile):
        virname  = os.path.basename(repFile).split('_')[0]
    else:
        virname  = os.path.splitext(os.path.basename(repFile))[0]
    version  = args.versionDB
    db       = os.path.splitext(os.path.basename(repFile))[0].split('_')[-1]
    outDir  = os.path.dirname(repFile)

    repList = load_rep_db(repFile)
    outPrefix = write_ariba_db(repList, virname, version, outDir, db)
    prepareref(outPrefix)


def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='Write virDB for ariba - Version ' + version())
    parser.add_argument('-rep', '--repFile', dest="repFile", default='/usr/local/readmapper-v0.1/dbPLM/I1-5_pilL.fas', help="Fasta repDB")
    parser.add_argument('-v', '--versionDB', dest="versionDB", default="1", help="Version of repDB [1]")
    parser.add_argument('-V', '--version', action='version', version='write_virDB-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
