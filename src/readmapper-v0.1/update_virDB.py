#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/06/2017

import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os


def load_vir_db(inpfile):
    virDic = {}
    f = open(inpfile, 'r')
    for n, line in enumerate(f):
        line = line.strip()
        if n == 0:
            header = line.split('\t')
        elif line != '' and n > 0:
            line = line.split('\t')
            data = dict(zip(header,line))
            for key in header:
                if key not in data.keys():
                    data[key]=''
            key = data['gene_id']
            virDic[key] = data
    f.close()
    print '\nLoading of %s done!' % inpfile
    print 'Number of records: %i' % len(virDic.keys())
    return virDic


def write_ariba_db(virDic, db_name, db_ver, outdir, db=['all']):
    keys = virDic.keys()
    records = []
    txt = ''
    keys.sort()
    for key in keys:
        data = virDic[key]
        if data['genus'].lower() in db or db == ['all']:
            gene_name = data['gene_name']
            dna_seq, prot_seq = dna_translate(Seq(data['VF_sequence']))

            ID   = '%s::%s' % (key, gene_name.replace(' ','_'))
            record = SeqRecord(id=ID, name=ID, description='', seq=Seq(dna_seq))
            records.append(record)

            if prot_seq != '':
                mol_type = 1
            else:
                mol_type = 0
            genus = data['genus']
            species = data['species']
            strain = data['strain']
            vf_name = data['VF_Name']
            vf_accession = data['VF_Accession']
            comment = 'genus:%s,species:%s,strain:%s,gene_name:%s,vf_name:%s,vf_accession:%s' % \
                      (genus, species, strain, gene_name, vf_name, vf_accession)

            search_type = 0
            txt = txt + '%s\t%s\t%s\t%s\t%s\t%s\n' % (ID, mol_type, search_type, '.', '.', comment)

    db = '-'.join(db)
    if db == 'escherichia-shigella-salmonella-klebsiella-enterobacter-serratia-yersinia-citrobacter-proteus-pantoea':
        db = 'enterob'

    outPrefix = '%s_%s_%s' % (os.path.join(outdir, db_name), db_ver, db)
    SeqIO.write(records, open(outPrefix + '.fa','w'), 'fasta')
    print 'Number of selected records: %i' % len(records)

    f = open(outPrefix + '.tsv','w')
    f.write(txt)
    f.close()
    return outPrefix


def prepareref(outPrefix):

    txt = '#!/bin/bash\n'
    txt = txt + 'ariba prepareref -f %s.fa -m %s.tsv %s\n' % (outPrefix, outPrefix, outPrefix)

    f = open('update_virDB.sh','w')
    f.write(txt)
    f.close()

    os.system('chmod a+x update_virDB.sh')
    os.system('./update_virDB.sh')
    os.system('rm ./update_virDB.sh')


def dna_translate(dna_sequence):

    try:
        prot_sequence = dna_sequence.translate(table='Bacterial', cds=True)
    except:
        try:
            prot_sequence = dna_sequence[:-1].translate(table='Bacterial', cds=True)
            dna_sequence = dna_sequence[:-1]
        except:
            try:
                prot_sequence = dna_sequence.reverse_complement().translate(table='Bacterial', cds=True)
                dna_sequence = dna_sequence.reverse_complement()
            except:
                try:
                    prot_sequence = dna_sequence[:-1].reverse_complement().translate(table='Bacterial', cds=True)
                    dna_sequence = dna_sequence[:-1].reverse_complement()
                except:
                    prot_sequence = ''
    return str(dna_sequence), str(prot_sequence)


def main(args):

    virFile  = os.path.abspath(args.virFile)
    virname  = os.path.splitext(os.path.basename(virFile))[0]
    version  = args.versionDB
    subsetDB = args.subDB
    if subsetDB == '':
        subsetDB = ['all']
    else:
        subsetDB = subsetDB.split(',')

    outDir  = os.path.dirname(virFile)
    virDic = load_vir_db(virFile)
    outPrefix = write_ariba_db(virDic, virname, version, outDir, subsetDB)
    prepareref(outPrefix)


def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='Write virDB for ariba - Version ' + version())
    parser.add_argument('-vir', '--armFile', dest="virFile", default='/usr/local/readmapper-v0.1/dbVIR/virDB_A.csv', help="CSV virDB")
    parser.add_argument('-sub', '--subDB', dest="subDB", default='escherichia,shigella,salmonella,klebsiella,enterobacter,serratia,yersinia,citrobacter,proteus,pantoea', help="Subset of armDB as comas-separated genus list")
    parser.add_argument('-v', '--versionDB', dest="versionDB", default="1", help="Version of virDB [1]")
    parser.add_argument('-V', '--version', action='version', version='write_virDB-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
