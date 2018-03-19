#!/usr/bin/python
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from prepare_mapping import read_sample_file
from prepare_mapping import read_setting
import pandas as pd


def load_repDB(inpfile):
    print 'database used: %s' % inpfile
    gene_repDic = {}
    f = open(inpfile, 'r')
    for n, rec in enumerate(SeqIO.parse(f, 'fasta')):
        gene_repDic[rec.id] = rec
    f.close()

    print '\nLoading of %s done!' % inpfile
    print 'Number of rep sequence: %i' % len(gene_repDic.keys())
    return gene_repDic


def load_repRes(tsvFile, genFile, seqFile):

    if os.path.exists(genFile):
        genFile  = convert_Fagz_file(genFile)
    else:
        genFile = os.path.splitext(genFile)[0]
    genDic = read_fa_file(genFile)

    if os.path.exists(seqFile):
        seqFile = convert_Fagz_file(seqFile)
    else:
        seqFile = os.path.splitext(seqFile)[0]
    seqDic =  read_fa_file(seqFile)

    resDic = {}
    f = open(tsvFile)
    for n,line in enumerate(f):
        if n == 0:
            header = line.strip().split('\t')
        else:
            dtDic = dict(zip(header, line[:-1].split('\t')))
            ref_name = dtDic['ref_name']
            #gene     = dtDic['gene'] #1 gene 0 non-coding
            
            found = '0'
            for key in genDic.keys():
                if dtDic['ctg'] in key:
                    found = '1'
                    dna_rec = genDic[key]
                    break
            if found == '0':
                for key in seqDic.keys():
                    if dtDic['ref_name'] in key:
                        dna_rec = seqDic[key]
                        break

            ref_len = float(dtDic['ref_len'])
            ref_base_assembled = int(dtDic['ref_base_assembled'])
            cov = (ref_base_assembled / ref_len) * 100
            dtDic['pc_coverage'] = '%.2f' % cov

            update = '1'
            if ref_name in resDic.keys():
                if cov > float(resDic[ref_name.split('__')[0]]['pc_coverage']):
                    update = '1'
                else:
                    update = '0'

            if update == '1':
                key = ref_name.split('__')[0]
                resDic[key] = {'pc_identity':dtDic['pc_ident'],'pc_coverage':dtDic['pc_coverage'],
                                'gene':dtDic['gene'], 'flag':dtDic['flag'], 'ctg_name':dtDic['ctg'],
                                #'ctg_len':dtDic['ctg_len'],
                                'ref_name': dtDic['ref_name'],
                                'mean_depth':dtDic['ctg_cov'], 'dna_sequence':dna_rec}
    f.close()
    return resDic


def convert_Fagz_file(fagzFile):
    cmd = 'gunzip %s' % fagzFile
    os.system(cmd)
    faFile = os.path.splitext(fagzFile)[0]
    return faFile


def read_fa_file(faFile):
    faDic = {}
    for rec in SeqIO.parse(open(faFile), 'fasta'):
       faDic[rec.id] = rec
    return faDic


def filter_results(resDic, resuFile, passcov=70, passid=65):

    f = open(resuFile,'a')
    del_keys = []
    for key in resDic.keys():
        # filter target coverage < passcov
        if float(resDic[key]['pc_coverage']) < passcov:
            del_keys.append(key)
            f.write('Filter pc_coverage < %.2f\n:' % passcov)
            f.write('%s\n%s\n\n' % (key, resDic[key]))
            # Alarm with id >= 90 and passid > 40
            if float(resDic[key]['pc_identity']) >= passid \
                    and float(resDic[key]['pc_coverage']) > 50:
                print '\n'
                print '############################################################'
                print '###  WARNING PUTATIVE MIS-DETECTION: Low coverage alert  ###'
                print '############################################################'
                print ''
                print 'Record:', key, '\tCoverage:', resDic[key]['pc_coverage'], \
                    '\tIdentity:', resDic[key]['pc_identity'], '\tMean depth:', resDic[key]['mean_depth']
                print ''
                #print resDic[key]
                #print ''
                print 'The record will be deleted in the final results'
                print ''

        # filter id percentage < passid
        if float(resDic[key]['pc_identity']) < passid:
            del_keys.append(key)
            f.write('Filter pc_identity < %.2f\n:' % passid)
            f.write('%s\n%s\n\n' % (key, resDic[key]))
            # Alarm if id >= 40 and cov >=75
            if float(resDic[key]['pc_identity']) >= 40 \
                    and float(resDic[key]['pc_coverage']) >= 80:
                print '\n'
                print '############################################################'
                print '###  WARNING PUTATIVE MIS-DETECTION: Low identity alert  ###'
                print '############################################################'
                print ''
                print 'Record:', key, '\tCoverage:', resDic[key]['pc_coverage'], \
                    '\tIdentity:', resDic[key]['pc_identity'], '\tMean depth:', resDic[key]['mean_depth']
                print ''
                #print resDic[key]
                #print ''
                print 'The record will be deleted in the final results'
                print ''

        if float(resDic[key]['mean_depth']) <= 15 and float(resDic[key]['pc_identity']) >= passid \
                and float(resDic[key]['pc_coverage']) > passcov:
            print '\n'
            print '####################################################################'
            print '###  WARNING PUTATIVE MIS-DETECTION: < 15 sequencing depth alert  ###'
            print '####################################################################'
            print ''
            print 'Record:', key, '\tCoverage:', resDic[key]['pc_coverage'], \
                '\tIdentity:', resDic[key]['pc_identity'], '\tMean depth:', resDic[key]['mean_depth']
            print ''
            #print resDic[key]
            #print ''
            print 'The record will be kept in the final results'
            print ''
    f.close()
    # Deletion of records
    #print 'Filtered results (%s): %i' % (species, len(del_keys))
    for key in list(set(del_keys)):
        del resDic[key]
    return resDic


def write_csv_result(resDic, sample_id, outdir, dtbasename):

    f = open(os.path.join(outdir, 'results_%s.csv' % dtbasename), 'w')
    f.write('sample_id\treplicon\tmean_depth\tpc_coveraget\tpc_identity\tSequence\n')
    rep_list = resDic.keys()
    rep_list.sort()
    for rep_id in rep_list:
        mean_depth = resDic[rep_id]['mean_depth']
        pc_coverage = resDic[rep_id]['pc_coverage']
        pc_identity = resDic[rep_id]['pc_identity']
        dna_sequence = str(resDic[rep_id]['dna_sequence'].seq)
        txt = '%s\t%s\t%s\t%s\t%s\t%s\n' % (sample_id, rep_id, mean_depth, pc_coverage, pc_identity, dna_sequence)
        f.write(txt)
    f.close()


def write_summary_result(resDic, outdir, dtbasename, sampleID):

    writer = pd.ExcelWriter(os.path.join(outdir, 'summary_results_%s.xlsx' % dtbasename))
    rep_list = resDic.keys()
    rep_list.sort()
    data = {}
    for rep in rep_list:
        data[rep] = '%s [id:%.2f;cv:%.2f;dp:%.2f]' % (rep, float(resDic[rep]['pc_identity']),float(resDic[rep]['pc_coverage']),float(resDic[rep]['mean_depth']))
    df = pd.DataFrame(data, index=[sampleID])
    df.to_excel(writer, 'replicons', index=True)
    writer.save()

def main(args):
    sampleID    = args.sampleID
    samplefile  = args.sampleFile
    settingfile = args.settingFile
    dtbasetype  = args.dtbase
    wkdir       = args.wkDir

    if samplefile == '':
        samplefile = os.path.join(wkdir, 'sample.csv')
    if wkdir == '':
        wkdir = os.path.dirname(samplefile)
    outdir = os.path.join(wkdir, sampleID)

    setDic = read_setting(settingfile)
    sampleDic, sampleList = read_sample_file(samplefile)
    species = sampleDic[sampleID]
    set_species = setDic[species.lower()]
    dtbasename = '%s_%s' % (set_species[dtbasetype][0], set_species[dtbasetype][1])
    dtbasedir = setDic['dbDir'][dtbasetype][0]
    dtbasefile = os.path.join(dtbasedir, dtbasename + '.fa')
    #clusterfile = os.path.join(dtbasedir, '_'.join(dtbasename.split('_')[:2]) + '.clu')

    tsvFile = os.path.join(outdir, dtbasename, 'report.tsv')
    genFile = os.path.join(outdir, dtbasename, 'assembled_genes.fa.gz')
    seqFile = os.path.join(outdir, dtbasename, 'assembled_seqs.fa.gz')

    resDic = load_repRes(tsvFile, genFile, seqFile)
    print 'Number of results: %i' % len(resDic)
    gene_virDic = load_repDB(dtbasefile)

    #Aurelien: j ai rajoute un path defini pour ce fichier pour eviter qu il soit dans le rep courant
    resuFile = os.path.join(outdir, dtbasename, 'filtered_out_rep_results.txt')    

    resDic = filter_results(resDic, resuFile, 80, 95)
    print 'Number of parsed results: %i' % len(resDic)

    write_csv_result(resDic, sampleID, outdir, dtbasename)
    write_summary_result(resDic, outdir, dtbasename, sampleID)


def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='parse_rep_detection - Version ' + version())
    parser.add_argument('-s', '--sampleID', dest="sampleID", help='Sample ID')
    parser.add_argument('-sf','--sampleFile', dest="sampleFile", default='', help='Sample file')
    parser.add_argument('-d', '--dtbase', dest="dtbase", default='rep', help="Database type")
    parser.add_argument('-wd','--wkDir',      dest="wkDir", default='/home/bacteriologie/ariba/anses/mcr-3', help="Working directory")
    parser.add_argument('-st','--settingFile',dest="settingFile",default='/usr/local/readmapper-v0.1/setting.txt', help="Setting file")
    parser.add_argument('-v', '--verbose',    dest="verbose",    default="0",help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version',    action='version', version='parse_rep_detection-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
