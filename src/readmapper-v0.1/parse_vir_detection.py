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


def load_virDB(inpfile):
    print 'database used: %s' % inpfile
    gene_virDic = {}
    vf_virDic = {}
    f = open(inpfile, 'r')
    for n, line in enumerate(f):
        line = line.strip()
        if n == 0:
            header = line.strip().split('\t')
        elif line != '' and n > 0:
            line = line.split('\t')[:-1]
            data = dict(zip(header, line))

            vf  = data['VF_Accession']
            strain  = '%s__%s__%s' % (data['genus'], data['species'], data['strain'].replace(' ','_'))

            gene_virDic[data['gene_id']] = data
            try:
                vf_virDic[vf][strain].append(data['gene_id'])
            except KeyError:
                if vf not in vf_virDic.keys():
                    vf_virDic[vf] = {strain : [data['gene_id']]}
                else:
                    vf_virDic[vf][strain] = [data['gene_id']]
    f.close()

    print '\nLoading of %s done!' % inpfile
    print 'Number of virulence genes: %i' % len(gene_virDic.keys())
    print 'Number of virulence factors: %i' % len(vf_virDic)
    return gene_virDic, vf_virDic


def load_cluDB(clusterfile):
    cluDic = {}
    for line in open(clusterfile):
        line = line.strip().split('\t')

        key = line[0].split('::')[0]
        gene_id = line[1].split('::')[0]

        if len(line) > 2:
            if key not in cluDic.keys():
                cluDic[key] = {gene_id:[]}
            if gene_id not in cluDic[key].keys():
                cluDic[key][gene_id] = []

            homologs = line[2].split(',')
            for homolog in homologs:
                gene_id2 = homolog.split('::')[0]
                cluDic[key][gene_id].append(gene_id2)

    return cluDic


def load_virRes(tsvFile, genFile, seqFile):

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
            gene     = dtDic['gene'] #1 gene 0 non-coding
            
            found = '0'
            for key in genDic.keys():
                if dtDic['ctg'] in key:
                    found = '1'
                    dna_rec = genDic[key]
                    prot_rec = str(translate_dna(dna_rec.seq))
                    break
            if found == '0':
                for key in seqDic.keys():
                    if dtDic['ref_name'] in key:
                        found = '1'
                        dna_rec = seqDic[key]
                        prot_rec = '.'
                        break

            if found == '0':
                print '\nSequence not found:', dtDic['ctg']
                print dtDic
                exit(1)
            
            ref_len = float(dtDic['ref_len'])
            ref_base_assembled = int(dtDic['ref_base_assembled'])
            cov = (ref_base_assembled / ref_len) * 100
            dtDic['pc_coverage'] = '%.2f' % cov

            update = '1'
            if ref_name in resDic.keys():
                if cov > float(resDic[ref_name]['pc_coverage']):
                    update = '1'
                else:
                    update = '0'

            if update == '1':
                key = ref_name.split('__')[0]
                resDic[key] = {'pc_identity':dtDic['pc_ident'],'pc_coverage':dtDic['pc_coverage'],
                                'gene':dtDic['gene'], 'flag':dtDic['flag'], 'ctg_name':dtDic['ctg'],
                                #'ctg_len':dtDic['ctg_len'],
                                'ref_name': dtDic['ref_name'],
                                'mean_depth':dtDic['ctg_cov'], 'dna_sequence':dna_rec, 'prot_sequence':prot_rec}
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


def filter_results(resDic, passcov=70, passid=70):
    f = open('filtered_out_vir_results.txt','a')
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


def score_results(resDic, gene_virDic, vf_virDic, cluDic):

    vf_list = []
    for key in resDic.keys():
        vf_list.append(gene_virDic[key]['VF_Accession'])
        #print key, gene_virDic[key]['VF_Accession']
    vf_list = list(set(vf_list))

    scoreDic = {}
    for vf in vf_list:
        strain_list = vf_virDic[vf].keys()
        strain_list.sort()
        for key in strain_list:
            total_gene = 0
            positive_gene = 0
            gene_list = vf_virDic[vf][key]
            gene_list.sort()
            geneDic = {}
            for gene_id in gene_list:
                total_gene += 1
                try:
                    pc_coverage = resDic[gene_id]['pc_coverage']
                    pc_identity = resDic[gene_id]['pc_identity']
                    positive_gene += 1
                except KeyError:
                    pc_coverage = '.'
                    pc_identity = '.'

                geneDic[gene_id] = {'pc_coverage': pc_coverage, 'pc_identity': pc_identity}

            if positive_gene > 0:
                if vf in scoreDic.keys():
                    scoreDic[vf][key]={'gene_dic':geneDic}
                    scoreDic[vf][key]['score']={'positive_gene':positive_gene, 'total_gene':total_gene}
                else:
                    scoreDic[vf]={key:{'gene_dic':geneDic}}
                    scoreDic[vf][key]['score']={'positive_gene':positive_gene, 'total_gene':total_gene}    

    print_scoreDic(scoreDic, vf_virDic, gene_virDic)

    print '\nUpdated results:'
    import copy
    updatedDic = copy.deepcopy(scoreDic)
    vf_list = scoreDic.keys()
    vf_list.sort()

    max_pos_vf = {}
    max_tgene_vf = {}
    for vf in vf_list:
        strain_list = scoreDic[vf].keys()
        strain_list.sort()

        for strain in strain_list:
            positive_gene = scoreDic[vf][strain]['score']['positive_gene']
            total_gene = scoreDic[vf][strain]['score']['total_gene']
            if vf not in max_pos_vf.keys():
                max_pos_vf[vf] = positive_gene
                max_tgene_vf[vf] = total_gene
            else:
                if positive_gene > max_pos_vf[vf]:
                    max_pos_vf[vf] = positive_gene
                if total_gene > max_tgene_vf[vf]:
                    max_tgene_vf[vf] = total_gene
            gene_list = vf_virDic[vf][strain]
            gene_list.sort()
            txt = ''

            for gene_id in gene_list:
                vf_name = gene_virDic[gene_id]['VF_Name']
                gene_name = gene_virDic[gene_id]['gene_name']


                try:
                    homolog_gene_ids = cluDic[vf][gene_id]
                except KeyError:
                    homolog_gene_ids = []

                homolog_gene_id = ''
                for homolog_gene_id in homolog_gene_ids:
                    if homolog_gene_id in resDic.keys():
                        #print homolog_gene_id, resDic[homolog_gene_id]
                        vf2 = gene_virDic[homolog_gene_id]['VF_Accession']
                        vf2_name = gene_virDic[homolog_gene_id]['VF_Name']
                        #print scoreDic[vf2]
                        gene2_name = gene_virDic[homolog_gene_id]['gene_name']
                        genus2 = gene_virDic[homolog_gene_id]['genus']
                        species2 = gene_virDic[homolog_gene_id]['species']
                        strain2 = gene_virDic[homolog_gene_id]['strain'].replace(' ','_')
                        key_strain2 = '%s__%s__%s' % (genus2, species2, strain2)
                        #positive_gene2 = scoreDic[vf2][key_strain2]['score']['positive_gene']
                        #total_gene2 = scoreDic[vf2][key_strain2]['score']['total_gene']
                        break
                    else:
                        homolog_gene_id = ''

                pc_coverage = str(scoreDic[vf][strain]['gene_dic'][gene_id]['pc_coverage'])
                pc_identity = str(scoreDic[vf][strain]['gene_dic'][gene_id]['pc_identity'])

                if pc_coverage == '.' and pc_identity == '.' and homolog_gene_id != '':
                    updatedDic[vf][strain]['gene_dic'][gene_id]['pc_coverage'] = scoreDic[vf2][key_strain2]['gene_dic'][homolog_gene_id]['pc_coverage']
                    updatedDic[vf][strain]['gene_dic'][gene_id]['pc_identity'] = scoreDic[vf2][key_strain2]['gene_dic'][homolog_gene_id]['pc_identity']
                    updatedDic[vf][strain]['gene_dic'][gene_id]['homolog'] = 'VF_id:%s,VF_name:%s,genus:%s,species:%s,strain:%s,gene_id:%s,gene_name:%s' \
                                                                       % (vf2, vf2_name, genus2, species2, strain2, homolog_gene_id, gene2_name)
                    pc_coverage = str(updatedDic[vf][strain]['gene_dic'][gene_id]['pc_coverage'])
                    pc_identity = str(updatedDic[vf][strain]['gene_dic'][gene_id]['pc_identity'])
                    updatedDic[vf][strain]['score']['positive_gene'] = updatedDic[vf][strain]['score']['positive_gene'] + 1
                    gene_name = '%s=%s' % (gene_name, gene2_name)

                txt = txt + ' gene: %s [%s/%s]' % (gene_name, pc_identity, pc_coverage)

            txt = txt[1:]
            print 'VF_id: %s\t%s\tVF_name: %s\tScore: %i/%i\t%s' % (vf, strain, vf_name, positive_gene, total_gene, txt)

    rmDic = {}
    kpDic = {}
    for vf in updatedDic.keys():
        #max_tgene = max_tgene_vf[vf]
        max_pos = max_pos_vf[vf]
        for strain in updatedDic[vf].keys():
            #total_gene = scoreDic[vf][strain]['score']['total_gene']
            positive_gene = scoreDic[vf][strain]['score']['positive_gene']
            #if total_gene < max_tgene:
            if positive_gene < max_pos:
                if vf not in rmDic.keys():
                    rmDic[vf] = [strain]
                else:
                    rmDic[vf].append(strain)
            else:
                if vf not in kpDic.keys():
                    kpDic[vf] = [strain]
                else:
                    kpDic[vf].append(strain)
    for vf in kpDic.keys():
        if len(kpDic[vf]) > 1:
            max_tgene = max_tgene_vf[vf]
            #max_pos = max_pos_vf[vf]
            for strain in kpDic[vf]:
                total_gene = scoreDic[vf][strain]['score']['total_gene']
                #positive_gene = scoreDic[vf][strain]['score']['positive_gene']
                if total_gene < max_tgene:
                #if positive_gene < max_pos:
                    if vf not in rmDic.keys():
                        rmDic[vf] = [strain]
                    else:
                        rmDic[vf].append(strain)

    for vf in rmDic.keys():
        for strain in rmDic[vf]:
            del updatedDic[vf][strain]
            if updatedDic[vf] == {}:
                print 'Deletion of vf %s' % vf
                del updatedDic[vf]

    print '\nPurged results:'
    print_scoreDic(updatedDic,vf_virDic,gene_virDic)

    return updatedDic


def translate_dna(dnaObject, table='Bacterial', cds=True):
    try:
        prot = dnaObject.translate(table=table, cds=cds)
    except:
        try:
            prot = dnaObject.reverse_complement().translate(table=table, cds=cds)
        except:
            prot = ''
    return prot


def print_scoreDic(scoreDic, vf_virDic, gene_virDic):
    vf_list = scoreDic.keys()
    vf_list.sort()
    for i,vf in enumerate(vf_list):
        strain_list = scoreDic[vf].keys()
        strain_list.sort()
        for strain in strain_list:
            positive_gene = scoreDic[vf][strain]['score']['positive_gene']
            total_gene = scoreDic[vf][strain]['score']['total_gene']
            score = '%i/%i' % (positive_gene,total_gene)
            gene_list = vf_virDic[vf][strain]
            gene_list.sort()
            txt = ''
            for gene_id in gene_list:
                vf_name = gene_virDic[gene_id]['VF_Name']
                gene_name = gene_virDic[gene_id]['gene_name']

                pc_coverage = str(scoreDic[vf][strain]['gene_dic'][gene_id]['pc_coverage'])
                pc_identity = str(scoreDic[vf][strain]['gene_dic'][gene_id]['pc_identity'])


                if 'homolog' in scoreDic[vf][strain]['gene_dic'][gene_id].keys():
                    homolog = scoreDic[vf][strain]['gene_dic'][gene_id]['homolog']
                    txt = txt + 'gene: %s [%s: %s/%s]\t' % (gene_name, homolog, pc_identity, pc_coverage)
                else:
                    txt = txt + 'gene: %s [%s/%s]\t' % (gene_name, pc_identity, pc_coverage)

            txt = txt[:-1]
            print '%s\t%s\t%s\t%s\t%s' % (vf, strain, vf_name, score, txt)


def write_csv_result(resDic, vf_virDic, gene_virDic, outdir, dtbasename):

    f = open(os.path.join(outdir, 'results_%s.csv' % dtbasename), 'w')
    f.write('VF_Accession\tStrain\tVF_Name\tScore\tGenes')
    vf_list = resDic.keys()
    vf_list.sort()
    for vf in vf_list:
        strain_list = resDic[vf].keys()
        strain_list.sort()
        for strain in strain_list:
            positive_gene = resDic[vf][strain]['score']['positive_gene']
            total_gene = resDic[vf][strain]['score']['total_gene']
            score = '%i/%i' % (positive_gene, total_gene)

            gene_list = vf_virDic[vf][strain]
            gene_list.sort()
            txt = ''
            for gene_id in gene_list:
                vf_name = gene_virDic[gene_id]['VF_Name']
                gene_name = gene_virDic[gene_id]['gene_name']

                pc_coverage = str(resDic[vf][strain]['gene_dic'][gene_id]['pc_coverage'])
                pc_identity = str(resDic[vf][strain]['gene_dic'][gene_id]['pc_identity'])

                if 'homolog' in resDic[vf][strain]['gene_dic'][gene_id].keys():
                    homolog = resDic[vf][strain]['gene_dic'][gene_id]['homolog']
                    txt = txt + 'gene:%s,homolog:%s,pc_id:%s,pc_cv:%s\t' % (gene_name, homolog, pc_identity, pc_coverage)
                else:
                    txt = txt + 'gene:%s,pc_id:%s,pc_cv:%s\t' % (gene_name, pc_identity, pc_coverage)
            txt = '%s\t%s\t%s\t%s\t%s\n' % (vf, strain, vf_name, score, txt)
            f.write(txt)
    f.close()


def write_summary_result(resDic, vf_virDic, gene_virDic, outdir, dtbasename, sampleID):
    from collections import OrderedDict
    writer = pd.ExcelWriter(os.path.join(outdir, 'summary_results_%s.xlsx' % dtbasename))
    vf_list = resDic.keys()
    vf_list.sort()
    for vf in vf_list:
        strain_list = resDic[vf].keys()
        strain_list.sort()
        for strain in strain_list:
            #positive_gene = resDic[vf][strain]['score']['positive_gene']
            #total_gene = resDic[vf][strain]['score']['total_gene']
            #score = '%i/%i' % (positive_gene, total_gene)

            gene_list = vf_virDic[vf][strain]
            gene_list.sort()
            data = OrderedDict()
            #data['sample_id'] = sampleID
            data['Strain'] = strain
            for gene_id in gene_list:
                vf_name = gene_virDic[gene_id]['VF_Name']
                gene_name = gene_virDic[gene_id]['gene_name']

                pc_identity = str(resDic[vf][strain]['gene_dic'][gene_id]['pc_identity'])
                if 'homolog' in resDic[vf][strain]['gene_dic'][gene_id].keys():
                    homologs = resDic[vf][strain]['gene_dic'][gene_id]['homolog'].split(',')
                    for item in homologs:
                        if 'species' in item:
                            hspecies = item.split(':')[1]
                        if 'genus' in item:
                            hgenus = item.split(':')[1]
                        if 'strain' in item:
                            hstrain = item.split(':')[1]
                    hstrain = '%s__%s__%s' % (hgenus, hspecies, hstrain)

                    pc_identity = pc_identity + ' [%s]' % hstrain

                data[gene_name] = pc_identity

            df = pd.DataFrame(data, index=[sampleID])
            vf_name = vf_name.replace('/',' ')
            if len(vf_name) >= 31:
                vf_name=vf_name[:30]
            df.to_excel(writer, vf_name, index=True, index_label='sample_id')
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
    dtbasefile = os.path.join(dtbasedir, '_'.join(dtbasename.split('_')[:2]) + '.csv')
    clusterfile = os.path.join(dtbasedir, '_'.join(dtbasename.split('_')[:2]) + '.clu')

    tsvFile = os.path.join(outdir, dtbasename, 'report.tsv')
    genFile = os.path.join(outdir, dtbasename, 'assembled_genes.fa.gz')
    seqFile = os.path.join(outdir, dtbasename, 'assembled_seqs.fa.gz')

    resDic = load_virRes(tsvFile, genFile, seqFile)
    print 'Number of results: %i' % len(resDic)
    gene_virDic, vf_virDic = load_virDB(dtbasefile)
    cluDic = load_cluDB(clusterfile)
    resDic = filter_results(resDic)
    print 'Number of parsed results: %i' % len(resDic)
    resDic = score_results(resDic, gene_virDic, vf_virDic, cluDic)

    write_csv_result(resDic, vf_virDic, gene_virDic, outdir, dtbasename)
    write_summary_result(resDic, vf_virDic, gene_virDic, outdir, dtbasename, sampleID)


def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='parse__detection - Version ' + version())
    parser.add_argument('-s', '--sampleID', dest="sampleID", default='27841', help='Sample ID')
    parser.add_argument('-sf','--sampleFile', dest="sampleFile", default='', help='Sample file')
    parser.add_argument('-d', '--dtbase', dest="dtbase", default='vir', help="Database name")
    parser.add_argument('-wd','--wkDir',      dest="wkDir", default='/home/bacteriologie/ariba/anses/mcr-3', help="Working directory")
    parser.add_argument('-st','--settingFile',dest="settingFile",default='/usr/local/readmapper-v0.1/setting.txt', help="Setting file")
    parser.add_argument('-v', '--verbose',    dest="verbose",    default="0",help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version',    action='version', version='parse_detection-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
