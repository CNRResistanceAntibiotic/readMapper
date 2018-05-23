#!/usr/bin/python3
import os
from Bio import SeqIO


def gunzip_file(gz_file):
    cmd = 'gunzip {0}'.format(gz_file)
    os.system(cmd)
    return os.path.splitext(gz_file)[0]


def read_fasta_file(fa_file):
    faDic = {}
    for rec in SeqIO.parse(open(fa_file), 'fasta'):
        faDic[rec.id] = rec
    return faDic


def filter_results(res_dic, resu_file, passcov=70, passid=65):
    with open(resu_file, 'a') as f:
        del_keys = []
        for key in res_dic.keys():

            # filter target coverage < passcov
            if float(res_dic[key]['pc_coverage']) < passcov:
                del_keys.append(key)
                f.write('Filter pc_coverage < {0}\n:'.format(round(passcov, 2)))
                f.write('{0}\n{1}\n\n'.format(key, res_dic[key]))

                # Alarm with id >= 90 and passid > 40
                if float(res_dic[key]['pc_identity']) >= passid \
                        and float(res_dic[key]['pc_coverage']) > 50:
                    print('\n')
                    print('############################################################')
                    print('###  WARNING PUTATIVE MIS-DETECTION: Low coverage alert  ###')
                    print('############################################################')
                    print('')
                    print('Record:', key, '\tCoverage:', res_dic[key]['pc_coverage'],
                          '\tIdentity:', res_dic[key]['pc_identity'], '\tMean depth:', res_dic[key]['mean_depth'])
                    print('')
                    print('The record will be deleted in the final results')
                    print('')

            # filter id percentage < passid
            if float(res_dic[key]['pc_identity']) < passid:
                del_keys.append(key)
                f.write('Filter pc_identity < {0}\n:'.format(round(passid, 2)))
                f.write('{0}\n{1}\n\n'.format(key, res_dic[key]))
                # Alarm if id >= 40 and cov >=75
                if float(res_dic[key]['pc_identity']) >= 40 \
                        and float(res_dic[key]['pc_coverage']) >= 80:
                    print('\n')
                    print('############################################################')
                    print('###  WARNING PUTATIVE MIS-DETECTION: Low identity alert  ###')
                    print('############################################################')
                    print('')
                    print('Record:', key, '\tCoverage:', res_dic[key]['pc_coverage'],
                          '\tIdentity:', res_dic[key]['pc_identity'], '\tMean depth:', res_dic[key]['mean_depth'])
                    print('')
                    print('The record will be deleted in the final results')
                    print('')

            if float(res_dic[key]['mean_depth']) <= 15 and float(res_dic[key]['pc_identity']) >= passid \
                    and float(res_dic[key]['pc_coverage']) > passcov:
                print('\n')
                print('####################################################################')
                print('###  WARNING PUTATIVE MIS-DETECTION: < 15 sequencing depth alert  ###')
                print('####################################################################')
                print('')
                print('Record:', key, '\tCoverage:', res_dic[key]['pc_coverage'],
                      '\tIdentity:', res_dic[key]['pc_identity'], '\tMean depth:', res_dic[key]['mean_depth'])
                print('')
                print('The record will be kept in the final results')
                print('')

    # Deletion of records
    for key in list(set(del_keys)):
        del res_dic[key]
    return res_dic


def translate_dna(dna_object, table='Bacterial', cds=True):
    prot = ""
    try:
        try:
            prot = dna_object.translate(table=table, cds=cds)
        except Exception as e:
            # print(e)
            prot = ''
        try:
            prot = dna_object.reverse_complement().translate(table=table, cds=cds)
        except Exception as e:
            # print(e)
            prot = ''
    finally:

        return prot


def load_arm_db(inp_file):
    arm_dic = {}
    with open(inp_file, 'r') as f:
        header = ""
        for n, line in enumerate(f):
            line = line.strip()
            if n == 0:
                header = line.split('\t')
            elif line != '' and n > 0:
                line = line.split('\t')
                data = dict(zip(header, line))
                for key in header:
                    if key not in data.keys():
                        data[key] = ''
                key = data['key']
                arm_dic[key] = data

    print('\nLoading of {0} done!'.format(inp_file))
    print('Number of records: {0}'.format(len(arm_dic.keys())))
    return arm_dic
