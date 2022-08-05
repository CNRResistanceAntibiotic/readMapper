#!/usr/bin/python3
import os
from Bio import SeqIO


def gunzip_file(gz_file):
    cmd = f'gunzip {gz_file}'
    os.system(cmd)
    return os.path.splitext(gz_file)[0]


def read_fasta_file(fa_file):
    fa_dic = {}
    for rec in SeqIO.parse(open(fa_file), 'fasta'):
        fa_dic[rec.id] = rec
    return fa_dic


def filter_results(res_dic, resu_file, pass_cov=70, pass_id=65):
    log_message = ""
    with open(resu_file, 'a') as f:
        del_keys = []
        for key_dic in res_dic.keys():
            # filter target coverage < pass_cov
            if float(res_dic[key_dic]['pc_coverage']) < pass_cov:
                del_keys.append(key_dic)
                f.write(f'Filter pc_coverage < {round(pass_cov, 2)}\n:')
                f.write(f'{key_dic}\n{res_dic[key_dic]}\n\n')

                # Alarm with id >= 90 and pass_id > 40
                if float(res_dic[key_dic]['pc_identity']) >= pass_id \
                        and float(res_dic[key_dic]['pc_coverage']) > 50:
                    log_message = log_message + '\n############################################################'
                    log_message = log_message + '###  WARNING PUTATIVE MIS-DETECTION: Low coverage alert  ###'
                    log_message = log_message + '############################################################\n'
                    log_message = log_message + f'Record: {key_dic}\tCoverage: {res_dic[key_dic]["pc_coverage"]}\t' \
                                                f'Identity: {res_dic[key_dic]["pc_identity"]}\t' \
                                                f'Mean depth: {res_dic[key_dic]["mean_depth"]}'
                    log_message = log_message + '\nThe record will be deleted in the final results\n'

            # filter id percentage < pass_id
            if float(res_dic[key_dic]['pc_identity']) < pass_id:
                del_keys.append(key_dic)
                f.write(f'Filter pc_identity < {round(pass_id, 2)}\n:')
                f.write('{0}\n{1}\n\n'.format(key_dic, res_dic[key_dic]))
                # Alarm if id >= 40 and cov >=75
                if float(res_dic[key_dic]['pc_identity']) >= 40 \
                        and float(res_dic[key_dic]['pc_coverage']) >= 80:
                    log_message = log_message + '\n############################################################'
                    log_message = log_message + '###  WARNING PUTATIVE MIS-DETECTION: Low identity alert  ###'
                    log_message = log_message + '############################################################\n'
                    log_message = log_message + 'Record: {0}\tCoverage: {1}\tIdentity: {2}\tMean depth: {3}' \
                        .format(key_dic, res_dic[key_dic]['pc_coverage'], res_dic[key_dic]['pc_identity'],
                                res_dic[key_dic]['mean_depth'])
                    log_message = log_message + '\nThe record will be deleted in the final results\n'

            if float(res_dic[key_dic]['mean_depth']) <= 15 and float(res_dic[key_dic]['pc_identity']) >= pass_id \
                    and float(res_dic[key_dic]['pc_coverage']) > pass_cov:
                log_message = log_message + '\n####################################################################'
                log_message = log_message + '###  WARNING PUTATIVE MIS-DETECTION: < 15 sequencing depth alert  ###'
                log_message = log_message + '####################################################################\n'
                log_message = log_message + 'Record: {0}\tCoverage: {1}\tIdentity: {2}\tMean depth: {3}' \
                    .format(key_dic, res_dic[key_dic]['pc_coverage'], res_dic[key_dic]['pc_identity'],
                            res_dic[key_dic]['mean_depth'])
                log_message = log_message + '\nThe record will be kept in the final results\n'

    # Deletion of records
    for key in list(set(del_keys)):
        del res_dic[key]
    return res_dic, log_message


def translate_dna(dna_object, table='Bacterial', cds=True):
    prot = ""
    try:
        try:
            prot = dna_object.translate(table=table, cds=cds)
        except Exception as e:
            # print("Not found in direct strand", dna_object.translate(table=table, cds=cds), flush=True)

            # try in reverse complement strand
            try:
                prot = dna_object.reverse_complement().translate(table=table, cds=cds)
            except Exception as e:
                # print("Not found in direct and reverse complement strand", e, flush=True)
                prot = ''
    finally:
        return prot


def load_vir_arm_db(inp_file):
    log_message = ""
    arm_dic = {}
    with open(inp_file, 'r') as f:
        header = ""
        n = 0
        for line in f.readlines():
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
            n += 1

    log_message = log_message + f"\nLoading of {inp_file} done!\nNumber of records: {len(arm_dic.keys())}\n"
    return arm_dic, log_message
