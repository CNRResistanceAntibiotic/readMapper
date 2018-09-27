#!/usr/bin/python3
import os
import argparse
from Bio import SeqIO
import pandas as pd

from readmapper.prepare_mapping import read_sample_file, read_setting_file
from readmapper.parser_readmapper.utils_parser import gunzip_file, read_fasta_file, filter_results


def load_rep_db(inp_file):
    log_message = ""
    log_message = log_message + "database used: {0}\n".format(inp_file)
    gene_rep_dic = {}
    with open(inp_file, 'r') as f:
        for n, rec in enumerate(SeqIO.parse(f, 'fasta')):
            gene_rep_dic[rec.id] = rec

    log_message = log_message + "\nLoading of {0} done!\n".format(inp_file)
    log_message = log_message + "Number of rep sequence: {0}\n".format(len(gene_rep_dic.keys()))
    return gene_rep_dic, log_message


def load_rep_res(tsv_file, gen_file, seq_file):
    if os.path.exists(gen_file):
        gen_file = gunzip_file(gen_file)
    else:
        gen_file = os.path.splitext(gen_file)[0]
    gen_dic = read_fasta_file(gen_file)

    if os.path.exists(seq_file):
        seq_file = gunzip_file(seq_file)
    else:
        seq_file = os.path.splitext(seq_file)[0]
    seq_dic = read_fasta_file(seq_file)

    res_dic = {}

    with open(tsv_file) as f:
        dna_rec = ""
        header = ""
        for n, line in enumerate(f):
            if n == 0:
                header = line.strip().split('\t')
            else:
                dtDic = dict(zip(header, line[:-1].split('\t')))
                ref_name = dtDic['ref_name']

                found = '0'
                for key in gen_dic.keys():
                    if dtDic['ctg'] in key:
                        found = '1'
                        dna_rec = gen_dic[key]
                        break
                if found == '0':
                    for key in seq_dic.keys():
                        if dtDic['ref_name'] in key:
                            dna_rec = seq_dic[key]
                            break

                ref_len = float(dtDic['ref_len'])
                ref_base_assembled = int(dtDic['ref_base_assembled'])
                cov = (ref_base_assembled / ref_len) * 100
                dtDic['pc_coverage'] = '{0}'.format(round(cov, 2))

                update = '1'
                if ref_name in res_dic.keys():
                    if cov > float(res_dic[ref_name.split('__')[0]]['pc_coverage']):
                        update = '1'
                    else:
                        update = '0'

                if update == '1':
                    key = ref_name.split('__')[0]
                    res_dic[key] = {'pc_identity': dtDic['pc_ident'], 'pc_coverage': dtDic['pc_coverage'],
                                    'gene': dtDic['gene'], 'flag': dtDic['flag'], 'ctg_name': dtDic['ctg'],
                                    # 'ctg_len':dtDic['ctg_len'],
                                    'ref_name': dtDic['ref_name'],
                                    'mean_depth': dtDic['ctg_cov'], 'dna_sequence': dna_rec}

    return res_dic


def write_csv_result(res_dic, sample_id, out_dir, dt_basename):
    with open(os.path.join(out_dir, 'results_{0}.csv'.format(dt_basename)), 'w') as f:
        f.write('sample_id\treplicon\tmean_depth\tpc_coveraget\tpc_identity\tSequence\n')
        rep_list = list(res_dic.keys())
        rep_list.sort()
        for rep_id in rep_list:
            mean_depth = res_dic[rep_id]['mean_depth']
            pc_coverage = res_dic[rep_id]['pc_coverage']
            pc_identity = res_dic[rep_id]['pc_identity']
            dna_sequence = str(res_dic[rep_id]['dna_sequence'].seq)
            txt = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(sample_id, rep_id, mean_depth, pc_coverage, pc_identity,
                                                          dna_sequence)
            f.write(txt)


def write_summary_result(res_dic, out_dir, dt_basename, sample_id):
    writer = pd.ExcelWriter(os.path.join(out_dir, 'summary_results_{0}.xlsx'.format(dt_basename)))
    rep_list = list(res_dic.keys())
    rep_list.sort()
    data = {}
    for rep in rep_list:
        data[rep] = '{0} [id:{1};cv:{2};dp:{3}]'.format(rep, round(float(res_dic[rep]['pc_identity'])),
                                                        round(float(res_dic[rep]['pc_coverage'])),
                                                        round(float(res_dic[rep]['mean_depth'])))
    df = pd.DataFrame(data, index=[sample_id])
    df.to_excel(writer, 'replicons', index=True)
    writer.save()


def pre_main(args):
    sample_id = args.sampleID
    sample_file = args.sampleFile
    setting_file = args.settingFile
    dt_base_type = args.dtbase
    wk_dir = args.wkDir

    # execution main
    main(sample_id, sample_file, setting_file, dt_base_type, wk_dir)


def main(sample_id, sample_file, setting_file, dt_base_type, wk_dir, subgroup):

    log_message = ""

    if sample_file == '':
        sample_file = os.path.join(wk_dir, 'sample.csv')
    if wk_dir == '':
        wk_dir = os.path.dirname(sample_file)
    out_dir = os.path.join(wk_dir, sample_id)

    set_dic = read_setting_file(setting_file)
    sample_dic, sample_list = read_sample_file(sample_file)
    species = sample_dic[sample_id]
    set_species = set_dic[species.lower()]
    dt_basename = '{0}_{1}'.format(*set_species[dt_base_type], subgroup)

    tsv_file = os.path.join(out_dir, dt_basename, 'report.tsv')
    gen_file = os.path.join(out_dir, dt_basename, 'assembled_genes.fa.gz')
    seq_file = os.path.join(out_dir, dt_basename, 'assembled_seqs.fa.gz')

    res_dic = load_rep_res(tsv_file, gen_file, seq_file)
    log_message = log_message + "Number of results: {0}\n".format(len(res_dic))

    resu_file = os.path.join(out_dir, dt_basename, 'filtered_out_rep_results.txt')

    res_dic = filter_results(res_dic, resu_file, 80, 95)
    log_message = log_message + "Number of parsed results: {0}\n".format(len(res_dic))

    write_csv_result(res_dic, sample_id, out_dir, dt_basename)
    write_summary_result(res_dic, out_dir, dt_basename, sample_id)

    return log_message


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='parse_rep_detection - Version ' + version())
    parser.add_argument('-s', '--sampleID', dest="sampleID", help='Sample ID')
    parser.add_argument('-sf', '--sampleFile', dest="sampleFile", default='', help='Sample file')
    parser.add_argument('-d', '--dtbase', dest="dtbase", default='rep', help="Database type")
    parser.add_argument('-wd', '--wkDir', dest="wkDir", default='/home/bacteriologie/ariba/anses/mcr-3',
                        help="Working directory")
    parser.add_argument('-st', '--settingFile', dest="settingFile", default='/usr/local/readmapper-v0.1/setting.txt',
                        help="Setting file")
    parser.add_argument('-db', '--databasePath', dest="databasePath", default='',
                        help="Database directory path")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0",
                        help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='parse_rep_detection-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
