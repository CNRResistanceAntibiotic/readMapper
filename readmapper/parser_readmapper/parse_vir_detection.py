#!/usr/bin/python3
import os
import argparse
from readmapper.prepare_mapping import read_sample_file, read_setting_file
from readmapper.parser_readmapper.utils_parser import gunzip_file, read_fasta_file, filter_results, translate_dna
from collections import OrderedDict
import pandas as pd
import copy


def load_vir_db(inp_file):
    print('database used: {0}'.format(inp_file))
    gene_virDic = {}
    vf_virDic = {}
    with open(inp_file, 'r') as f:
        header = ""
        for n, line in enumerate(f):
            line = line.strip()
            if n == 0:
                header = line.strip().split('\t')
            elif line != '' and n > 0:
                line = line.split('\t')[:-1]
                data = dict(zip(header, line))

                vf = data['VF_Accession']
                strain = '{0}__{1}__{2}'.format(data['genus'], data['species'], data['strain'].replace(' ', '_'))

                gene_virDic[data['gene_id']] = data
                try:
                    vf_virDic[vf][strain].append(data['gene_id'])
                except KeyError:
                    if vf not in vf_virDic.keys():
                        vf_virDic[vf] = {strain: [data['gene_id']]}
                    else:
                        vf_virDic[vf][strain] = [data['gene_id']]

    print('\nLoading of {0} done!'.format(inp_file))
    print('Number of virulence genes: {0}'.format(len(gene_virDic.keys())))
    print('Number of virulence factors: {0}'.format(len(vf_virDic)))
    return gene_virDic, vf_virDic


def load_clu_db(cluster_file):
    clu_dic = {}
    with open(cluster_file) as f:
        for line in f:
            line = line.strip().split('\t')

            key = line[0].split('::')[0]
            gene_id = line[1].split('::')[0]

            if len(line) > 2:
                if key not in clu_dic.keys():
                    clu_dic[key] = {gene_id: []}
                if gene_id not in clu_dic[key].keys():
                    clu_dic[key][gene_id] = []

                homologs = line[2].split(',')
                for homolog in homologs:
                    gene_id2 = homolog.split('::')[0]
                    clu_dic[key][gene_id].append(gene_id2)

    return clu_dic


def load_vir_res(tsv_file, gen_file, seq_file):
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
        prot_rec = ""
        header = ""
        for n, line in enumerate(f):
            if n == 0:
                header = line.strip().split('\t')
            else:
                dt_dic = dict(zip(header, line[:-1].split('\t')))
                ref_name = dt_dic['ref_name']
                # gene = dtDic['gene']  # 1 gene 0 non-coding

                found = '0'
                for key in gen_dic.keys():
                    if dt_dic['ctg'] in key:
                        found = '1'
                        dna_rec = gen_dic[key]
                        prot_rec = str(translate_dna(dna_rec.seq))
                        break
                if found == '0':
                    for key in seq_dic.keys():
                        if dt_dic['ref_name'] in key:
                            found = '1'
                            dna_rec = seq_dic[key]
                            prot_rec = '.'
                            break

                if found == '0':
                    print('\nSequence not found:', dt_dic['ctg'])
                    print(dt_dic)
                    exit(1)

                ref_len = float(dt_dic['ref_len'])
                ref_base_assembled = int(dt_dic['ref_base_assembled'])
                cov = (ref_base_assembled / ref_len) * 100
                dt_dic['pc_coverage'] = '{0}'.format(round(cov, 2))

                update = '1'
                if ref_name in res_dic.keys():
                    if cov > float(res_dic[ref_name]['pc_coverage']):
                        update = '1'
                    else:
                        update = '0'

                if update == '1':
                    key = ref_name.split('__')[0]
                    res_dic[key] = {'pc_identity': dt_dic['pc_ident'], 'pc_coverage': dt_dic['pc_coverage'],
                                    'gene': dt_dic['gene'], 'flag': dt_dic['flag'], 'ctg_name': dt_dic['ctg'],
                                    # 'ctg_len':dtDic['ctg_len'],
                                    'ref_name': dt_dic['ref_name'],
                                    'mean_depth': dt_dic['ctg_cov'], 'dna_sequence': dna_rec, 'prot_sequence': prot_rec}

    return res_dic


def score_results(res_dic, gene_vir_dic, vf_vir_dic, clu_dic):
    """
    This function filter the score results of the virulence detection
    :param res_dic: The dict of results of virulence detection
    :param gene_vir_dic: The dict of gene by virulence reference
    :param vf_vir_dic: The dict of virulence reference by strain
    :param clu_dic:
    :return: The dict of score filtered and updated
    """
    vf_list = []
    for key in res_dic.keys():
        vf_list.append(gene_vir_dic[key]['VF_Accession'])
    vf_list = list(set(vf_list))

    score_dic = {}
    for vf in vf_list:
        strain_list = list(vf_vir_dic[vf].keys())
        strain_list.sort()
        for key in strain_list:
            total_gene = 0
            positive_gene = 0
            gene_list = vf_vir_dic[vf][key]
            gene_list.sort()
            gene_dic = {}
            for gene_id in gene_list:
                total_gene += 1
                try:
                    pc_coverage = res_dic[gene_id]['pc_coverage']
                    pc_identity = res_dic[gene_id]['pc_identity']
                    positive_gene += 1
                except KeyError:
                    pc_coverage = '.'
                    pc_identity = '.'

                gene_dic[gene_id] = {'pc_coverage': pc_coverage, 'pc_identity': pc_identity}

            if positive_gene > 0:
                if vf in score_dic.keys():
                    score_dic[vf][key] = {'gene_dic': gene_dic}
                    score_dic[vf][key]['score'] = {'positive_gene': positive_gene, 'total_gene': total_gene}
                else:
                    score_dic[vf] = {key: {'gene_dic': gene_dic}}
                    score_dic[vf][key]['score'] = {'positive_gene': positive_gene, 'total_gene': total_gene}

    # Display the result of virulence dectection
    print_score_dic(score_dic, vf_vir_dic, gene_vir_dic)

    print('\nUpdated results:')

    updated_dic = copy.deepcopy(score_dic)
    vf_list = list(score_dic.keys())
    vf_list.sort()

    max_pos_vf = {}
    max_t_gene_vf = {}

    vf2 = ""
    key_strain2 = ""
    for vf in vf_list:
        strain_list = list(score_dic[vf].keys())
        strain_list.sort()

        for strain in strain_list:
            positive_gene = score_dic[vf][strain]['score']['positive_gene']
            total_gene = score_dic[vf][strain]['score']['total_gene']
            if vf not in max_pos_vf.keys():
                max_pos_vf[vf] = positive_gene
                max_t_gene_vf[vf] = total_gene
            else:
                if positive_gene > max_pos_vf[vf]:
                    max_pos_vf[vf] = positive_gene
                if total_gene > max_t_gene_vf[vf]:
                    max_t_gene_vf[vf] = total_gene
            gene_list = vf_vir_dic[vf][strain]
            gene_list.sort()

            txt = ''
            vf_name = ''

            for gene_id in gene_list:
                vf_name = gene_vir_dic[gene_id]['VF_Name']
                gene_name = gene_vir_dic[gene_id]['gene_name']

                try:
                    homolog_gene_ids = clu_dic[vf][gene_id]
                except KeyError:
                    homolog_gene_ids = []

                homolog_gene_id = ''
                vf2_name = ''
                gene2_name = ''
                genus2 = ''
                species2 = ''
                strain2 = ''

                for homolog_gene_id in homolog_gene_ids:
                    if homolog_gene_id in res_dic.keys():
                        vf2 = gene_vir_dic[homolog_gene_id]['VF_Accession']
                        vf2_name = gene_vir_dic[homolog_gene_id]['VF_Name']
                        gene2_name = gene_vir_dic[homolog_gene_id]['gene_name']
                        genus2 = gene_vir_dic[homolog_gene_id]['genus']
                        species2 = gene_vir_dic[homolog_gene_id]['species']
                        strain2 = gene_vir_dic[homolog_gene_id]['strain'].replace(' ', '_')
                        key_strain2 = '{0}__{1}__{2}'.format(genus2, species2, strain2)
                        break
                    else:
                        homolog_gene_id = ''

                pc_coverage = str(score_dic[vf][strain]['gene_dic'][gene_id]['pc_coverage'])
                pc_identity = str(score_dic[vf][strain]['gene_dic'][gene_id]['pc_identity'])

                if pc_coverage == '.' and pc_identity == '.' and homolog_gene_id != '':
                    updated_dic[vf][strain]['gene_dic'][gene_id]['pc_coverage'] = \
                        score_dic[vf2][key_strain2]['gene_dic'][homolog_gene_id]['pc_coverage']
                    updated_dic[vf][strain]['gene_dic'][gene_id]['pc_identity'] = \
                        score_dic[vf2][key_strain2]['gene_dic'][homolog_gene_id]['pc_identity']
                    updated_dic[vf][strain]['gene_dic'][gene_id]['homolog'] = \
                        'VF_id:{0},VF_name:{1},genus:{2},species:{3},strain:{4},gene_id:{5},gene_name:{6}' \
                            .format(vf2, vf2_name, genus2, species2, strain2, homolog_gene_id, gene2_name)
                    pc_coverage = str(updated_dic[vf][strain]['gene_dic'][gene_id]['pc_coverage'])
                    pc_identity = str(updated_dic[vf][strain]['gene_dic'][gene_id]['pc_identity'])
                    updated_dic[vf][strain]['score']['positive_gene'] += 1
                    gene_name = '{0}={1}'.format(gene_name, gene2_name)

                txt = txt + ' gene: {0} [{1}/{2}]'.format(gene_name, pc_identity, pc_coverage)

            txt = txt[1:]
            print('UPDATED : VF_id: {0}\t{1}\tVF_name: {2}\tScore: {3}/{4}\t{5}'.format(vf, strain, vf_name,
                                                                                        positive_gene,
                                                                                        total_gene, txt))

    rm_dic = {}
    kp_dic = {}

    # get strain to remove of final results
    for vf in updated_dic.keys():
        max_pos = max_pos_vf[vf]
        for strain in updated_dic[vf].keys():
            positive_gene = score_dic[vf][strain]['score']['positive_gene']
            if positive_gene < max_pos:
                if vf not in rm_dic.keys():
                    rm_dic[vf] = [strain]
                else:
                    rm_dic[vf].append(strain)
            else:
                if vf not in kp_dic.keys():
                    kp_dic[vf] = [strain]
                else:
                    kp_dic[vf].append(strain)
    # get strain to remove of final results
    for vf in kp_dic.keys():
        if len(kp_dic[vf]) > 1:
            max_t_gene = max_t_gene_vf[vf]
            for strain in kp_dic[vf]:
                total_gene = score_dic[vf][strain]['score']['total_gene']
                if total_gene < max_t_gene:
                    if vf not in rm_dic.keys():
                        rm_dic[vf] = [strain]
                    else:
                        rm_dic[vf].append(strain)

    # remove the strain unwanted
    for vf in rm_dic.keys():
        for strain in rm_dic[vf]:
            del updated_dic[vf][strain]
            if updated_dic[vf] == {}:
                print('DELETE VF {0}'.format(vf))
                del updated_dic[vf]

    print('\nPurged results:')
    print_score_dic(updated_dic, vf_vir_dic, gene_vir_dic)

    print("\nFiltering finish !\n")

    return updated_dic


def print_score_dic(score_dic, vf_vir_dic, gene_vir_dic):
    """
    This function prints the results of the virulence detection
    :param score_dic: the dict with score result records
    :param vf_vir_dic: the dict with virulence reference by species
    :param gene_vir_dic: the dict of gene ref by virulence reference
    :return: nothing
    """
    vf_list = list(score_dic.keys())
    vf_list.sort()
    for i, vf in enumerate(vf_list):
        strain_list = list(score_dic[vf].keys())
        strain_list.sort()
        for strain in strain_list:
            positive_gene = score_dic[vf][strain]['score']['positive_gene']
            total_gene = score_dic[vf][strain]['score']['total_gene']
            score = '{0}/{1}'.format(positive_gene, total_gene)
            gene_list = vf_vir_dic[vf][strain]
            gene_list.sort()

            txt = ''
            vf_name = ''

            for gene_id in gene_list:
                vf_name = gene_vir_dic[gene_id]['VF_Name']
                gene_name = gene_vir_dic[gene_id]['gene_name']

                pc_coverage = str(score_dic[vf][strain]['gene_dic'][gene_id]['pc_coverage'])
                pc_identity = str(score_dic[vf][strain]['gene_dic'][gene_id]['pc_identity'])

                if 'homolog' in score_dic[vf][strain]['gene_dic'][gene_id].keys():
                    homolog = score_dic[vf][strain]['gene_dic'][gene_id]['homolog']
                    txt = txt + 'gene: {0} [{1}: {2}/{3}]\t'.format(gene_name, homolog, pc_identity, pc_coverage)
                else:
                    txt = txt + 'gene: {0} [{1}/{2}]\t'.format(gene_name, pc_identity, pc_coverage)

            txt = txt[:-1]
            print('RESULTS : {0}\t{1}\t{2}\t{3}\t{4}'.format(vf, strain, vf_name, score, txt))


def write_csv_result(res_dic, vf_vir_dic, gene_vir_dic, out_dir, dt_basename):
    with open(os.path.join(out_dir, 'results_{0}.csv'.format(dt_basename)), 'w') as f:
        f.write('VF_Accession\tStrain\tVF_Name\tScore\tGenes')
        vf_list = list(res_dic.keys())
        vf_list.sort()
        for vf in vf_list:
            strain_list = list(res_dic[vf].keys())
            strain_list.sort()
            for strain in strain_list:
                positive_gene = res_dic[vf][strain]['score']['positive_gene']
                total_gene = res_dic[vf][strain]['score']['total_gene']
                score = '{0}/{0}'.format(positive_gene, total_gene)

                gene_list = vf_vir_dic[vf][strain]
                gene_list.sort()

                txt = ''
                vf_name = ''

                for gene_id in gene_list:
                    vf_name = gene_vir_dic[gene_id]['VF_Name']
                    gene_name = gene_vir_dic[gene_id]['gene_name']

                    pc_coverage = str(res_dic[vf][strain]['gene_dic'][gene_id]['pc_coverage'])
                    pc_identity = str(res_dic[vf][strain]['gene_dic'][gene_id]['pc_identity'])

                    if 'homolog' in res_dic[vf][strain]['gene_dic'][gene_id].keys():
                        homolog = res_dic[vf][strain]['gene_dic'][gene_id]['homolog']
                        txt = txt + 'gene:{0},homolog:{1},pc_id:{2},pc_cv:{3}\t'.format(
                            gene_name, homolog, pc_identity, pc_coverage)
                    else:
                        txt = txt + 'gene:{0},pc_id:{1},pc_cv:{2}\t'.format(gene_name, pc_identity, pc_coverage)
                txt = '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(vf, strain, vf_name, score, txt)
                f.write(txt)


def write_summary_result(res_dic, vf_vir_dic, gene_vir_dic, out_dir, dt_basename, sample_id):
    """

    :param res_dic:
    :param vf_vir_dic:
    :param gene_vir_dic:
    :param out_dir:
    :param dt_basename:
    :param sample_id:
    :return:
    """

    writer = pd.ExcelWriter(os.path.join(out_dir, 'summary_results_{0}.xlsx'.format(dt_basename)))
    vf_list = list(res_dic.keys())
    vf_list.sort()
    vf_sheet_hash = {}

    for vf in vf_list:
        strain_list = list(res_dic[vf].keys())
        strain_list.sort()
        for strain in strain_list:
            # positive_gene = res_dic[vf][strain]['score']['positive_gene']
            # total_gene = res_dic[vf][strain]['score']['total_gene']
            # score = '%i/%i' % (positive_gene, total_gene)

            gene_list = vf_vir_dic[vf][strain]
            gene_list.sort()
            data = OrderedDict()
            # data['sample_id'] = sample_id
            data['Strain'] = strain
            vf_name = ''

            for gene_id in gene_list:
                vf_name = gene_vir_dic[gene_id]['VF_Name']
                gene_name = gene_vir_dic[gene_id]['gene_name']

                pc_identity = str(res_dic[vf][strain]['gene_dic'][gene_id]['pc_identity'])
                if 'homolog' in res_dic[vf][strain]['gene_dic'][gene_id].keys():
                    homologs = res_dic[vf][strain]['gene_dic'][gene_id]['homolog'].split(',')

                    h_genus = ''
                    h_species = ''
                    h_strain = ''

                    for item in homologs:
                        if 'species' in item:
                            h_species = item.split(':')[1]
                        if 'genus' in item:
                            h_genus = item.split(':')[1]
                        if 'strain' in item:
                            h_strain = item.split(':')[1]
                    h_strain = '{0}__{1}__{2}'.format(h_genus, h_species, h_strain)

                    pc_identity = pc_identity + ' [{0}]'.format(h_strain)

                data[gene_name] = pc_identity

                # push data to common vf_sheet_hash
                if vf_name in vf_sheet_hash:
                    data_hash = vf_sheet_hash.get(vf_name)
                    vf_sheet_hash[vf_name.lower()] = {**data_hash, **data}
                else:
                    vf_sheet_hash[vf_name.lower()] = data

    for vf_name in vf_sheet_hash:
        data = vf_sheet_hash.get(vf_name)
        df = pd.DataFrame(data, index=[sample_id])
        vf_name = vf_name.replace('/', ' ')

        if len(vf_name) >= 31:
            vf_name = vf_name[:30]

        df.to_excel(writer, sheet_name=vf_name, index=True, index_label='sample_id')

    writer.save()


def pre_main(args):
    sample_id = args.sampleID
    sample_file = args.sampleFile
    setting_file = args.settingFile
    dt_base_type = args.dtbase
    wk_dir = args.wkDir
    db_path = args.databasePath

    # execution main
    main(sample_id, sample_file, setting_file, dt_base_type, wk_dir, db_path)


def main(sample_id, sample_file, setting_file, dt_base_type, wk_dir, db_path):
    if sample_file == '':
        sample_file = os.path.join(wk_dir, 'sample.csv')
    if wk_dir == '':
        wk_dir = os.path.dirname(sample_file)
    out_dir = os.path.join(wk_dir, sample_id)

    set_dic = read_setting_file(setting_file)
    sample_dic, sample_list = read_sample_file(sample_file)
    species = sample_dic[sample_id]
    set_species = set_dic[species.lower()]
    dt_basename = '{0}_{1}'.format(set_species[dt_base_type][0], set_species[dt_base_type][1])
    dt_base_file = os.path.join(db_path + "/dbVIR/", '_'.join(dt_basename.split('_')[:2]) + '.csv')
    cluster_file = os.path.join(db_path + "/dbVIR/", '_'.join(dt_basename.split('_')[:2]) + '.clu')

    tsv_file = os.path.join(out_dir, dt_basename, 'report.tsv')
    gen_file = os.path.join(out_dir, dt_basename, 'assembled_genes.fa.gz')
    seq_file = os.path.join(out_dir, dt_basename, 'assembled_seqs.fa.gz')

    res_dic = load_vir_res(tsv_file, gen_file, seq_file)

    print('Number of results: {0}'.format(len(res_dic)))

    gene_vir_dic, vf_vir_dic = load_vir_db(dt_base_file)
    clu_dic = load_clu_db(cluster_file)

    resu_file = os.path.join(out_dir, dt_basename, 'filtered_out_vir_results.txt')
    res_dic = filter_results(res_dic, resu_file, passcov=70, passid=70)

    print('Number of parsed results: {0}'.format(len(res_dic)))

    res_dic = score_results(res_dic, gene_vir_dic, vf_vir_dic, clu_dic)

    write_csv_result(res_dic, vf_vir_dic, gene_vir_dic, out_dir, dt_basename)
    write_summary_result(res_dic, vf_vir_dic, gene_vir_dic, out_dir, dt_basename, sample_id)


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='parse__detection - Version ' + version())
    parser.add_argument('-s', '--sampleID', dest="sampleID", default='27841', help='Sample ID')
    parser.add_argument('-sf', '--sampleFile', dest="sampleFile", default='', help='Sample file')
    parser.add_argument('-d', '--dtbase', dest="dtbase", default='vir', help="Database name")
    parser.add_argument('-wd', '--wkDir', dest="wkDir", default='/home/bacteriologie/ariba/anses/mcr-3',
                        help="Working directory")
    parser.add_argument('-st', '--settingFile', dest="settingFile", default='/usr/local/readmapper-v0.1/setting.txt',
                        help="Setting file")
    parser.add_argument('-db', '--databasePath', dest="databasePath", default='',
                        help="Database directory path")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0",
                        help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='parse_detection-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
