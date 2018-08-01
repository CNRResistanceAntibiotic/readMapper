#!/usr/bin/python3
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from readmapper.prepare_mapping import read_sample_file, read_setting_file
from readmapper.parser_readmapper.utils_parser import gunzip_file, read_fasta_file, translate_dna, load_arm_db
import pandas as pd


def load_arm_res(tsv_file, gen_file, seq_file):
    if os.path.exists(gen_file):
        gen_file = gunzip_file(gen_file)
    else:
        gen_file = os.path.splitext(gen_file)[0]
    genDic = read_fasta_file(gen_file)

    if os.path.exists(seq_file):
        seq_file = gunzip_file(seq_file)
    else:
        seq_file = os.path.splitext(seq_file)[0]
    seq_dic = read_fasta_file(seq_file)

    res_dic = {}
    with open(tsv_file) as f:
        header = ""
        prot_rec = ""
        dna_rec = ""
        for n, line in enumerate(f):
            if n == 0:
                header = line.strip().split('\t')
            else:
                dt_dic = dict(zip(header, line[:-1].split('\t')))
                ref_name = dt_dic['ref_name']
                # gene = dtDic['gene']  # 1 gene 0 non-coding

                found = '0'
                for key in genDic.keys():
                    if dt_dic['ctg'] in key:
                        found = '1'
                        dna_rec = genDic[key]
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
                    print('Not found:', dt_dic['ctg'])
                    print(dt_dic)

                ref_len = float(dt_dic['ref_len'])
                ref_base_assembled = int(dt_dic['ref_base_assembled'])
                cov = (ref_base_assembled / ref_len) * 100
                dt_dic['pc_coverage'] = '{0}'.format(round(cov, 2))

                ref_nt = dt_dic['ref_nt']
                ctg_nt = dt_dic['ctg_nt']
                nucleotide_change = '{0}->{1}'.format(ref_nt, ctg_nt)
                if nucleotide_change == '.->.':
                    nucleotide_change = '.'
                smtls_nts = dt_dic['smtls_nts']
                smtls_nts_depth = dt_dic['smtls_nts_depth']
                smtls_total_depth = dt_dic['smtls_total_depth']

                sub_depth = parse_nt_info(ref_nt, ctg_nt, smtls_nts, smtls_nts_depth, smtls_total_depth)

                if dt_dic['ref_ctg_effect'] != '.':
                    if dt_dic['known_var_change'] == dt_dic['ref_ctg_change']:
                        known_change = dt_dic['ref_ctg_change']
                        unknown_change = '.'
                    else:
                        known_change = '.'
                        unknown_change = dt_dic['ref_ctg_change']
                else:
                    known_change, unknown_change = '.', '.'

                if ref_name in res_dic.keys():
                    indexes = list(res_dic[ref_name]['mutations'].keys())
                    indexes.sort()
                    index = indexes[-1] + 1
                    if dt_dic['var_type'] == 'HET':
                        res_dic[ref_name]['Warning_HET'] = 'warning:heterozygote'
                    res_dic[ref_name]['mutations'][index] = {'unknown_change': unknown_change,
                                                             'ref_ctg_change': dt_dic['ref_ctg_change'],
                                                             'known_change': known_change,
                                                             'searched_change': dt_dic['known_var_change'],
                                                             'has_known_var': dt_dic['has_known_var'],
                                                             'variant_type': dt_dic['var_type'],
                                                             'variant_seq_type': dt_dic['var_seq_type'],
                                                             'ref_ctg_effect': dt_dic['ref_ctg_effect'],
                                                             'summary_substitution_depth': sub_depth,
                                                             # 'ref_start':dtDic['ref_start'],
                                                             # 'ref_end':dtDic['ref_end'],
                                                             # 'ctg_nt':dtDic['ctg_nt']
                                                             'nucleotide_change': nucleotide_change,
                                                             # 'ctg_start':dtDic['ctg_start'],
                                                             # 'ctg_end':dtDic['ctg_end'],
                                                             'detected_nucleotide': dt_dic['smtls_nts'],
                                                             'depth_by_detected_nucleotide': dt_dic['smtls_nts_depth'],
                                                             'total_depth_by_nucleotide': dt_dic['smtls_total_depth']
                                                             }
                else:
                    res_dic[ref_name] = {'pc_identity': dt_dic['pc_ident'], 'pc_coverage': dt_dic['pc_coverage'],
                                         'gene': dt_dic['gene'],
                                         'var_only': dt_dic['var_only'],
                                         'flag': dt_dic['flag'],
                                         'ctg_name': dt_dic['ctg'],
                                         # 'ctg_len':dtDic['ctg_len'],
                                         'mean_depth': dt_dic['ctg_cov'],
                                         'known_var': dt_dic['known_var'],
                                         'dna_sequence': dna_rec,
                                         'prot_sequence': prot_rec,
                                         'mutations': {
                                             n: {'unknown_change': unknown_change,
                                                 'ref_ctg_change': dt_dic['ref_ctg_change'],
                                                 'known_change': known_change,
                                                 'the_searched_change': dt_dic['known_var_change'],
                                                 'has_known_var': dt_dic['has_known_var'],
                                                 'variant_type': dt_dic['var_type'],
                                                 'variant_seq_type': dt_dic['var_seq_type'],
                                                 'ref_ctg_effect': dt_dic['ref_ctg_effect'],
                                                 'summary_substitution_depth': sub_depth,
                                                 # 'ref_start':dtDic['ref_start'],
                                                 # 'ref_end':dtDic['ref_end'],
                                                 # 'ref_nt':dtDic['ref_nt'],
                                                 # 'ctg_start':dtDic['ctg_start'],
                                                 # 'ctg_end':dtDic['ctg_end'],
                                                 # 'ctg_nt':dtDic['ctg_nt'],
                                                 'nucleotide_change': nucleotide_change,
                                                 'detected_nucleotide': dt_dic['smtls_nts'],
                                                 'depth_by_detected_nucleotide': dt_dic['smtls_nts_depth'],
                                                 'total_depth_by_nucleotide': dt_dic['smtls_total_depth']
                                                 }}}

                    if dt_dic['var_type'] == 'HET':
                        res_dic[ref_name]['Warning_HET'] = 'warning:heterozygote'

                    if float(dt_dic['pc_coverage']) < 95.0:
                        res_dic[ref_name]['Warning_COV'] = 'warning:truncation'

    return res_dic


def parse_nt_info(ref_nt, ctg_nt, smtls_nts, smtls_nts_depth, smtls_total_depth):
    txt = ''
    for i in range(0, len(smtls_total_depth.split(';'))):
        nt = smtls_nts.split(';')[i]
        depth = smtls_nts_depth.split(';')[i]
        total_depth = smtls_total_depth.split(';')[i]

        if ',' not in nt:
            # perc = int(depth) / float(total_depth)
            txt = txt + '|{0}[{1}/{2}]'.format(nt, depth, total_depth)
        else:
            for j in range(0, len(nt.split(','))):
                n = nt.split(',')[j]
                d = depth.split(',')[j]
                # perc = int(d) / float(total_depth)
                if j == 0:
                    txt = txt + '|{0}[{1}/{2}]'.format(n, d, total_depth)
                else:
                    txt = txt + ',{0}[{1}/{2}]'.format(n, d, total_depth)
    txt = txt[1:]

    txt = '{0}->{1} {2}'.format(ref_nt, ctg_nt, txt)
    if txt == '.->. .[./.]':
        txt = '.'

    return txt


def filter_results(res_dic, species, resu_file, passcov=80, passid=80):
    with open(resu_file, 'a') as f:
        del_keys = []
        del_n_muts = []
        del_no_muts = []
        for key in res_dic.keys():

            # filter target coverage < passcov
            if float(res_dic[key]['pc_coverage']) < passcov:
                del_keys.append(key)
                f.write('Filter pc_coverage < {0}\n:'.format(passcov))
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
                f.write('Filter pc_identity < {0}\n:'.format(passid))
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

            # filter
            if res_dic[key]['var_only'] == '1':
                # filter SNP detection in wrong taxonomy
                del_keys = taxon_snp(del_keys, key, species)
                f.write('Filter wrong taxonomy\n:')
                f.write('{0}\n{1}\n\n'.format(key, res_dic[key]))
                # filter SNP search with no SNP
                del_no_muts, del_keys = no_change(del_no_muts, del_keys, res_dic, key)

            # filter nts SNP if not searched
            del_n_muts = unknown_synonymous(del_n_muts, res_dic, key)

        # Deletion SNP in records
        for key, n in del_n_muts:
            f.write('Filter nts SNP if not searched\n:')
            f.write('{0}\n{1}\n{2}\n\n'.format(key, res_dic[key], res_dic[key]['mutations'][n]))
            del res_dic[key]['mutations'][n]

        for key, n in del_no_muts:
            try:
                del res_dic[key]['mutations'][n]
                if key != '':
                    f.write('Filter SNP search with no SNP\n:')
                    f.write('{0}\n{1}\n{2}\n\n'.format(key, res_dic[key], res_dic[key]['mutations'][n]))
            except Exception as e:
                # exception to caracterise
                a = 1
                # print(e)

    # Deletion of records
    for key in list(set(del_keys)):
        del res_dic[key]
    return res_dic


def unknown_synonymous(del_n_muts, res_dic, key):
    for n in res_dic[key]['mutations'].keys():
        if res_dic[key]['mutations'][n]['variant_seq_type'] == 'n' and \
                res_dic[key]['mutations'][n]['known_change'] == '.':
            del_n_muts.append([key, n])
    return del_n_muts


def no_change(del_no_muts, del_keys, res_dic, key):
    no_change_found = True
    for n in res_dic[key]['mutations'].keys():
        if res_dic[key]['mutations'][n]['known_change'] == '.' \
                and res_dic[key]['mutations'][n]['unknown_change'] == '.':
            del_no_muts.append([key, n])
        else:
            no_change_found = False
    if no_change_found:
        del_keys.append(key)

    return del_no_muts, del_keys


def taxon_snp(del_keys, key, species, sep=','):
    pattern = re.compile('[0-9]+[.][0-9]+[_]{2}.+[_]{2}([A-Za-z0-9, _]+)[_]')
    entero_bacteriaceae = ['citrobacter freundii', 'enterobacter cloacae', 'enterobacter aerogenes',
                           'enterobacter asburiae', 'enterobacter hormaechei', 'escherichia coli',
                           'hafnia alvei', 'klebsiella pneumoniae', 'salmonella enteritidis']
    enterobacter_cloacae_complex = ['enterobacter cloacae', 'enterobacter asburiae', 'enterobacter hormaechei',
                                    'enterobacter kobei', 'enterobacter mori', 'enterobacter ludwigii',
                                    'enterobacter xiangfangensis']

    match = pattern.match(key)
    if match:
        taxons = match.group(1).split(sep)
        for n, taxon in enumerate(taxons):
            taxons[n] = taxon.lower().replace('_', ' ')
        if 'enterobacteriaceae' in taxons:
            taxons = taxons + entero_bacteriaceae
        if 'enterobacter cloacae complex' in taxons:
            taxons = taxons + enterobacter_cloacae_complex

        if species.lower() not in taxons:
            del_keys.append(key)

    return del_keys


def check_allele(res_dic, dt_base_file):
    armDic = load_arm_db(dt_base_file)

    for res_key in res_dic.keys():
        arm_key_1 = res_key.split('__')[0]
        data = res_dic[res_key]

        res_seq = str(data['dna_sequence'].seq)
        arm_seq = str(armDic[arm_key_1]['dna_sequence'])

        if res_seq == arm_seq:
            res_dic = add_res(res_dic, res_key, armDic, arm_key_1)
        else:
            found = True
            if data['gene'] == '1':
                res_seq = translate_dna(data['dna_sequence'].seq)
                arm_seq = translate_dna(Seq(armDic[arm_key_1]['dna_sequence']))
                if str(res_seq) == str(arm_seq) and res_seq != '':
                    res_dic = update_res(res_dic, res_key, armDic, arm_key_1)
                elif res_seq != '':
                    for arm_key_2 in armDic.keys():
                        arm_seq = translate_dna(Seq(armDic[arm_key_2]['dna_sequence']))
                        if str(arm_seq) == str(res_seq):
                            found = False
                            res_dic = update_res(res_dic, res_key, armDic, arm_key_2)
                            break
                if found:
                    res_dic = add_res(res_dic, res_key, armDic, arm_key_1)

            else:
                res_seq = data['dna_sequence'].seq
                for arm_key_2 in armDic.keys():
                    arm_seq = armDic[arm_key_2]['dna_sequence']
                    if str(arm_seq) == str(res_seq):
                        found = False
                        res_dic = update_res(res_dic, res_key, armDic, arm_key_2)
                        break
                if found:
                    res_dic = add_res(res_dic, res_key, armDic, arm_key_1)

    return res_dic


def add_res(res_dic, res_key, arm_dic, arm_key):
    res_dic[res_key]['designation'] = arm_dic[arm_key]['entry_name']
    res_dic[res_key]['keyDB'] = arm_key
    res_dic[res_key]['comment'] = arm_dic[arm_key]['comment']
    res_dic[res_key]['alternative_designation'] = arm_dic[arm_key]['alternative_names']
    res_dic[res_key]['searched_dna_changes'] = arm_dic[arm_key]['dna_snp']
    res_dic[res_key]['searched_prot_changes'] = arm_dic[arm_key]['prot_snp']
    res_dic[res_key]['ATB_groups'] = arm_dic[arm_key]['function_grp_names']
    res_dic[res_key]['mechanism_groups'] = arm_dic[arm_key]['mechanism_names']
    res_dic[res_key]['related_groups'] = arm_dic[arm_key]['cluster90_grp_name']
    return res_dic


def update_res(res_dic, res_key, arm_dic, arm_key):
    res_dic[res_key]['designation'] = arm_dic[arm_key]['entry_name']
    res_dic[res_key]['keyDB'] = arm_key
    res_dic[res_key]['comment'] = arm_dic[arm_key]['comment']
    res_dic[res_key]['alternative_designation'] = arm_dic[arm_key]['alternative_names']
    res_dic[res_key]['searched_dna_changes'] = arm_dic[arm_key]['dna_snp']
    res_dic[res_key]['searched_prot_changes'] = arm_dic[arm_key]['prot_snp']
    res_dic[res_key]['ATB_groups'] = arm_dic[arm_key]['function_grp_names']
    res_dic[res_key]['mechanism_groups'] = arm_dic[arm_key]['mechanism_names']
    res_dic[res_key]['related_groups'] = arm_dic[arm_key]['cluster90_grp_name']
    res_dic[res_key]['pc_identity'] = '100.00'
    res_dic[res_key]['pc_coverage'] = '100.00'
    res_dic[res_key]['pc_coverage'] = '100.00'
    res_dic[res_key]['variant_seq_type'] = '.'
    res_dic[res_key]['known_change'] = '.'
    res_dic[res_key]['unknown_change'] = '.'
    res_dic[res_key]['summary_substitution_depth'] = '.'
    res_dic[res_key]['nucleotide_change'] = '.'
    res_dic[res_key]['detected_nucleotide'] = '.'
    res_dic[res_key]['depth_by_detected_nucleotide'] = '.'
    res_dic[res_key]['total_depth_by_nucleotide'] = '.'
    return res_dic


def write_csv_result(res_dic, out_dir, dt_basename):
    nw_res_Dic = {}
    n = 0
    for key in res_dic.keys():
        for key1 in res_dic[key]['mutations'].keys():
            nw_res_Dic[n] = {}
            for key2 in res_dic[key]['mutations'][key1].keys():
                nw_res_Dic[n][key2] = res_dic[key]['mutations'][key1][key2]
            for key3 in res_dic[key].keys():
                if key3 != 'mutations':
                    if key3 == 'dna_sequence':
                        nw_res_Dic[n][key3] = str(res_dic[key][key3].seq)
                    else:
                        nw_res_Dic[n][key3] = res_dic[key][key3]
            n += 1

    arrays_items = ['keyDB', 'designation', 'pc_identity', 'pc_coverage', 'mean_depth', 'variant_seq_type',
                    'known_change', 'unknown_change',
                    # 'summary_substitution_depth',
                    'nucleotide_change', 'detected_nucleotide', 'depth_by_detected_nucleotide',
                    # 'total_depth_by_nucleotide',
                    'ATB_groups', 'mechanism_groups', 'related_groups',
                    'alternative_designation', 'searched_prot_changes', 'searched_dna_changes', 'comment',
                    'dna_sequence', 'prot_sequence']
    res_list = []
    sort_list = list(nw_res_Dic.keys())
    sort_list.sort()
    for i in sort_list:
        line_list = []
        for item in arrays_items:
            if item in ['pc_identity', 'pc_coverage']:
                line_list.append(float(nw_res_Dic[i][item]))
            else:
                if item == 'designation':
                    txt = ''
                    if 'Warning_HET' in nw_res_Dic[i].keys():
                        txt = nw_res_Dic[i]['Warning_HET'] + ','
                    if 'Warning_COV' in nw_res_Dic[i].keys():
                        txt = txt + nw_res_Dic[i]['Warning_COV'] + ','
                    txt = txt + nw_res_Dic[i][item]
                    line_list.append(txt)
                else:
                    line_list.append(nw_res_Dic[i][item])
        res_list.append(line_list)

    df = pd.DataFrame.from_records(res_list, columns=arrays_items)
    df.sort_values(by=['pc_coverage', 'pc_identity', 'keyDB', 'designation', 'known_change', 'unknown_change'],
                   ascending=[False, False, True, True, False, False], inplace=True)

    writer = pd.ExcelWriter(os.path.join(out_dir, 'results_{0}.xlsx'.format(dt_basename)))
    df.to_excel(writer, 'sheet1', index=False)
    writer.save()
    df.to_csv(os.path.join(out_dir, 'results_{0}.csv'.format(dt_basename)), sep='\t', index=False)


def write_summary_result(res_dic, out_dir, dt_basename, sample_id):
    csv_dic = {}
    xlsDic = {}
    dnaRecords = []
    protRecords = []
    keys = list(res_dic.keys())
    keys.sort()
    for key in keys:
        sequence_ID = '{0}__ctg_X__{1}_locusX__{2}__{3}'.format(
            sample_id, sample_id, res_dic[key]['designation'], res_dic[key]['keyDB'])
        description = 'func:{0},mechanism:{1},id:{2},cov:{3},dep:{4}'.format(res_dic[key]['ATB_groups'],
                                                                             res_dic[key]['mechanism_groups'],
                                                                             res_dic[key]['pc_identity'],
                                                                             res_dic[key]['pc_coverage'],
                                                                             res_dic[key]['mean_depth'])

        col_name = '{0}::{1}::{2}'.format(res_dic[key]['ATB_groups'], res_dic[key]['designation'],
                                          res_dic[key]['keyDB'])

        val_name = ''
        if 'Warning_HET' in res_dic[key].keys():
            val_name = res_dic[key]['Warning_HET'] + ','
        if 'Warning_COV' in res_dic[key].keys():
            val_name = res_dic[key]['Warning_COV'] + ','
        val_name = val_name + 'id:{0},cov:{1},dep:{2}'.format(
            res_dic[key]['pc_identity'], res_dic[key]['pc_coverage'], res_dic[key]['mean_depth'])

        snp, sub = {}, {}
        if res_dic[key]['var_only'] == '1':
            for key2 in res_dic[key]['mutations'].keys():
                if res_dic[key]['mutations'][key2]['variant_seq_type'] == 'p':
                    var_sequence_type = 'pt'
                elif res_dic[key]['mutations'][key2]['variant_seq_type'] == 'n':
                    var_sequence_type = 'nt'
                else:
                    var_sequence_type = '?'

                if res_dic[key]['mutations'][key2]['known_change'] != '.':
                    try:
                        snp[var_sequence_type].append(res_dic[key]['mutations'][key2]['known_change'])
                    except KeyError:
                        snp[var_sequence_type] = [res_dic[key]['mutations'][key2]['known_change']]

                if res_dic[key]['mutations'][key2]['unknown_change'] != '.':
                    try:
                        sub[var_sequence_type].append(res_dic[key]['mutations'][key2]['unknown_change'])
                    except KeyError:
                        sub[var_sequence_type] = [res_dic[key]['mutations'][key2]['unknown_change']]
        snp_txt = ''
        if snp != {}:
            for key3 in snp.keys():
                snp_txt = snp_txt + ',snp:{0}[{1}]'.format('|'.join(snp[key3]), key3)
            description = description + snp_txt
            val_name = val_name + snp_txt

        if snp_txt == '' and res_dic[key]['var_only'] == '1':
            val_name = val_name + ',snp:None'

        sub_txt = ''
        if sub != {}:
            for key3 in sub.keys():
                sub_txt = sub_txt + ',sub:{0}[{1}]'.format('|'.join(sub[key3]), key3)
            description = description + sub_txt
            val_name = val_name + sub_txt

        record = False
        if res_dic[key]['var_only'] == '0':
            record = True
        elif snp != {}:
            record = True
        elif sub != {}:
            record = True

        if record:
            csv_dic[col_name] = val_name

            dna_sequence = SeqRecord(res_dic[key]['dna_sequence'].seq, id=sequence_ID, name=sequence_ID,
                                     description=description)
            dnaRecords.append(dna_sequence)
            if res_dic[key]['prot_sequence'] != '.' and res_dic[key]['prot_sequence'] != '.':
                prot_sequence = SeqRecord(Seq(res_dic[key]['prot_sequence']), id=sequence_ID, name=sequence_ID,
                                          description=description)
                protRecords.append(prot_sequence)

    SeqIO.write(dnaRecords, open(os.path.join(out_dir, 'results_{0}.fna'.format(dt_basename)), 'w'), 'fasta')
    SeqIO.write(protRecords, open(os.path.join(out_dir, 'results_{0}.faa'.format(dt_basename)), 'w'), 'fasta')

    df = pd.DataFrame(csv_dic, index=[sample_id, ])
    df.sort_index(axis=1, inplace=True)
    df.to_csv(os.path.join(out_dir, 'summary_results_{0}.csv'.format(dt_basename)), sep='\t', index=False)

    for key in csv_dic.keys():
        keys = key.split('::')
        try:
            xlsDic[keys[0]]['::'.join(keys[1:])] = csv_dic[key]
        except KeyError:
            xlsDic[keys[0]] = {'::'.join(keys[1:]): csv_dic[key]}

    sheet_names = list(xlsDic.keys())
    sheet_names.sort()
    writer = pd.ExcelWriter(os.path.join(out_dir, 'summary_results_{0}.xlsx'.format(dt_basename)))
    for n, sheet_name in enumerate(sheet_names):
        df = pd.DataFrame(xlsDic[sheet_name], index=[sample_id, ])
        df.sort_index(axis=1, inplace=True)
        df.to_excel(writer, sheet_name, index=False)
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
    sample_dic, sampleList = read_sample_file(sample_file)
    species = sample_dic[sample_id]
    set_species = set_dic[species.lower()]
    dt_basename = '{0}_{1}'.format(set_species[dt_base_type][0], set_species[dt_base_type][1])
    dt_base_file = os.path.join(db_path + "/dbARM/", dt_basename + '.csv')

    tsv_file = os.path.join(out_dir, dt_basename, 'report.tsv')
    gen_file = os.path.join(out_dir, dt_basename, 'assembled_genes.fa.gz')
    seq_file = os.path.join(out_dir, dt_basename, 'assembled_seqs.fa.gz')

    resu_file = os.path.join(out_dir, dt_basename, 'filtered_out_results.txt')

    res_dic = load_arm_res(tsv_file, gen_file, seq_file)

    res_dic = filter_results(res_dic, species, resu_file)

    res_dic = check_allele(res_dic, dt_base_file)

    write_csv_result(res_dic, out_dir, dt_basename)
    write_summary_result(res_dic, out_dir, dt_basename, sample_id)


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='parse__detection - Version ' + version())
    parser.add_argument('-s', '--sampleID', dest="sampleID", default='CNR1979', help='Sample ID')
    parser.add_argument('-sf', '--sampleFile', dest="sampleFile", default='sample.csv', help='Sample file')
    parser.add_argument('-d', '--dtbase', dest="dtbase", default='arm', help="Database name")
    parser.add_argument('-wd', '--wkDir', dest="wkDir", default='', help="Working directory")
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
