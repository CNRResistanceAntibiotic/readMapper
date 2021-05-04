#!/usr/bin/python3
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from readmapper_src.prepare_mapping import read_sample_file, read_setting_file
from readmapper_src.parser_readmapper.utils_parser import gunzip_file, read_fasta_file, translate_dna, load_vir_arm_db
import pandas as pd


def load_vir_res(tsv_file, gen_file, seq_file):
    log_message = ""
    log_message = log_message + f"{tsv_file}\n{gen_file}\n{seq_file}\n"
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
                    log_message = log_message + f"Not found:{dt_dic['ctg']}\n"
                    log_message = log_message + f"{dt_dic}\n"

                ref_len = float(dt_dic['ref_len'])
                ref_base_assembled = int(dt_dic['ref_base_assembled'])
                cov = (ref_base_assembled / ref_len) * 100
                dt_dic['pc_coverage'] = f'{round(cov, 2)}'

                ref_nt = dt_dic['ref_nt']
                ctg_nt = dt_dic['ctg_nt']
                nucleotide_change = f'{ref_nt}->{ctg_nt}'
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
            # p_depth = int(depth) / float(total_depth)
            txt = txt + f'|{nt}[{depth}/{total_depth}]'
        else:
            for j in range(0, len(nt.split(','))):
                n = nt.split(',')[j]
                d = depth.split(',')[j]
                # p_depth = int(d) / float(total_depth)
                if j == 0:
                    txt = txt + f'|{n}[{d}/{total_depth}]'
                else:
                    txt = txt + f',{n}[{d}/{total_depth}]'
    txt = txt[1:]

    txt = f'{ref_nt}->{ctg_nt} {txt}'
    if txt == '.->. .[./.]':
        txt = '.'

    return txt


def filter_results(res_dic, species, resu_file, pass_cov=80, pass_id=80):
    log_message = ""

    with open(resu_file, 'a') as f:
        del_keys = []
        del_n_muts = []
        del_no_muts = []
        for key in res_dic.keys():

            # filter target coverage < pass_cov
            if float(res_dic[key]['pc_coverage']) < pass_cov:
                del_keys.append(key)
                f.write(f'Filter pc_coverage < {pass_cov}\n:')
                f.write(f'{key}\n{res_dic[key]}\n\n')

                # Alarm with id >= 90 and pass_id > 40
                if float(res_dic[key]['pc_identity']) >= pass_id and float(res_dic[key]['pc_coverage']) > 50:
                    log_message = log_message + "\n############################################################\n"
                    log_message = log_message + "###  WARNING PUTATIVE MIS-DETECTION: Low coverage alert  ###\n"
                    log_message = log_message + "############################################################\n"
                    log_message = log_message + f"\nRecord: {key}\tCoverage:{res_dic[key]['pc_coverage']}\tIdentity:{res_dic[key]['pc_identity']}\tMean depth:{res_dic[key]['mean_depth']}\n"
                    log_message = log_message + "\nThe record will be deleted in the final results\n"

            # filter id percentage < pass_id
            if float(res_dic[key]['pc_identity']) < pass_id:
                del_keys.append(key)
                f.write(f'Filter pc_identity < {pass_id}\n:')
                f.write(f'{key}\n{res_dic[key]}\n\n')

                # Alarm if id >= 40 and cov >=75
                if float(res_dic[key]['pc_identity']) >= 40 and float(res_dic[key]['pc_coverage']) >= 80:
                    log_message = log_message + "\n############################################################\n"
                    log_message = log_message + "###  WARNING PUTATIVE MIS-DETECTION: Low identity alert  ###\n"
                    log_message = log_message + "############################################################\n"
                    log_message = log_message + f"\nRecord: {key}\tCoverage:{res_dic[key]['pc_coverage']}\tIdentity:{res_dic[key]['pc_identity']}\tMean depth:{res_dic[key]['mean_depth']}\n"
                    log_message = log_message + "\nThe record will be deleted in the final results\n"

            if float(res_dic[key]['mean_depth']) <= 15 and float(res_dic[key]['pc_identity']) >= pass_id \
                    and float(res_dic[key]['pc_coverage']) > pass_cov:
                log_message = log_message + "\n####################################################################\n"
                log_message = log_message + "###  WARNING PUTATIVE MIS-DETECTION: < 15 sequencing depth alert  ###\n"
                log_message = log_message + "####################################################################\n"
                log_message = log_message + f"\nRecord: {key}\tCoverage:{res_dic[key]['pc_coverage']}\tIdentity:{res_dic[key]['pc_identity']}\tMean depth:{res_dic[key]['mean_depth']}\n"
                log_message = log_message + "\nThe record will be kept in the final results\n"

            # filter
            if res_dic[key]['var_only'] == '1':
                # filter SNP detection in wrong taxonomy
                del_keys = taxon_snp(del_keys, key, species)
                f.write('Filter wrong taxonomy\n:')
                f.write(f'{key}\n{res_dic[key]}\n\n')
                # filter SNP search with no SNP
                del_no_muts, del_keys = no_change(del_no_muts, del_keys, res_dic, key)

            # filter nts SNP if not searched
            del_n_muts = unknown_synonymous(del_n_muts, res_dic, key)

        # Deletion SNP in records
        for key, n in del_n_muts:
            f.write('Filter nts SNP if not searched\n:')
            f.write(f'{key}\n{res_dic[key]}\n{res_dic[key]["mutations"][n]}\n\n')
            del res_dic[key]['mutations'][n]

        for key, n in del_no_muts:
            try:
                del res_dic[key]['mutations'][n]
                if key != '':
                    f.write('Filter SNP search with no SNP\n:')
                    f.write(f'{key}\n{res_dic[key]}\n{res_dic[key]["mutations"][n]}\n\n')
            except Exception as e:
                log_message = log_message + f"Warning: Problem to deletion SNP. {e}"

    # Deletion of records
    for key in list(set(del_keys)):
        del res_dic[key]

    return res_dic, log_message


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
    vir_dic, log_message = load_vir_arm_db(dt_base_file)
    perfect_nucl_match = perfect_prot_match = 0
    # explore each gene result of ariba
    for res_key in res_dic.keys():
        vir_key_1 = res_key.split('__')[0]
        data = res_dic[res_key]
        res_seq_nucl = str(data['dna_sequence'].seq)
        vir_seq_nucl = str(vir_dic[vir_key_1]['dna_sequence'])
        # test if the ariba result matches perfectly with an existing nucl sequence in database
        if res_seq_nucl == vir_seq_nucl:
            perfect_nucl_match += 1
            res_dic = add_res(res_dic, res_key, vir_dic, vir_key_1)
            continue
        # lot of case where the truly result is the translated prot sequence against the nucl sequence in database
        else:
            found = True
            # case of gene coding (cds) in ariba result
            if data['gene'] == '1':
                res_seq_prot = translate_dna(data['dna_sequence'].seq)
                vir_seq_prot = translate_dna(Seq(vir_dic[vir_key_1]['dna_sequence']))
                # if the ariba result in prot match perfectly with the prot translation of an existing nucl sequence
                # in database
                if str(res_seq_prot) == str(vir_seq_prot) and res_seq_prot != '':
                    res_dic = update_res(res_dic, res_key, vir_dic, vir_key_1)
                # if the match doesnt appear
                elif res_seq_prot != '':
                    # for each existing entry in the database
                    for vir_key_2 in vir_dic.keys():
                        # translate into protein an entry of the nucl sequence in database
                        vir_seq_prot = translate_dna(Seq(vir_dic[vir_key_2]['dna_sequence']))
                        # if the ariba result in prot match perfectly with the prot translation of an existing nucl
                        # sequence in database
                        if str(vir_seq_prot) == str(res_seq_prot):
                            found = False
                            res_dic = update_res(res_dic, res_key, vir_dic, vir_key_2)
                            break
                # if any match founded. the ariba result will be conserved
                if found:
                    res_dic = add_res(res_dic, res_key, vir_dic, vir_key_1)
            # case of gene no coding (dna) in ariba result
            else:
                res_seq_prot = data['dna_sequence'].seq
                for vir_key_2 in vir_dic.keys():
                    vir_seq_prot = vir_dic[vir_key_2]['dna_sequence']
                    if str(vir_seq_prot) == str(res_seq_prot):
                        found = False
                        res_dic = update_res(res_dic, res_key, vir_dic, vir_key_2)
                        break
                if found:
                    res_dic = add_res(res_dic, res_key, vir_dic, vir_key_1)

    log_message = log_message + f"Number of results: {len(res_dic)}\n"
    if res_dic:
        log_message = log_message + f"Number of perfect nucleotide match of ariba result to" \
                                    f" {os.path.basename(dt_base_file)} database: {perfect_nucl_match} -" \
                                    f" {round((perfect_nucl_match / len(res_dic) * 100), 2)}% of total results\n"
        log_message = log_message + f"Number of perfect protein match of ariba result to {os.path.basename(dt_base_file)}" \
                                    f" database: {perfect_prot_match} -" \
                                    f" {round((perfect_prot_match / len(res_dic) * 100), 2)}% of total results\n"

    return res_dic, log_message


def add_res(res_dic, res_key, vir_dic, vir_key):
    res_dic[res_key]['designation'] = vir_dic[vir_key]['entry_name']
    res_dic[res_key]['keyDB'] = vir_key
    res_dic[res_key]['comment'] = vir_dic[vir_key]['comment']
    res_dic[res_key]['alternative_designation'] = vir_dic[vir_key]['alternative_names']
    res_dic[res_key]['searched_dna_changes'] = vir_dic[vir_key]['dna_snp']
    res_dic[res_key]['searched_prot_changes'] = vir_dic[vir_key]['prot_snp']
    res_dic[res_key]['ATB_groups'] = vir_dic[vir_key]['antibiotic_family']
    res_dic[res_key]['mechanism_groups'] = vir_dic[vir_key]['mechanism_names']
    return res_dic


def update_res(res_dic, res_key, vir_dic, vir_key):
    res_dic[res_key]['designation'] = vir_dic[vir_key]['entry_name']
    res_dic[res_key]['keyDB'] = vir_key
    res_dic[res_key]['comment'] = vir_dic[vir_key]['comment']
    res_dic[res_key]['alternative_designation'] = vir_dic[vir_key]['alternative_names']
    res_dic[res_key]['searched_dna_changes'] = vir_dic[vir_key]['dna_snp']
    res_dic[res_key]['searched_prot_changes'] = vir_dic[vir_key]['prot_snp']
    res_dic[res_key]['ATB_groups'] = vir_dic[vir_key]['antibiotic_family']
    res_dic[res_key]['mechanism_groups'] = vir_dic[vir_key]['mechanism_names']
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
    nw_res_dict = {}
    n = 0
    for key in res_dic.keys():
        for key1 in res_dic[key]['mutations'].keys():
            nw_res_dict[n] = {}
            for key2 in res_dic[key]['mutations'][key1].keys():
                nw_res_dict[n][key2] = res_dic[key]['mutations'][key1][key2]
            for key3 in res_dic[key].keys():
                if key3 != 'mutations':
                    if key3 == 'dna_sequence':
                        nw_res_dict[n][key3] = str(res_dic[key][key3].seq)
                    else:
                        nw_res_dict[n][key3] = res_dic[key][key3]
            n += 1

    arrays_items = ['keyDB', 'designation', 'pc_identity', 'pc_coverage', 'mean_depth', 'variant_seq_type',
                    'known_change', 'unknown_change',
                    'nucleotide_change', 'detected_nucleotide', 'depth_by_detected_nucleotide',
                    'ATB_groups', 'mechanism_groups',
                    'alternative_designation', 'searched_prot_changes', 'searched_dna_changes', 'comment',
                    'dna_sequence', 'prot_sequence']
    res_list = []
    sort_list = list(nw_res_dict.keys())
    sort_list.sort()
    for i in sort_list:
        line_list = []
        for item in arrays_items:
            if item in ['pc_identity', 'pc_coverage']:
                line_list.append(float(nw_res_dict[i][item]))
            else:
                if item == 'designation':
                    txt = ''
                    if 'Warning_HET' in nw_res_dict[i].keys():
                        txt = nw_res_dict[i]['Warning_HET'] + ','
                    if 'Warning_COV' in nw_res_dict[i].keys():
                        txt = txt + nw_res_dict[i]['Warning_COV'] + ','
                    txt = txt + nw_res_dict[i][item]
                    line_list.append(txt)
                else:
                    line_list.append(nw_res_dict[i][item])
        res_list.append(line_list)

    df = pd.DataFrame.from_records(res_list, columns=arrays_items)
    df.sort_values(by=['pc_coverage', 'pc_identity', 'keyDB', 'designation', 'known_change', 'unknown_change'],
                   ascending=[False, False, True, True, False, False], inplace=True)

    writer = pd.ExcelWriter(os.path.join(out_dir, f'results_{dt_basename}.xlsx'))
    df.to_excel(writer, 'sheet1', index=False)
    writer.save()
    df.to_csv(os.path.join(out_dir, f'results_{dt_basename}.tsv'), sep='\t', index=False)


def write_summary_result(res_dic, out_dir, dt_basename, sample_id):
    csv_dic = {}
    xls_dic = {}
    dna_records = []
    prot_records = []
    keys = list(res_dic.keys())
    keys.sort()
    for key in keys:
        sequence_id = f'{sample_id}__ctg_X__{sample_id}_locusX__{res_dic[key]["designation"]}__{res_dic[key]["keyDB"]}'
        description = f'func:{res_dic[key]["ATB_groups"]},mechanism:{res_dic[key]["mechanism_groups"]},' \
                      f'id:{res_dic[key]["pc_identity"]},cov:{res_dic[key]["pc_coverage"]},' \
                      f'dep:{res_dic[key]["mean_depth"]}'

        col_name = f'{res_dic[key]["ATB_groups"]}::{res_dic[key]["designation"]}::{res_dic[key]["keyDB"]}'

        val_name = ''
        if 'Warning_HET' in res_dic[key].keys():
            val_name = res_dic[key]['Warning_HET'] + ','
        if 'Warning_COV' in res_dic[key].keys():
            val_name = res_dic[key]['Warning_COV'] + ','
        val_name = val_name + f'id:{res_dic[key]["pc_identity"]},cov:{res_dic[key]["pc_coverage"]},' \
                              f'dep:{res_dic[key]["mean_depth"]}'

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
                snp_txt = snp_txt + f',snp:{"|".join(snp[key3])}[{key3}]'
            description = description + snp_txt
            val_name = val_name + snp_txt

        if snp_txt == '' and res_dic[key]['var_only'] == '1':
            val_name = val_name + ',snp:None'

        sub_txt = ''
        if sub != {}:
            for key3 in sub.keys():
                sub_txt = sub_txt + f',sub:{"|".join(sub[key3])}[{key3}]'
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

            dna_sequence = SeqRecord(res_dic[key]['dna_sequence'].seq, id=sequence_id, name=sequence_id,
                                     description=description)
            dna_records.append(dna_sequence)
            if res_dic[key]['prot_sequence'] != '.' and res_dic[key]['prot_sequence'] != '.':
                prot_sequence = SeqRecord(Seq(res_dic[key]['prot_sequence']), id=sequence_id, name=sequence_id,
                                          description=description)
                prot_records.append(prot_sequence)

    SeqIO.write(dna_records, open(os.path.join(out_dir, f'results_{dt_basename}.fna'), 'w'), 'fasta')
    SeqIO.write(prot_records, open(os.path.join(out_dir, f'results_{dt_basename}.faa'), 'w'), 'fasta')

    df = pd.DataFrame(csv_dic, index=[sample_id, ])
    df.sort_index(axis=1, inplace=True)
    df.to_csv(os.path.join(out_dir, f'summary_results_{dt_basename}.csv'), sep='\t', index=False)

    for key in csv_dic.keys():
        keys = key.split('::')
        try:
            xls_dic[keys[0]]['::'.join(keys[1:])] = csv_dic[key]
        except KeyError:
            xls_dic[keys[0]] = {'::'.join(keys[1:]): csv_dic[key]}

    sheet_names = list(xls_dic.keys())
    sheet_names.sort()
    writer = pd.ExcelWriter(os.path.join(out_dir, f'summary_results_{dt_basename}.xlsx'))
    for n, sheet_name in enumerate(sheet_names):
        df = pd.DataFrame(xls_dic[sheet_name], index=[sample_id, ])
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
    subgroup = args.subGroup

    # execution main
    main(sample_id, sample_file, setting_file, dt_base_type, wk_dir, db_path, subgroup)


def main(sample_id, sample_file, setting_file, dt_base_type, wk_dir, db_path, subgroup):
    log_message = ""
    if sample_file == '':
        sample_file = os.path.join(wk_dir, 'sample.csv')
    if wk_dir == '':
        wk_dir = os.path.dirname(sample_file)
    out_dir = os.path.join(wk_dir, sample_id)
    set_dic = read_setting_file(setting_file)
    sample_dic, sample_list = read_sample_file(sample_file)
    species = sample_dic[sample_id]
    dt_basename = dt_base_file = ""
    if species.lower() in set_dic:
        set_species = set_dic[species.lower()]
        dt_basename_pre_split = str(*set_species[dt_base_type]).split("_")
        dt_basename = f'{dt_basename_pre_split[0]}_ariba_{dt_basename_pre_split[1]}_{subgroup}'
        dt_base_file = os.path.join(db_path, "dbVIR", "subsets", "{0}_{1}.tsv".format(*set_species[dt_base_type], subgroup))
    else:
        # if no set its "all"
        for file in os.listdir(out_dir):
            if "virDB_ariba" in file:
                dt_basename = file
                dt_basename_split = dt_basename.split("_")
        dt_base_file = os.path.join(db_path, "dbVIR", "subsets", "virDB_{0}_all.tsv".format(dt_basename_split[2]))

    tsv_file = os.path.join(out_dir, dt_basename, 'report.tsv')
    gen_file = os.path.join(out_dir, dt_basename, 'assembled_genes.fa.gz')
    seq_file = os.path.join(out_dir, dt_basename, 'assembled_seqs.fa.gz')
    resu_file = os.path.join(out_dir, dt_basename, 'filtered_out_results.txt')

    res_dic = load_vir_res(tsv_file, gen_file, seq_file)

    res_dic, log_message_tmp = filter_results(res_dic, species, resu_file)
    log_message = log_message + log_message_tmp

    res_dic, log_message_tmp = check_allele(res_dic, dt_base_file)
    log_message = log_message + log_message_tmp

    write_csv_result(res_dic, out_dir, dt_basename)
    write_summary_result(res_dic, out_dir, dt_basename, sample_id)

    return log_message


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='parse__detection - Version ' + version())
    parser.add_argument('-s', '--sampleID', dest="sampleID", default='CNR1979', help='Sample ID')
    parser.add_argument('-sf', '--sampleFile', dest="sampleFile", default='sample.csv', help='Sample file')
    parser.add_argument('-d', '--dtbase', dest="dtbase", default='vir', help="Database name")
    parser.add_argument('-wd', '--wkDir', dest="wkDir", default='', help="Working directory")
    parser.add_argument('-st', '--settingFile', dest="settingFile", default='/usr/local/readmapper-v0.1/setting.txt',
                        help="Setting file")
    parser.add_argument('-sg', '--subGroup', dest="subGroup", default='default value in setting.txt',
                        help="Sub group of gene")
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
