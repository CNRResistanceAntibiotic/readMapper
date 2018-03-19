#!/usr/bin/python
import os
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from update_armDB import load_arm_db
from prepare_mapping import read_sample_file
from prepare_mapping import read_setting
import pandas as pd

def load_armRes(tsvFile, genFile, seqFile):

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
                print 'Not found:', dtDic['ctg']
                print dtDic
            
            ref_len = float(dtDic['ref_len'])
            ref_base_assembled = int(dtDic['ref_base_assembled'])
            cov = (ref_base_assembled / ref_len) * 100
            dtDic['pc_coverage'] = '%.2f' % cov

            ref_nt = dtDic['ref_nt']
            ctg_nt = dtDic['ctg_nt']
            nucleotide_change = '%s->%s' % (ref_nt, ctg_nt)
            if nucleotide_change == '.->.':
                nucleotide_change = '.'
            smtls_nts = dtDic['smtls_nts']
            smtls_nts_depth = dtDic['smtls_nts_depth']
            smtls_total_depth = dtDic['smtls_total_depth']

            subdepth = parse_nt_info(ref_nt, ctg_nt, smtls_nts, smtls_nts_depth, smtls_total_depth)

            if dtDic['ref_ctg_effect'] != '.':
                if dtDic['known_var_change'] == dtDic['ref_ctg_change']:
                    known_change = dtDic['ref_ctg_change']
                    unknown_change = '.'
                else:
                    known_change = '.'
                    unknown_change = dtDic['ref_ctg_change']
            else:
                known_change, unknown_change = '.', '.'


            if ref_name in resDic.keys():
                indexes = resDic[ref_name]['mutations'].keys()
                indexes.sort()
                index = indexes[-1] + 1
                if dtDic['var_type'] == 'HET':
                    resDic[ref_name]['Warning_HET'] =  'warning:heterozygote'   
                resDic[ref_name]['mutations'][index] = {'unknown_change':unknown_change, 'ref_ctg_change':dtDic['ref_ctg_change'],
                                   'known_change':known_change, 'searched_change':dtDic['known_var_change'],
                                   'has_known_var': dtDic['has_known_var'],
                                   'variant_type':dtDic['var_type'],'variant_seq_type':dtDic['var_seq_type'],
                                   'ref_ctg_effect':dtDic['ref_ctg_effect'], 'summary_substitution_depth': subdepth,
                                   #'ref_start':dtDic['ref_start'],
                                   #'ref_end':dtDic['ref_end'],
                                   #'ctg_nt':dtDic['ctg_nt']
                                   'nucleotide_change':nucleotide_change,
                                   #'ctg_start':dtDic['ctg_start'],
                                   #'ctg_end':dtDic['ctg_end'],
                                   'detected_nucleotide':dtDic['smtls_nts'],
                                   'depth_by_detected_nucleotide':dtDic['smtls_nts_depth'],
                                   'total_depth_by_nucleotide':dtDic['smtls_total_depth']
                                   }
            else:
                resDic[ref_name] = {'pc_identity':dtDic['pc_ident'],'pc_coverage':dtDic['pc_coverage'],
                                'gene':dtDic['gene'],
                                'var_only':dtDic['var_only'],
                                'flag':dtDic['flag'],
                                'ctg_name':dtDic['ctg'],
                                #'ctg_len':dtDic['ctg_len'],
                                'mean_depth':dtDic['ctg_cov'],
                                'known_var':dtDic['known_var'],
                                'dna_sequence':dna_rec,
                                'prot_sequence':prot_rec,
                                'mutations':{n:{'unknown_change':unknown_change, 'ref_ctg_change':dtDic['ref_ctg_change'],
                                   'known_change':known_change, 'the_searched_change':dtDic['known_var_change'],
                                   'has_known_var': dtDic['has_known_var'],
                                   'variant_type':dtDic['var_type'],'variant_seq_type':dtDic['var_seq_type'],
                                   'ref_ctg_effect':dtDic['ref_ctg_effect'],'summary_substitution_depth': subdepth,
                                   #'ref_start':dtDic['ref_start'],
                                   #'ref_end':dtDic['ref_end'],
                                   #'ref_nt':dtDic['ref_nt'],
                                   #'ctg_start':dtDic['ctg_start'],
                                   #'ctg_end':dtDic['ctg_end'],
                                   #'ctg_nt':dtDic['ctg_nt'],
                                   'nucleotide_change': nucleotide_change,
                                   'detected_nucleotide':dtDic['smtls_nts'],
                                   'depth_by_detected_nucleotide':dtDic['smtls_nts_depth'],
                                   'total_depth_by_nucleotide': dtDic['smtls_total_depth']
                                   }}}
                if dtDic['var_type'] == 'HET':
                    resDic[ref_name]['Warning_HET'] = 'warning:heterozygote'  
                    #print ref_name, resDic[ref_name]['Warning_HET']
                #print ref_name, float(dtDic['pc_coverage'])
                if float(dtDic['pc_coverage']) < 95.0:
                    #print 'ok'
                    resDic[ref_name]['Warning_COV'] =  'warning:truncation'
                    #print ref_name, resDic[ref_name]['Warning_COV']
    f.close()
    return resDic


def parse_nt_info(ref_nt, ctg_nt, smtls_nts, smtls_nts_depth, smtls_total_depth):

    txt = ''
    for i in range(0,len(smtls_total_depth.split(';'))):
        nt = smtls_nts.split(';')[i]
        depth = smtls_nts_depth.split(';')[i]
        total_depth = smtls_total_depth.split(';')[i]

        if ',' not in nt:
            # perc = int(depth) / float(total_depth)
            txt = txt + '|%s[%s/%s]' % (nt, depth, total_depth)
        else:
            for j in range(0,len(nt.split(','))):
                n = nt.split(',')[j]
                d = depth.split(',')[j]
                # perc = int(d) / float(total_depth)
                if j == 0:
                    txt = txt + '|%s[%s/%s]' % (n, d, total_depth)
                else:
                    txt = txt + ',%s[%s/%s]' % (n, d, total_depth)
    txt = txt[1:]

    txt = '%s->%s %s' % (ref_nt, ctg_nt, txt)
    if txt == '.->. .[./.]':
        txt = '.'

    return txt


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


def filter_results(resDic, species, resuFile,  passcov=80, passid=80):

    

    f = open(resuFile,'a')
    del_keys = []
    del_n_muts  = []
    del_no_muts = []
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


        # filter
        if resDic[key]['var_only'] == '1':
            #filter SNP detection in wrong taxonomy
            del_keys = taxon_snp(del_keys, key, species)
            f.write('Filter wrong taxonomy\n:')
            f.write('%s\n%s\n\n' % (key, resDic[key]))
            #filter SNP search with no SNP
            del_no_muts, del_keys = no_change(del_no_muts, del_keys, resDic, key)

        # filter nts SNP if not searched
        del_n_muts = unknown_synonymous(del_n_muts, resDic, key)


    # Deletion SNP in records
    #print 'Filtered %i unknown nucleotide SNP' % len(del_n_muts)
    for key, n in del_n_muts:
        f.write('Filter nts SNP if not searched\n:')
        f.write('%s\n%s\n%s\n\n' % (key, resDic[key], resDic[key]['mutations'][n]))
        del resDic[key]['mutations'][n]

    for key, n in del_no_muts:
        try:
            del resDic[key]['mutations'][n]
            if key != '':
                f.write('Filter SNP search with no SNP\n:')
                f.write('%s\n%s\n%s\n\n' % (key, resDic[key], resDic[key]['mutations'][n]))
        except:
            x = 0
    f.close()
    # Deletion of records
    #print 'Filtered results (%s): %i' % (species, len(del_keys))
    for key in list(set(del_keys)):
        del resDic[key]
    return resDic


def unknown_synonymous(del_n_muts, resDic, key):
    for n in resDic[key]['mutations'].keys():
        if resDic[key]['mutations'][n]['variant_seq_type'] == 'n' and  resDic[key]['mutations'][n]['known_change'] == '.':
            del_n_muts.append([key, n])
    return del_n_muts


def no_change(del_no_muts, del_keys, resDic, key):
    no_change_found = True
    for n in resDic[key]['mutations'].keys():
        if resDic[key]['mutations'][n]['known_change'] == '.' and resDic[key]['mutations'][n]['unknown_change'] == '.' :
            del_no_muts.append([key, n])
        else:
            no_change_found = False
    if no_change_found == True:
        del_keys.append(key)

    return del_no_muts, del_keys


def taxon_snp(del_keys, key, species, sep=','):
    pattern = re.compile('[0-9]+[.]{1}[0-9]+[_]{2}.+[_]{2}([A-Za-z0-9, _]+)[_]{1}')
    enterobacteriaceae  = ['citrobacter freundii','enterobacter cloacae','enterobacter aerogenes',
                           'enterobacter asburiae','enterobacter hormaechei','escherichia coli',
                           'hafnia alvei','klebsiella pneumoniae','salmonella enteritidis']
    enterobacter_cloacae_complex = ['enterobacter cloacae','enterobacter asburiae','enterobacter hormaechei',
                                    'enterobacter kobei','enterobacter mori','enterobacter ludwigii',
                                    'enterobacter xiangfangensis']
    match = pattern.match(key)
    if match:
        taxons = match.group(1).split(sep)
        for n,taxon in enumerate(taxons):
            taxons[n] = taxon.lower().replace('_',' ')
        if 'enterobacteriaceae' in taxons:
            taxons    = taxons + enterobacteriaceae
        if 'enterobacter cloacae complex' in taxons:
            taxons    = taxons + enterobacter_cloacae_complex

        if species.lower() not in taxons:
            del_keys.append(key)

    return del_keys


def check_allele(resDic, dtbasefile):
    armDic = load_arm_db(dtbasefile)

    for reskey in resDic.keys():
        iarmkey = reskey.split('__')[0]
        data = resDic[reskey]

        resSeq = str(data['dna_sequence'].seq)
        armSeq = str(armDic[iarmkey]['dna_sequence'])


        if resSeq == armSeq:
            resDic = add_res(resDic, reskey, armDic, iarmkey)
            #print '100% DNA :', reskey, iarmkey
        else:
            found = True
            if data['gene'] == '1':
                resSeq = translate_dna(data['dna_sequence'].seq)
                armSeq = translate_dna(Seq(armDic[iarmkey]['dna_sequence']))
                if str(resSeq) == str(armSeq) and resSeq != '':
                    resDic = update_res(resDic, reskey, armDic, iarmkey)
                    #print '100% PROT:', reskey, iarmkey
                elif resSeq != '':
                    for armkey in armDic.keys():
                        armSeq = translate_dna(Seq(armDic[armkey]['dna_sequence']))
                        if str(armSeq) == str(resSeq):
                            found = False
                            resDic = update_res(resDic, reskey, armDic, armkey)
                            #print '100% NEW PROT', reskey, armkey, armDic[armkey]['entry_name']
                            break
                if found == True:
                    resDic = add_res(resDic, reskey, armDic, iarmkey)

            else:
                resSeq = data['dna_sequence'].seq
                for armkey in armDic.keys():
                    armSeq = armDic[armkey]['dna_sequence']
                    if str(armSeq) == str(resSeq):
                        found = False
                        resDic = update_res(resDic, reskey, armDic, armkey)
                        break
                if found == True:
                    resDic = add_res(resDic, reskey, armDic, iarmkey)

    return resDic


def add_res(resDic, reskey, armDic, armkey):
    resDic[reskey]['designation'] = armDic[armkey]['entry_name']
    resDic[reskey]['keyDB'] = armkey
    resDic[reskey]['comment'] = armDic[armkey]['comment']
    resDic[reskey]['alternative_designation'] = armDic[armkey]['alternative_names']
    resDic[reskey]['searched_dna_changes'] = armDic[armkey]['dna_snp']
    resDic[reskey]['searched_prot_changes'] = armDic[armkey]['prot_snp']
    resDic[reskey]['ATB_groups'] = armDic[armkey]['function_grp_names']
    resDic[reskey]['mechanism_groups'] = armDic[armkey]['mechanism_names']
    resDic[reskey]['related_groups'] = armDic[armkey]['cluster90_grp_name']
    return resDic


def update_res(resDic, reskey, armDic, armkey):
    resDic[reskey]['designation'] = armDic[armkey]['entry_name']
    resDic[reskey]['keyDB'] = armkey
    resDic[reskey]['comment'] = armDic[armkey]['comment']
    resDic[reskey]['alternative_designation'] = armDic[armkey]['alternative_names']
    resDic[reskey]['searched_dna_changes'] = armDic[armkey]['dna_snp']
    resDic[reskey]['searched_prot_changes'] = armDic[armkey]['prot_snp']
    resDic[reskey]['ATB_groups'] = armDic[armkey]['function_grp_names']
    resDic[reskey]['mechanism_groups'] = armDic[armkey]['mechanism_names']
    resDic[reskey]['related_groups'] = armDic[armkey]['cluster90_grp_name']
    resDic[reskey]['pc_identity'] = '100.00'
    resDic[reskey]['pc_coverage'] = '100.00'
    resDic[reskey]['pc_coverage'] = '100.00'
    resDic[reskey]['variant_seq_type'] = '.'
    resDic[reskey]['known_change'] = '.'
    resDic[reskey]['unknown_change'] = '.'
    resDic[reskey]['summary_substitution_depth'] = '.'
    resDic[reskey]['nucleotide_change'] = '.'
    resDic[reskey]['detected_nucleotide'] = '.'
    resDic[reskey]['depth_by_detected_nucleotide'] = '.'
    resDic[reskey]['total_depth_by_nucleotide'] = '.'
    return resDic


def translate_dna(dnaObject, table='Bacterial', cds=True):
    try:
        prot = dnaObject.translate(table=table, cds=cds)
    except:
        try:
            prot = dnaObject.reverse_complement().translate(table=table, cds=cds)
        except:
            prot = ''
    return prot


def write_csv_result(resDic, outdir, dtbasename):

    nwresDic = {}
    n = 0
    for key in resDic.keys():
        for key1 in resDic[key]['mutations'].keys():
            nwresDic[n] = {}
            for key2 in resDic[key]['mutations'][key1].keys():
                nwresDic[n][key2] = resDic[key]['mutations'][key1][key2]
            for key3 in resDic[key].keys():
                if key3 != 'mutations':
                    if key3 == 'dna_sequence':
                        nwresDic[n][key3]=str(resDic[key][key3].seq)
                    else:
                        nwresDic[n][key3]=resDic[key][key3]
            n += 1

    arrays_items = ['keyDB', 'designation', 'pc_identity', 'pc_coverage', 'mean_depth','variant_seq_type',
                    'known_change', 'unknown_change',
                    #'summary_substitution_depth',
                    'nucleotide_change', 'detected_nucleotide', 'depth_by_detected_nucleotide',
                    #'total_depth_by_nucleotide',
                    'ATB_groups','mechanism_groups','related_groups',
                    'alternative_designation', 'searched_prot_changes', 'searched_dna_changes', 'comment',
                    'dna_sequence','prot_sequence']
    resList = []
    sortlist = nwresDic.keys()
    sortlist.sort()
    for i in sortlist:
        line_list = []
        for item in arrays_items:
            if item in ['pc_identity', 'pc_coverage']:
                line_list.append(float(nwresDic[i][item]))
            else:
                if item == 'designation':
                    txt = ''
                    if 'Warning_HET' in nwresDic[i].keys():
                        txt = nwresDic[i]['Warning_HET'] + ','
                    if 'Warning_COV' in nwresDic[i].keys():
                        txt = txt + nwresDic[i]['Warning_COV'] + ','
                    txt = txt + nwresDic[i][item]
                    line_list.append(txt)
                else:
                    line_list.append(nwresDic[i][item])
        resList.append(line_list)

    df = pd.DataFrame.from_records(resList, columns=arrays_items)
    df.sort_values(by=['pc_coverage','pc_identity','keyDB','designation','known_change','unknown_change'],
                   ascending = [False, False, True, True, False, False], inplace=True)

    writer = pd.ExcelWriter(os.path.join(outdir, 'results_%s.xlsx' % dtbasename))
    df.to_excel(writer, 'sheet1', index=False)
    writer.save()
    df.to_csv(os.path.join(outdir, 'results_%s.csv' % dtbasename), sep='\t', index=False)


def write_summary_result(resDic, outdir, dtbasename, sampleID):
    csvDic = {}
    xlsDic = {}
    dnaRecords = []
    protRecords = []
    keys = resDic.keys()
    keys.sort()
    for key in keys:
        sequence_ID = '%s__ctg_X__%s_locusX__%s__%s' % (sampleID, sampleID, resDic[key]['designation'], resDic[key]['keyDB'])
        description = 'func:%s,mechanism:%s,id:%s,cov:%s,dep:%s' % \
                      (resDic[key]['ATB_groups'], resDic[key]['mechanism_groups'], resDic[key]['pc_identity'],
                       resDic[key]['pc_coverage'], resDic[key]['mean_depth'])

        colname = '%s::%s::%s' % (resDic[key]['ATB_groups'], resDic[key]['designation'], resDic[key]['keyDB'])

        valname = ''
        if 'Warning_HET' in resDic[key].keys():
            valname = resDic[key]['Warning_HET'] + ','
        if 'Warning_COV' in resDic[key].keys():
            valname = resDic[key]['Warning_COV'] + ','
        valname = valname + 'id:%s,cov:%s,dep:%s' % (resDic[key]['pc_identity'], resDic[key]['pc_coverage'], resDic[key]['mean_depth'])

        snp, sub = {},{}
        if resDic[key]['var_only'] == '1':
            for key2 in resDic[key]['mutations'].keys():
                if resDic[key]['mutations'][key2]['variant_seq_type'] == 'p':
                    var_sequence_type = 'pt'
                elif resDic[key]['mutations'][key2]['variant_seq_type'] == 'n':
                    var_sequence_type = 'nt'
                else:
                    var_sequence_type = '?'

                if resDic[key]['mutations'][key2]['known_change'] != '.':
                    try:
                        snp[var_sequence_type].append(resDic[key]['mutations'][key2]['known_change'])
                    except KeyError:
                        snp[var_sequence_type] = [resDic[key]['mutations'][key2]['known_change']]

                if resDic[key]['mutations'][key2]['unknown_change'] != '.':
                    try:
                        sub[var_sequence_type].append(resDic[key]['mutations'][key2]['unknown_change'])
                    except KeyError:
                        sub[var_sequence_type] = [resDic[key]['mutations'][key2]['unknown_change']]
        snp_txt = ''
        if snp != {}:
            for key3 in snp.keys():
                snp_txt = snp_txt + ',snp:%s[%s]' % ('|'.join(snp[key3]),key3)
            description = description + snp_txt
            valname = valname + snp_txt

        if snp_txt == '' and resDic[key]['var_only'] == '1':
            valname = valname + ',snp:None'

        sub_txt = ''
        if sub != {}:
            for key3 in sub.keys():
                sub_txt = sub_txt + ',sub:%s[%s]' % ('|'.join(sub[key3]),key3)
            description = description + sub_txt
            valname = valname + sub_txt

        record = False
        if resDic[key]['var_only'] == '0':
            record = True
        elif snp != {}:
            record = True
        elif sub != {}:
            record = True

        if record == True:
            csvDic[colname] = valname

            dna_sequence = SeqRecord(resDic[key]['dna_sequence'].seq, id=sequence_ID, name=sequence_ID, description=description)
            dnaRecords.append(dna_sequence)
            if resDic[key]['prot_sequence'] != '.' and  resDic[key]['prot_sequence'] != '.':
                prot_sequence = SeqRecord(Seq(resDic[key]['prot_sequence']), id=sequence_ID, name=sequence_ID, description=description)
                protRecords.append(prot_sequence)

    SeqIO.write(dnaRecords, open(os.path.join(outdir, 'results_%s.fna' % dtbasename), 'w'),'fasta')
    SeqIO.write(protRecords, open(os.path.join(outdir, 'results_%s.faa' % dtbasename), 'w'),'fasta')

    df = pd.DataFrame(csvDic, index=[sampleID,])
    df.sort_index(axis=1, inplace=True)
    df.to_csv(os.path.join(outdir, 'summary_results_%s.csv' % dtbasename), sep='\t', index=False)


    for key in csvDic.keys():
        keys = key.split('::')
        try:
            xlsDic[keys[0]]['::'.join(keys[1:])] = csvDic[key]
        except KeyError:
            xlsDic[keys[0]]={'::'.join(keys[1:]):csvDic[key]}

    sheetnames = xlsDic.keys()
    sheetnames.sort()
    writer = pd.ExcelWriter(os.path.join(outdir, 'summary_results_%s.xlsx' % dtbasename))
    for n, sheetname in enumerate(sheetnames):
        df = pd.DataFrame(xlsDic[sheetname], index=[sampleID,])
        df.sort_index(axis=1, inplace=True)
        df.to_excel(writer, sheetname, index=False)
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
    dtbasefile = os.path.join(dtbasedir, dtbasename + '.csv')

    tsvFile = os.path.join(outdir, dtbasename, 'report.tsv')
    genFile = os.path.join(outdir, dtbasename, 'assembled_genes.fa.gz')
    seqFile = os.path.join(outdir, dtbasename, 'assembled_seqs.fa.gz')

    #Aurelien: j ai rajouter ce path la qui etait mal renseigne
    resuFile = os.path.join(outdir, dtbasename, 'filtered_out_results.txt')


    resDic = load_armRes(tsvFile, genFile, seqFile)
    #print 'Number of result: %i' % len(resDic)
    resDic = filter_results(resDic, species, resuFile)
    #print 'Number of filtered result: %i' % len(resDic)
    resDic = check_allele(resDic, dtbasefile)

    #print ''
    #for key in resDic.keys():
    #    print key, resDic[key]
    #    for key2 in resDic[key]['mutations'].keys():
    #        print key2, resDic[key]['mutations'][key2]
    #    print ''

    write_csv_result(resDic, outdir, dtbasename)
    write_summary_result(resDic, outdir, dtbasename, sampleID)


def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='parse__detection - Version ' + version())
    parser.add_argument('-s', '--sampleID', dest="sampleID", default='CNR1979', help='Sample ID')
    parser.add_argument('-sf','--sampleFile', dest="sampleFile", default='sample.csv', help='Sample file')
    parser.add_argument('-d', '--dtbase', dest="dtbase", default='arm', help="Database name")
    parser.add_argument('-wd','--wkDir',      dest="wkDir",      default='', help="Working directory")
    parser.add_argument('-st','--settingFile',dest="settingFile",default='/usr/local/readmapper-v0.1/setting.txt', help="Setting file")
    parser.add_argument('-v', '--verbose',    dest="verbose",    default="0",help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version',    action='version', version='parse_detection-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
