#!/usr/bin/python3
import os
import glob
import argparse
import collections
import pandas as pd
import datetime

from readmapper.write_docx import read_mlst_results_csv_file, read_summary_arm_results_csv_file


def read_sample_file(samplefile, sep='\t'):
    if os.path.exists(samplefile):
        sampleDic = collections.OrderedDict()
        with open(samplefile) as f:
            for line in f:
                line = line.strip()
                if line != '':
                    print(line)
                    line = line.split(sep)
                    sampleDic[line[0]] = line[1]
        return sampleDic
    else:
        print('\nSample file {0} not found\n'.format(samplefile))
        exit(1)


def read_rep_results_csv_file(filename, sep='\t'):
    global header
    if os.path.exists(filename):
        repDic = []
        with open(filename, 'r') as f:
            for n, line in enumerate(f):
                line = line.strip().split(sep)
                if n == 0:
                    header = line
                elif n == 1:
                    data = line
                    dataDic = dict(zip(header, data))
                    repDic[dataDic['']] = dataDic
        return repDic
    else:
        print('\nNo armDB result file {0}\n'.format(filename))
        exit(1)


def write_merged_mlst_dic(merged_mlst_dic, wk_dir, initial):
    if merged_mlst_dic.keys():
        outfile = os.path.join(wk_dir, 'merged_results_mlst_{0}_{1}.xlsx'.format(datetime.date.today(), initial))
        writer = pd.ExcelWriter(outfile)
        for n, sheet_name in enumerate(merged_mlst_dic.keys()):
            df = pd.DataFrame(merged_mlst_dic[sheet_name])
            df.sort_index(axis=1, inplace=True)
            col_names_2 = list(set(df.columns.tolist()) - {'sample_id', 'mlst_name', 'ST'})
            col_names_2.sort()
            col_names = ['sample_id', 'ST'] + col_names_2
            df = df[col_names]
            df.to_excel(writer, sheet_name, index=False)
        writer.save()
        return outfile
    else:
        print('No MLST data to merge\n')
        return ''


def write_merged_arm_dic(merged_arm_dic, merged_mlst_dic, sample_dic, wk_dir, initial):
    atb_dic = {}
    for sampleID in merged_arm_dic.keys():
        for atb in merged_arm_dic[sampleID].keys():
            data_dic = merged_arm_dic[sampleID][atb]
            if atb not in atb_dic.keys():
                atb_dic[atb] = {}
                for sample_ID in merged_arm_dic.keys():
                    ST, found = '', '0'
                    for mlst_name in merged_mlst_dic.keys():
                        for mlst_dic in merged_mlst_dic[mlst_name]:
                            if sample_ID == mlst_dic['sample_id']:
                                ST = mlst_dic['ST']
                                found = '1'
                                break
                        if found == '1':
                            break
                    atb_dic[atb][sample_ID] = {'sample_id': sample_ID, 'species': sample_dic[sample_ID],
                                               'mlst': 'ST-{0}'.format(ST)}
            for item in data_dic.keys():
                key = item.split('_')[0]
                txt = ''
                for item2 in ['id', 'cov', 'dep', 'snp', 'sub']:
                    if item2 in data_dic[item].keys():
                        txt = txt + '{0}:{1},'.format(item2, data_dic[item][item2])
                if 'warning' in data_dic[item].keys():
                    # print item, data_dic[item]['warning']
                    txt = 'WARNING! ' + txt
                txt = txt[:-1]
                value = item.replace('_', ' ') + ' [{0}]'.format(txt)

                atb_dic[atb][sampleID].update({key: value})

    outfile = os.path.join(wk_dir, 'merged_results_arm_{0}_{1}.xlsx'.format(datetime.date.today(), initial))
    writer = pd.ExcelWriter(outfile)
    atb_list = list(atb_dic.keys())
    atb_list.sort()
    for n, sheet_name in enumerate(atb_list):
        records = []
        sampleIDs = list(atb_dic[sheet_name].keys())
        sampleIDs.sort()
        for sampleID in sampleIDs:
            records.append(atb_dic[sheet_name][sampleID])
        df = pd.DataFrame.from_dict(records)
        df.sort_index(axis=1, inplace=True)
        col_names_2 = list(set(df.columns.tolist()) - {'sample_id', 'species', 'mlst'})
        col_names_2.sort()
        col_names = ['sample_id', 'species', 'mlst'] + col_names_2
        df = df[col_names]
        df.to_excel(writer, sheet_name, index=False)
    writer.save()
    return outfile


def write_merged_vir_xls(vir_file_list, wk_dir, initial):
    VF_list = []
    for vir_file in vir_file_list:
        xl = pd.ExcelFile(vir_file)
        VF_list = VF_list + xl.sheet_names
    VF_list = list(set(VF_list))

    outfile = os.path.join(wk_dir, 'merged_results_vir_{0}_{1}.xlsx'.format(datetime.date.today(), initial))
    writer = pd.ExcelWriter(outfile)
    VF_list.sort()
    for VF in VF_list:
        df = pd.DataFrame()
        for vir_file in vir_file_list:
            xl = pd.ExcelFile(vir_file)
            if VF in xl.sheet_names:
                df = df.append(xl.parse(VF))
        df.to_excel(writer, VF, index=False, index_label='sample_id')
    writer.save()

    return outfile


def write_merged_rep_xls(rep_file_list, wk_dir, initial):
    outfile = os.path.join(wk_dir, 'merged_results_rep_{0}_{1}.xlsx'.format(datetime.date.today(), initial))
    writer = pd.ExcelWriter(outfile)

    df = pd.DataFrame()
    for rep_file in rep_file_list:
        xl = pd.ExcelFile(rep_file)
        df = df.append(xl.parse('replicons'))
    df.to_excel(writer, 'replicons', index=True, index_label='sample_id')
    writer.save()

    return outfile


def pre_main(args):

    wk_dir = os.path.abspath(args.workdir)
    initial = args.initial
    samplefile = args.samplefile

    # execution main
    main(wk_dir, initial, samplefile)


def main(wk_dir, initial, samplefile):

    if samplefile == '':
        samplefile = os.path.join(wk_dir, 'sample.csv')

    print('\nData from sample file {0}'.format(samplefile))

    sampleDic = read_sample_file(samplefile, sep='\t')
    print('Number of samples: {0}\n'.format(len(sampleDic.keys())))

    merged_mlstDic = collections.OrderedDict()
    merged_armDic = collections.OrderedDict()
    vir_file_list = []
    rep_file_list = []
    for data_type in ['mlst', 'arm', 'vir', 'rep']:
        for sampleID in sampleDic.keys():
            if data_type == 'mlst':
                try:
                    ST_filename = glob.glob(os.path.join(wk_dir, sampleID, 'mlst_report.tsv'))[0]
                    stDic = read_mlst_results_csv_file(ST_filename, sep='\t')
                except IndexError:
                    stDic = ''
                    print('For {0} no MLST data found'.format(sampleID))
                if stDic != '':
                    mlst_name = stDic['mlst_name']
                    try:
                        merged_mlstDic[mlst_name].append(stDic)
                    except KeyError:
                        merged_mlstDic[mlst_name] = [stDic]

            elif data_type == 'arm':
                try:
                    arm_filename = glob.glob(os.path.join(wk_dir, sampleID, 'summary_results_armDB_*.csv'))[0]
                    armDic = read_summary_arm_results_csv_file(arm_filename, sep='\t')
                except IndexError:
                    armDic = ''
                    print('For {0} no ARM data found'.format(sampleID))
                if armDic != '':
                    merged_armDic[sampleID] = armDic

            elif data_type == 'vir':
                try:
                    vir_file_list = vir_file_list + [
                        glob.glob(os.path.join(wk_dir, sampleID, 'summary_results_virDB_*.xlsx'))[0]]
                    # virDBname   = os.path.splitext(os.path.basename(virfilename).replace('summary_results_',''))[0]
                except IndexError:
                    print('For {0} no VIR data found'.format(sampleID))

            elif data_type == 'rep':
                try:
                    rep_file_list = rep_file_list + [
                        glob.glob(os.path.join(wk_dir, sampleID, 'summary_results_repDB_*.xlsx'))[0]]
                    # repDBname   = os.path.splitext(os.path.basename(repfilename).replace('summary_results_',''))[0]
                except IndexError:
                    print('For {0} no REP data found'.format(sampleID))

    print('')
    mlst_file = write_merged_mlst_dic(merged_mlstDic, wk_dir, initial)
    if mlst_file != '':
        print('MLST results were written in {0}\n'.format(mlst_file))

    arm_file = write_merged_arm_dic(merged_armDic, merged_mlstDic, sampleDic, wk_dir, initial)
    print('ARM results were written in {0}\n'.format(arm_file))

    if rep_file_list:
        rep_file = write_merged_rep_xls(rep_file_list, wk_dir, initial)
        print('REP results were written in {0}\n'.format(rep_file))

    if vir_file_list:
        vir_file = write_merged_vir_xls(vir_file_list, wk_dir, initial)
        print('VIR results were written in {0}\n'.format(vir_file))


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='write_merged_xlsx.py - Version ' + version())
    parser.add_argument('-wd', '--workdir', dest="workdir", default='/home/bacteriologie/ariba/anses/mcr-3',
                        help='working directory')
    parser.add_argument('-sf', '--sampleFile', dest="samplefile",
                        default='/home/bacteriologie/ariba/anses/mcr-3/sample.csv', help='sample file')
    parser.add_argument('-in', '--initial', dest="initial", default='RBO', help="Initials of the user")
    parser.add_argument('-V', '--version', action='version', version='write_merged_xlsx-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
