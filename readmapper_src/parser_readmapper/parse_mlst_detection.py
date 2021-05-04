#!/usr/bin/python3
import os
import argparse
from collections import OrderedDict as Odict
from readmapper_src.prepare_mapping import read_sample_file, read_setting_file
import pandas as pd


def load_mlst_res(mlst_file, res_dic, sep='\t'):
    with open(mlst_file, 'r') as f:
        header = ""
        for n, line in enumerate(f):
            line = line.strip().split(sep)
            if n == 0:
                header = line
            elif n == 1:
                for i in range(0, len(header)):
                    if i == 0:
                        res_dic[header[i]] = line[i]
                    else:
                        res_dic[f'gene{i}'] = f'{header[i]}:{line[i]}'
    return res_dic


def write_csv_result(res_dic, out_dir, subgroup):

    mlst_name = res_dic['mlst_name']
    df = pd.DataFrame(res_dic, index=[res_dic['sample_id'], ])
    writer = pd.ExcelWriter(os.path.join(out_dir, 'mlst_report_{}.xlsx'.format(subgroup)))
    df.to_excel(writer, mlst_name, index=False)
    writer.save()

    df.to_csv(os.path.join(out_dir, f'mlst_report_{subgroup}.tsv'), sep='\t', index=False)


def pre_main(args):
    sample_id = args.sampleID
    sample_file = args.sampleFile
    setting_file = args.settingFile
    dt_base_type = args.dtbase
    wk_dir = args.wkDir
    subgroup = args.subGroup

    # execution main
    main(sample_id, sample_file, setting_file, dt_base_type, wk_dir, subgroup)


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
    log_message = log_message + f"\nDetected species: {species}\n"
    set_species = set_dic[species.lower()]
    dt_basename = "{0}_{1}".format(*set_species[dt_base_type], subgroup)

    tsv_file = os.path.join(out_dir, dt_basename, 'mlst_report.tsv')

    res_dic = Odict()
    res_dic['sample_id'] = sample_id
    res_dic['mlst_name'] = dt_basename.split('_')[1]

    res_dic = load_mlst_res(tsv_file, res_dic)

    write_csv_result(res_dic, out_dir, subgroup)

    return log_message


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='parse__detection - Version ' + version())
    parser.add_argument('-s', '--sampleID', dest="sampleID", default='CNR1977', help='Sample ID')
    parser.add_argument('-sf', '--sampleFile', dest="sampleFile", default='sample.csv', help='Sample file')
    parser.add_argument('-d', '--dtbase', dest="dtbase", default='mlst', help="Database name")
    parser.add_argument('-wd', '--wkDir', dest="wkDir", default='', help="Working directory")
    parser.add_argument('-st', '--settingFile', dest="settingFile", default='/usr/local/readmapper-v0.1/setting.txt',
                        help="Setting file")
    parser.add_argument('-db', '--databasePath', dest="databasePath", default='',
                        help="Database directory path")
    parser.add_argument('-sg', '--subGroup', dest="subGroup", default='default value in setting.txt',
                        help="Sub group of gene")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0",
                        help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='parse_detection-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
