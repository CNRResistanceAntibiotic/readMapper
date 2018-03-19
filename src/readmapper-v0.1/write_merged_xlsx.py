#!/usr/bin/python
import os
import glob
import argparse
import collections
import pandas as pd
from write_docx import read_mlst_results_csv_file
from write_docx import read_summary_arm_results_csv_file
import datetime


def read_sample_file(samplefile, sep='\t'):
    if os.path.exists(samplefile) == True:
        sampleDic = collections.OrderedDict()
        f = open(samplefile)
        for line in f:
            line = line.strip()
            if line != '':
                print line
                line = line.split(sep)
                sampleDic[line[0]] = line[1]
        f.close()
        return sampleDic
    else:
        print '\nSample file %s not found\n' % samplefile
        exit(1)

def read_rep_results_csv_file(filename, sep='\t'):
    if os.path.exists(filename) == True:
        repDic = []
        f = open(filename,'r')
        for n,line in enumerate(f):
            line = line.strip().split(sep)
            if n == 0:
                header = line
            elif n == 1:
                data = line
                dataDic = dict(zip(header, data))
                repDic[dataDic['']] = dataDic
        f.close()
        return repDic
    else:
        print '\nNo armDB result file %s\n' % filename
        exit(1)

def write_merged_mlst_dic(merged_mlstDic, wkdir, initial):
    if merged_mlstDic.keys() != []:
        outfile = os.path.join(wkdir, 'merged_results_mlst_%s_%s.xlsx' % (datetime.date.today(), initial))
        writer = pd.ExcelWriter(outfile)
        for n, sheetname in enumerate(merged_mlstDic.keys()):
            df = pd.DataFrame(merged_mlstDic[sheetname])
            df.sort_index(axis=1, inplace=True)
            colnames_2 = list(set(df.columns.tolist()) - set(['sample_id', 'mlst_name', 'ST']))
            colnames_2.sort()
            colnames = ['sample_id', 'ST'] + colnames_2
            df = df[colnames]
            df.to_excel(writer, sheetname, index=False)
        writer.save()
        return outfile
    else:
        print 'No MLST data to merge\n'
        return ''


def write_merged_arm_dic(merged_armDic, merged_mlstDic, sampleDic, wkdir, initial):

    atb_dic = {}
    for sampleID in merged_armDic.keys():
        for atb in merged_armDic[sampleID].keys():
            datadic = merged_armDic[sampleID][atb]
            if atb not in atb_dic.keys():
                atb_dic[atb] = {}
                for sample_ID in merged_armDic.keys():
                    ST, found = '', '0'
                    for mlstname in merged_mlstDic.keys():
                        for mlst_dic in merged_mlstDic[mlstname]:
                            if sample_ID == mlst_dic['sample_id']:
                                ST = mlst_dic['ST']
                                found = '1'
                                break
                        if found == '1':
                            break
                    atb_dic[atb][sample_ID]={'sample_id': sample_ID, 'species': sampleDic[sample_ID], 'mlst':'ST-%s' % ST}
            for item in datadic.keys():
                key = item.split('_')[0]
                txt = ''
                for item2 in ['id','cov','dep','snp','sub']:
                    if item2 in datadic[item].keys():
                        txt = txt + '%s:%s,' % (item2, datadic[item][item2])
                if 'warning' in datadic[item].keys():
                    #print item, datadic[item]['warning']
                    txt = 'WARNING! ' +  txt
                txt = txt[:-1]
                value = item.replace('_',' ') + ' [%s]' % txt

                atb_dic[atb][sampleID].update({key: value})
                #print atb, sampleID, key, value

    #print atb_dic
    outfile = os.path.join(wkdir, 'merged_results_arm_%s_%s.xlsx' % (datetime.date.today(), initial))
    writer = pd.ExcelWriter(outfile)
    atb_list = atb_dic.keys()
    atb_list.sort()
    for n, sheetname in enumerate(atb_list):
        records = []
        sampleIDs = atb_dic[sheetname].keys()
        sampleIDs.sort()
        for sampleID in sampleIDs:
            records.append(atb_dic[sheetname][sampleID])
        df = pd.DataFrame.from_dict(records)
        df.sort_index(axis=1, inplace=True)
        colnames_2 = list(set(df.columns.tolist()) - set(['sample_id', 'species', 'mlst']))
        colnames_2.sort()
        colnames = ['sample_id', 'species', 'mlst'] + colnames_2
        df = df[colnames]
        df.to_excel(writer, sheetname, index=False)
    writer.save()
    return outfile


def write_merged_vir_xls(virfile_list, wkdir, initial):

    VF_list = []
    for virfile in virfile_list:
        xl = pd.ExcelFile(virfile)
        VF_list = VF_list + xl.sheet_names
    VF_list = list(set(VF_list))

    outfile = os.path.join(wkdir, 'merged_results_vir_%s_%s.xlsx' % (datetime.date.today(), initial))
    writer = pd.ExcelWriter(outfile)
    VF_list.sort()
    for VF in VF_list:
        df = pd.DataFrame()
        for virfile in virfile_list:
            xl = pd.ExcelFile(virfile)
            if VF in xl.sheet_names:
                df = df.append(xl.parse(VF))
        df.to_excel(writer, VF, index=False, index_label='sample_id')
    writer.save()

    return outfile

def write_merged_rep_xls(repfile_list, wkdir, initial):

    outfile = os.path.join(wkdir, 'merged_results_rep_%s_%s.xlsx' % (datetime.date.today(), initial))
    writer = pd.ExcelWriter(outfile)

    df = pd.DataFrame()
    for repfile in repfile_list:
        xl = pd.ExcelFile(repfile)
        df = df.append(xl.parse('replicons'))
    df.to_excel(writer, 'replicons', index=True, index_label='sample_id')
    writer.save()

    return outfile

def main(args):
    wkdir = os.path.abspath(args.workdir)
    initial = args.initial

    samplefile = args.samplefile
    if samplefile == '':
        samplefile  = os.path.join(wkdir, 'sample.csv')

    print '\nData from sample file %s' % samplefile

    sampleDic     = read_sample_file(samplefile, sep='\t')
    print 'Number of samples: %i\n' % len(sampleDic.keys())

    merged_mlstDic = collections.OrderedDict()
    merged_armDic = collections.OrderedDict()
    virfile_list = []
    repfile_list = []
    for datatype in ['mlst','arm','vir','rep']:
        for sampleID in sampleDic.keys():
            if datatype == 'mlst':
                try:
                    STfilename  = glob.glob(os.path.join(wkdir, sampleID, 'mlst_report.tsv'))[0]
                    stDic       = read_mlst_results_csv_file(STfilename, sep='\t')
                except IndexError:
                    stDic = ''
                    print 'For %s no MLST data found' % sampleID
                if stDic != '':
                    mlst_name = stDic['mlst_name']
                    try:
                        merged_mlstDic[mlst_name].append(stDic)
                    except KeyError:
                        merged_mlstDic[mlst_name] = [stDic]

            elif datatype == 'arm':
                try:
                    armfilename = glob.glob(os.path.join(wkdir, sampleID, 'summary_results_armDB_*.csv'))[0]
                    armDBname   = os.path.splitext(os.path.basename(armfilename).replace('summary_results_',''))[0]
                    armDic      = read_summary_arm_results_csv_file(armfilename, sep='\t')
                except IndexError:
                    armDic = ''
                    print 'For %s no ARM data found' % sampleID
                if armDic != '':
                    merged_armDic[sampleID] = armDic

            elif datatype == 'vir':
                try:
                    virfile_list = virfile_list + [glob.glob(os.path.join(wkdir, sampleID, 'summary_results_virDB_*.xlsx'))[0]]
                    #virDBname   = os.path.splitext(os.path.basename(virfilename).replace('summary_results_',''))[0]
                except IndexError:
                    print 'For %s no VIR data found' % sampleID

            elif datatype == 'rep':
                try:
                    repfile_list = repfile_list + [glob.glob(os.path.join(wkdir, sampleID, 'summary_results_repDB_*.xlsx'))[0]]
                    #repDBname   = os.path.splitext(os.path.basename(repfilename).replace('summary_results_',''))[0]
                except IndexError:
                    print 'For %s no REP data found' % sampleID

    print ''
    mlst_file = write_merged_mlst_dic(merged_mlstDic, wkdir, initial)
    if mlst_file != '':
        print 'MLST results were written in %s\n' % mlst_file

    arm_file = write_merged_arm_dic(merged_armDic, merged_mlstDic, sampleDic, wkdir, initial)
    print 'ARM results were written in %s\n' % arm_file

    if repfile_list != []:
        rep_file = write_merged_rep_xls(repfile_list, wkdir, initial)
        print 'REP results were written in %s\n' % rep_file

    if virfile_list != []:
        vir_file = write_merged_vir_xls(virfile_list, wkdir, initial)
        print 'VIR results were written in %s\n' % vir_file



def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='write_merged_xlsx.py - Version ' + version())
    parser.add_argument('-wd', '--workdir',    dest="workdir",    default='/home/bacteriologie/ariba/anses/mcr-3',    help='working directory')
    parser.add_argument('-sf', '--sampleFile', dest="samplefile", default='/home/bacteriologie/ariba/anses/mcr-3/sample.csv',    help='sample file')
    parser.add_argument('-in', '--initial',    dest="initial",    default='RBO', help="Initials of the user")
    parser.add_argument('-V',  '--version', action='version', version='write_merged_xlsx-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
	run()
