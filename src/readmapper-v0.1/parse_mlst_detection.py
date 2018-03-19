#!/usr/bin/python
import os
import argparse
from collections import OrderedDict as Odict
from prepare_mapping import read_sample_file
from prepare_mapping import read_setting
import pandas as pd

def load_mlstRes(mlstFile, resDic, sep='\t'):
    f = open(mlstFile,'r')
    for n,line in enumerate(f):
        line = line.strip().split(sep)
        if n == 0:
            header = line
        elif n == 1 :
            for i in range(0,len(header)):
                if i == 0:
                    resDic[header[i]] = line[i]
                else:
                    resDic['gene%i' % i] = '%s:%s' % (header[i], line[i])
    f.close()
    return resDic


def write_csv_result(resDic, outdir):

    header = '\t'.join(resDic.keys())
    data = '\t'.join(resDic.values())
    print ''
    print header
    print data

    mlst_name = resDic['mlst_name']
    df = pd.DataFrame(resDic, index=[resDic['sample_id'], ])
    writer = pd.ExcelWriter(os.path.join(outdir, 'mlst_report.xlsx'))
    df.to_excel(writer, mlst_name, index=False)
    writer.save()

    df.to_csv(os.path.join(outdir, 'mlst_report.tsv'), sep='\t', index=False)

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
    print '\nDetected species: %s' % species
    set_species = setDic[species.lower()]
    dtbasename = '%s_%s' % (set_species[dtbasetype][0], set_species[dtbasetype][1])

    tsvFile = os.path.join(outdir, dtbasename, 'mlst_report.tsv')

    resDic = Odict()
    resDic['sample_id'] = sampleID
    resDic['mlst_name'] = dtbasename.split('_')[1]
    resDic = load_mlstRes(tsvFile, resDic)


    write_csv_result(resDic, outdir)


def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='parse__detection - Version ' + version())
    parser.add_argument('-s', '--sampleID', dest="sampleID", default='CNR1977', help='Sample ID')
    parser.add_argument('-sf','--sampleFile', dest="sampleFile", default='sample.csv', help='Sample file')
    parser.add_argument('-d', '--dtbase', dest="dtbase", default='mlst', help="Database name")
    parser.add_argument('-wd','--wkDir',      dest="wkDir",      default='', help="Working directory")
    parser.add_argument('-st','--settingFile',dest="settingFile",default='/usr/local/readmapper-v0.1/setting.txt', help="Setting file")
    parser.add_argument('-v', '--verbose',    dest="verbose",    default="0",help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version',    action='version', version='parse_detection-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
