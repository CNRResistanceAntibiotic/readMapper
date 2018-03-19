#!/usr/bin/python
import os
import glob
import argparse
from update_armDB import load_arm_db
from update_armDB import write_ariba_db
import datetime


def read_setting(inpfile, sep='\t'):
    setDic = {}
    f = open(inpfile,'r')
    for line in f:
        line = line.strip()
        if line != '':
            key, data = line.split(sep)
            setDic[key] = {}
            dataList = data.split(',')
            for data in dataList:
                key2, data2 = data.split(':')
                #if '|' in data2:
                setDic[key][key2] = data2.split('|')
                #else:
                #    setDic[key][key2] = data2
    f.close()
    return setDic


def read_sample_file(inpfile, sep='\t'):
    sampleDic = {}
    sampleList = []
    f = open(inpfile, 'r')
    for line in f:
        if line.strip() != '':
            sampleName, species = line.strip().split(sep)
            sampleDic[sampleName] = species
            sampleList.append(sampleName)
    f.close()
    return sampleDic, sampleList


def mlst_calling(dbdir, sampleID, reads1, reads2, outdir, exeDir='/usr/local/bin'):
    dtbasename = os.path.basename(dbdir)
    outfile = os.path.join(outdir, '%s__calling__%s.sh' % (sampleID, dtbasename))
    outdir = os.path.join(outdir, sampleID, dtbasename)
    cmd = os.path.join(exeDir, 'ariba run %s %s %s %s ' % (os.path.join(dbdir,'ref_db'), reads1, reads2, outdir))
    with open(outfile,'w') as f:
        f.write('#!/bin/bash\n')
        f.write(cmd)

    cmd = 'chmod a+x %s' % outfile
    os.system(cmd)


def gene_calling(dbdir, sampleID, reads1, reads2, outdir, exeDir='/usr/local/bin'):
    dtbasename = os.path.basename(dbdir)
    outfile = os.path.join(outdir, '%s__calling__%s.sh' % (sampleID, dtbasename))
    outdir = os.path.join(outdir, sampleID, dtbasename)
    cmd = os.path.join(exeDir, 'ariba run %s %s %s %s\n' % (dbdir, reads1, reads2, outdir))
    #print cmd
    with open(outfile,'w') as f:
        f.write('#!/bin/bash\n')
        f.write('mkdir -p %s\n' % os.path.dirname(outdir))
        f.write(cmd)
    #f.write('cut -f1,4 %s | uniq' % outdir)
    cmd = 'chmod a+x %s' % outfile
    os.system(cmd)


def main(args):
    settingfile = args.setFile
    wkdir       = os.path.abspath(args.workDir)
    readsDir    = args.readsDir
    samplefile  = args.sampleFile
    force       = args.force
    initial     = args.initial
    #nucmer_min_id = args.nucmer_min_id


    if samplefile == '':
        samplefile  = os.path.join(wkdir, 'sample.csv')

    setDic = read_setting(settingfile)
       
    sampleDic, sampleList = read_sample_file(samplefile)

    for sampleID in sampleList:
        outDir = os.path.join(wkdir, sampleID)
        species = sampleDic[sampleID].lower()
        if os.path.exists(outDir) == False:
            os.makedirs(outDir)
        elif os.path.exists(outDir) == True and force == True:
            cmd = 'rm %s -Rf; mkdir -p %s' % (outDir, outDir)
            os.system(cmd)	
        f = open(os.path.join(outDir, 'sample.csv'), 'w')
        f.write('%s\t%s' % (sampleID, species))
        f.close()

        reads1, reads2 = glob.glob(os.path.join(readsDir, '%s_*.fastq*' % sampleID))
        #reads1, reads2 = trim_reads(reads1, reads2, method, outdir)
        for work in ['mlst','arm','rep','vir']:
            #mlstRes_file, armRes_file, repRes_file, virRes_file = '', '', '', ''
            if work == 'mlst':
                if 'mlst' in setDic[species]:
                    print 'Prepare ST detection for %s' % sampleID
                    mlstDir = setDic['dbDir']['mlst'][0]
                    mlstDB  = setDic[species][work]
                    mlstDB = os.path.join(mlstDir, '%s_%s' % (mlstDB[0], mlstDB[1]))
                    f = open('/usr/local/readmapper-v0.1/info/mlst_trace.log','a')
                    f.write('%s\t%s\t%s\t%s\n' % (datetime.date.today(), sampleID, mlstDB, initial))
                    f.close()
                    mlst_calling(mlstDB, sampleID, reads1, reads2, wkdir)

            elif work == 'arm':
                if 'arm' in setDic[species]:
                    print 'Prepare antibiotic resistance gene detection for %s' % sampleID
                    armDir  = setDic['dbDir']['arm'][0]
                    #armFile = os.path.join(armDir, setDic[species]['arm'][0])
                    db_name = setDic[species]['arm'][0]
                    #db_ver  = os.path.splitext(setDic[species]['arm'][0])[0].split('_')[-1]
                    db_subset = setDic[species]['arm'][1]
                    armPrefix = '%s_%s' % (os.path.join(armDir, db_name), db_subset)
                    if os.path.isdir(armPrefix) == False:
                        print 'Database directory %s not found' % armPrefix
                    else:
                    	f = open('/usr/local/readmapper-v0.1/info/arm_trace.log','a')
                    	f.write('%s\t%s\t%s\t%s\n' % (datetime.date.today(), sampleID, armPrefix, initial))
                    	f.close()
                        gene_calling(armPrefix, sampleID, reads1, reads2, wkdir)

            elif work == 'rep':
                if 'rep' in setDic[species]:
                    print 'Prepare replicon detection for %s' % sampleID
                    repDir  = setDic['dbDir']['rep'][0]
                    db_name = setDic[species]['rep'][0]
                    db_subset = setDic[species]['rep'][1]
                    repPrefix = '%s_%s' % (os.path.join(repDir, db_name), db_subset)
                    if os.path.isdir(repPrefix) == False:
                        print 'Database directory %s not found' % repPrefix
                    else:
                    	f = open('/usr/local/readmapper-v0.1/info/rep_trace.log','a')
                    	f.write('%s\t%s\t%s\t%s\n' % (datetime.date.today(), sampleID, repPrefix, initial))
                    	f.close()
                        gene_calling(repPrefix, sampleID, reads1, reads2, wkdir)

            elif work == 'vir':
                if 'vir' in setDic[species]:
                    print 'Prepare virulence detection for %s' % sampleID
                    virDir  = setDic['dbDir']['vir'][0]
                    db_name  = setDic[species]['vir'][0]
                    db_subset = setDic[species]['vir'][1]
                    virPrefix = '%s_%s' % (os.path.join(virDir, db_name), db_subset)
                    if os.path.isdir(virPrefix) == False:
                        print 'Database directory %s not found' % virPrefix
                    else:
                    	f = open('/usr/local/readmapper-v0.1/info/vir_trace.log','a')
                    	f.write('%s\t%s\t%s\t%s\n' % (datetime.date.today(), sampleID, virPrefix, initial))
                    	f.close()
                        gene_calling(virPrefix, sampleID, reads1, reads2, wkdir)


def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='prepare_mapping.py - Version ' + version())
    parser.add_argument('-sf', '--splFile', dest="sampleFile", default='/home/bacteriologie/ariba/test/sample.csv', help='fasta file of database to use')
    parser.add_argument('-rd', '--rdDir', dest="readsDir", default='/home/bacteriologie/ariba/test', help='fasta file of database to use')
    parser.add_argument('-wd', '--wkDir', dest="workDir",  default='/home/bacteriologie/ariba/test', help='tsv file of database to use')
    #parser.add_argument('-nmi', '--nucmer_min_id', dest="nucmer_min_id", default='90', help="nucmer_min_id [90]")
    parser.add_argument('-set', '--setFile', dest="setFile", default='/usr/local/readmapper-v0.1/setting.txt', help="setting file [setting.txt]")
    parser.add_argument('-F', '--force', dest="force", default=False, action='store_true', help="Overwrite output directory, if it already exists [False]")
    parser.add_argument('-in', '--initial', dest="initial", default="RBO", help="Initial of user")
    parser.add_argument('-V', '--version', action='version', version='prepare_mapping-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
