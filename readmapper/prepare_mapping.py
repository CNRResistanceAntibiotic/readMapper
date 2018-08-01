#!/usr/bin/python
import os
import glob
import argparse
import datetime


def read_setting_file(inp_file, sep='\t'):
    set_dic = {}

    with open(inp_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line != '':
                key, data = line.split(sep)
                set_dic[key] = {}
                dataList = data.split(',')
                for data in dataList:
                    key2, data2 = data.split(':')
                    set_dic[key][key2] = data2.split('|')
    return set_dic


def read_sample_file(inp_file, sep='\t'):
    sample_dic = {}
    sample_list = []

    with open(inp_file, 'r') as f:
        for line in f:
            if line.strip() != '':
                sampleName, species = line.strip().split(sep)
                sample_dic[sampleName] = species
                sample_list.append(sampleName)
    return sample_dic, sample_list


def mlst_calling(db_dir, sample_id, reads1, reads2, out_dir):
    dt_basename = os.path.basename(db_dir)
    outfile = os.path.join(out_dir, '{0}__calling__{1}.sh'.format(sample_id, dt_basename))
    out_dir = os.path.join(out_dir, sample_id, dt_basename)
    cmd = os.path.join('ariba run {0} {1} {2} {3} '.format(os.path.join(db_dir, 'ref_db'), reads1, reads2, out_dir))

    with open(outfile, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write(cmd)
    # cmd = 'chmod a+x {0}'.format(outfile)
    # os.system(cmd)


def gene_calling(db_dir, sample_id, reads1, reads2, out_dir):
    dt_basename = os.path.basename(db_dir)
    outfile = os.path.join(out_dir, '{0}__calling__{1}.sh'.format(sample_id, dt_basename))
    out_dir = os.path.join(out_dir, sample_id, dt_basename)
    cmd = os.path.join('ariba run {0} {1} {2} {3}\n'.format(db_dir, reads1, reads2, out_dir))
    with open(outfile, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('mkdir -p {0}\n'.format(os.path.dirname(out_dir)))
        f.write(cmd)
    # cmd = 'chmod a+x {0}'.format(outfile)

    # os.system(cmd)


def pre_main(args):

    setting_file = args.setFile
    wk_dir = os.path.abspath(args.workDir)
    reads_dir = args.readsDir
    samplefile = args.sampleFile
    if args.force == "True":
        print("\nForce the preparation : \n", flush=True)
        force = args.force
    else:
        force = ""
    initial = args.initial
    # nucmer_min_id = args.nucmer_min_id

    # execute main
    main(setting_file, wk_dir, reads_dir, samplefile, force, initial)


def main(setting_file, wk_dir, reads_dir, samplefile, force, initial):

    db_dir = os.path.abspath(os.path.join(setting_file, os.pardir))

    if samplefile == '':
        samplefile = os.path.join(wk_dir, 'sample.csv')

    set_dic = read_setting_file(setting_file)

    sample_dic, sample_list = read_sample_file(samplefile)

    for sample_id in sample_list:

        outDir = os.path.join(wk_dir, sample_id)
        species = sample_dic[sample_id].lower()

        if not os.path.exists(outDir):
            os.makedirs(outDir)
        elif os.path.exists(outDir) and force:
            cmd = 'rm {0} -Rf; mkdir -p {1}'.format(outDir, outDir)
            os.system(cmd)

        with open(os.path.join(outDir, 'sample.csv'), 'w') as f:
            f.write('{0}\t{1}'.format(sample_id, species))

        reads1, reads2 = glob.glob(os.path.join(reads_dir, '{0}_*.fastq*'.format(sample_id)))

        # reads1, reads2 = trim_reads(reads1, reads2, method, outdir)

        for work in ['mlst', 'arm', 'rep', 'vir']:

            # mlstRes_file, armRes_file, repRes_file, virRes_file = '', '', '', ''

            if work == 'mlst':
                if 'mlst' in set_dic[species]:

                    print('Prepare ST detection for {0}'.format(sample_id), flush=True)

                    mlst_db = set_dic[species][work]
                    mlst_db_path = db_dir + "/dbMLST/{0}_{1}".format(mlst_db[0], mlst_db[1])

                    if not os.path.isdir(mlst_db_path):
                        print('Database directory {0} not found'.format(mlst_db_path), flush=True)
                    else:
                        with open(db_dir + '/info/mlst_trace.log', 'a') as f:
                            f.write('{0}\t{1}\t{2}\t{3}\n'
                                    .format(datetime.date.today(), sample_id, mlst_db_path, initial))
                        # run the MLST calling
                        mlst_calling(mlst_db_path, sample_id, reads1, reads2, wk_dir)

            elif work == 'arm':
                if 'arm' in set_dic[species]:

                    print('Prepare antibiotic resistance gene detection for {0}'.format(sample_id), flush=True)

                    db_name = set_dic[species][work][0]
                    db_subset = set_dic[species][work][1]

                    arm_db_path = db_dir + "/dbARM/{0}_{1}".format(db_name, db_subset)

                    if not os.path.isdir(arm_db_path):
                        print('Database directory {0} not found'.format(arm_db_path), flush=True)
                    else:
                        with open(db_dir + '/info/arm_trace.log', 'a') as f:
                            f.write('{0}\t{1}\t{2}\t{3}\n'
                                    .format(datetime.date.today(), sample_id, arm_db_path, initial))

                        # run the ARM calling
                        gene_calling(arm_db_path, sample_id, reads1, reads2, wk_dir)

            elif work == 'rep':
                if 'rep' in set_dic[species]:

                    print('Prepare replicon detection for {0} '.format(sample_id), flush=True)

                    db_name = set_dic[species][work][0]
                    db_subset = set_dic[species][work][1]

                    rep_db_path = db_dir + "/dbREP/{0}_{1}".format(db_name, db_subset)

                    if not os.path.isdir(rep_db_path):
                        print('Database directory {0} not found'.format(rep_db_path), flush=True)
                    else:
                        with open(db_dir + '/info/rep_trace.log', 'a') as f:
                            f.write('{0}\t{1}\t{2}\t{3}\n'
                                    .format(datetime.date.today(), sample_id, rep_db_path, initial))

                        # run the REP calling
                        gene_calling(rep_db_path, sample_id, reads1, reads2, wk_dir)

            elif work == 'vir':
                if 'vir' in set_dic[species]:

                    print('Prepare virulence detection for {0}'.format(sample_id), flush=True)

                    db_name = set_dic[species][work][0]
                    db_subset = set_dic[species][work][1]

                    vir_db_path = db_dir + "/dbVIR/{0}_{1}".format(db_name, db_subset)

                    if not os.path.isdir(vir_db_path):
                        print('Database directory {0} not found'.format(vir_db_path), flush=True)
                    else:
                        with open(db_dir + '/info/vir_trace.log', 'a') as f:
                            f.write('{0}\t{1}\t{2}\t{3}\n'
                                    .format(datetime.date.today(), sample_id, vir_db_path, initial))

                        # run the VIR calling
                        gene_calling(vir_db_path, sample_id, reads1, reads2, wk_dir)


def version():
    return "1.0.1"


def run():
    parser = argparse.ArgumentParser(description='prepare_mapping.py - Version ' + version())
    parser.add_argument('-sf', '--splFile', dest="sampleFile", default='/home/bacteriologie/ariba/test/sample.csv',
                        help='fasta file of database to use')
    parser.add_argument('-rd', '--rdDir', dest="readsDir", default='/home/bacteriologie/ariba/test',
                        help='fasta file of database to use')
    parser.add_argument('-wd', '--wkDir', dest="workDir", default='/home/bacteriologie/ariba/test',
                        help='tsv file of database to use')
    # parser.add_argument('-nmi', '--nucmer_min_id', dest="nucmer_min_id", default='90', help="nucmer_min_id [90]")
    parser.add_argument('-set', '--setFile', dest="setFile", default='/usr/local/readmapper-v0.1/setting.txt',
                        help="setting file [setting.txt]")
    parser.add_argument('-F', '--force', dest="force", default=False, action='store_true',
                        help="Overwrite output directory, if it already exists [False]")
    parser.add_argument('-in', '--initial', dest="initial", default="RBO", help="Initial of user")
    parser.add_argument('-V', '--version', action='version', version='prepare_mapping-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
