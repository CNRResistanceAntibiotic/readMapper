#!/usr/bin/python3
import csv
import os
import glob
import argparse
import datetime

from Bio import SeqIO


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
                    data2_list = data2.split('|')
                    work = data2_list[0]
                    if "&" in data2:
                        set_dic[key][key2] = {work: data2_list[1].split("&")}
                    else:
                        set_dic[key][key2] = {work: [data2_list[1]]}

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


def gene_calling(db_dir, sample_id, reads1, reads2, out_dir):
    dt_basename = os.path.basename(db_dir)
    outfile = os.path.join(out_dir, '{0}__calling__{1}.sh'.format(sample_id, dt_basename))
    out_dir = os.path.join(out_dir, sample_id, dt_basename)
    cmd = os.path.join('ariba run {0} {1} {2} {3}\n'.format(db_dir, reads1, reads2, out_dir))
    with open(outfile, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('mkdir -p {0}\n'.format(os.path.dirname(out_dir)))
        f.write(cmd)


def get_all_subsets(all_tsv_file):
    db_subset_list = []

    with open(all_tsv_file, 'r', encoding="utf8") as f:

        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:

            subset_list = row['db_names'].split('|')

            for subset in subset_list:

                if subset not in db_subset_list:
                    db_subset_list .append(subset)

    return db_subset_list


def get_prepareref_seq_ariba(db_subset_list, db_arm_path, db_name, subsets_name):
    db_name_split = db_name.split("_")
    db_name = db_name_split[0]+"_ariba_"+db_name_split[1]

    seq_filename = os.path.join(db_arm_path, "{0}_{1}.fa".format(db_name, subsets_name))

    if not os.path.exists(seq_filename):

        records_hash = {}

        for subset in db_subset_list:

            subset_seq_file = os.path.join(db_arm_path, "{0}_{1}.fa".format(db_name, subset))
            for seq_record in SeqIO.parse(subset_seq_file, "fasta"):
                records_hash[seq_record.id] = seq_record
        # write output fasta
        with open(seq_filename, 'w') as handle:
            SeqIO.write(records_hash.values(), handle, 'fasta')

    return seq_filename


def get_prepareref_tsv_ariba(db_subset_list, db_arm_path, db_name, subsets_name):
    db_name_split = db_name.split("_")
    db_name = db_name_split[0]+"_ariba_"+db_name_split[1]

    tsv_filename = os.path.join(db_arm_path, "{0}_{1}.tsv".format(db_name, subsets_name))

    if not os.path.exists(tsv_filename):
        tsv_hash = {}

        for subset in db_subset_list:

            subset_tsv_file = os.path.join(db_arm_path, "{0}_{1}.tsv".format(db_name, subset))
            with open(subset_tsv_file) as tsv:
                reader = csv.reader(tsv, delimiter='\t')
                for row in reader:
                    id_prin = row[0]+'-'+row[3]
                    tsv_hash[id_prin] = row

        # write output tsv
        with open(tsv_filename, 'w') as tsv:
            writer = csv.writer(tsv, delimiter='\t')
            for key, row in tsv_hash.items():
                writer.writerow(row)

    return tsv_filename


def pre_main(args):
    setting_file = args.setFile
    wk_dir = os.path.abspath(args.workDir)
    reads_dir = args.readsDir
    samplefile = args.sampleFile
    subset_list = args.subset.split(',')
    if args.force == "True":
        print("\nForce the preparation : \n", flush=True)
        force = args.force
    else:
        force = ""
    initial = args.initial
    # nucmer_min_id = args.nucmer_min_id

    # execute main
    main(setting_file, wk_dir, reads_dir, samplefile, force, initial, subset_list)


def main(setting_file, wk_dir, reads_dir, samplefile, force, initial, subset_list):
    db_dir = os.path.abspath(os.path.join(setting_file, os.pardir))

    if samplefile == '':
        samplefile = os.path.join(wk_dir, 'sample.csv')

    set_dic = read_setting_file(setting_file)

    sample_dic, sample_list = read_sample_file(samplefile)

    for sample_id in sample_list:

        out_dir = os.path.join(wk_dir, sample_id)
        species = sample_dic[sample_id].lower()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        elif os.path.exists(out_dir) and force:
            cmd = 'rm {0} -Rf; mkdir -p {1}'.format(out_dir, out_dir)
            os.system(cmd)

        with open(os.path.join(out_dir, 'sample.csv'), 'w') as f:
            f.write('{0}\t{1}'.format(sample_id, species))

        reads1, reads2 = glob.glob(os.path.join(reads_dir, '{0}_*.fastq*'.format(sample_id)))

        for work in ['mlst', 'arm', 'rep', 'vir']:

            if work == 'mlst':
                if 'mlst' in set_dic[species]:

                    print('Prepare ST detection for {0}'.format(sample_id), flush=True)

                    mlst_db_hash = set_dic[species][work]

                    for key, value in mlst_db_hash.items():

                        # if multiple instruction like 2 or more schemas MLST
                        for subset in value:
                            mlst_db_path = os.path.join(db_dir, "dbMLST", "{0}_{1}".format(key, subset))

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

                    for db_name, db_subset_list in set_dic[species][work].items():

                        subsets_name = ""
                        # get user selected subsets
                        if subset_list:
                            db_subset_list = subset_list

                        for subset in db_subset_list:
                            if subset == 'all':
                                subsets_name = 'all'
                                db_subset_list = get_all_subsets(os.path.join(db_dir, 'dbARM', db_name + "_all.tsv"))
                                break
                            if not subsets_name:
                                subsets_name = subset
                                continue
                            subsets_name = subsets_name + "-" + subset

                        print("\nList ARM-DB Subset: {0}\n".format(subsets_name))

                        db_name_split = db_name.split("_")
                        db_name_ariba = "{0}_ariba_{1}".format(db_name_split[0], db_name_split[1])

                        print("\nVersion ARM-DB: {0}\n".format(db_name_split[1]))

                        arm_db_ariba_path = os.path.join(db_dir, "dbARM", "ariba", "{0}_{1}"
                                                         .format(db_name_ariba, subsets_name))

                        arm_subset_tsv_global_path = os.path.join(db_dir, "dbARM", "{0}_{1}.tsv"
                                                                  .format(db_name, subsets_name))

                        if not os.path.exists(arm_subset_tsv_global_path):

                            with open(arm_subset_tsv_global_path, 'w') as out:

                                writer = ""
                                pivot = 1

                                for subset in db_subset_list:
                                    arm_subset_tsv_path = os.path.join(db_dir, "dbARM",
                                                                       "{0}_{1}.tsv".format(db_name, subset))

                                    with open(arm_subset_tsv_path, 'r') as input_tsv:
                                        reader = csv.DictReader(input_tsv, delimiter='\t')
                                        if pivot:
                                            writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter='\t')
                                            writer.writeheader()
                                            pivot = 0

                                        for row in reader:
                                            writer.writerow(row)

                        elif not os.path.isdir(arm_db_ariba_path):
                            print('Database directory {0} not found'.format(arm_db_ariba_path), flush=True)
                            print('The ariba database construction started...', flush=True)

                            db_arm_path = os.path.abspath(os.path.join(db_dir, "dbARM", 'ariba'))
                            args_f = get_prepareref_seq_ariba(db_subset_list, db_arm_path, db_name, subsets_name)
                            args_m = get_prepareref_tsv_ariba(db_subset_list, db_arm_path, db_name, subsets_name)

                            cmd = 'ariba prepareref -f {0} -m {1} {2}'.format(args_f, args_m, arm_db_ariba_path)
                            print(cmd)
                            os.system(cmd)

                            print('The ariba database construction finished.', flush=True)

                        with open(os.path.join(db_dir, 'info', 'arm_trace.log'), 'a') as f:
                            f.write('{0}\t{1}\t{2}\t{3}\n'
                                    .format(datetime.date.today(), sample_id, arm_db_ariba_path, initial))

                        # run the ARM calling
                        gene_calling(arm_db_ariba_path, sample_id, reads1, reads2, wk_dir)

            elif work == 'rep':
                if 'rep' in set_dic[species]:

                    print('Prepare replicon detection for {0} '.format(sample_id), flush=True)

                    for db_name, db_subset_list in set_dic[species][work].items():

                        # if multiple instruction like 2 or more schemas MLST
                        for subset in db_subset_list:

                            rep_db_path = db_dir + "/dbREP/{0}_{1}".format(db_name, subset)

                            db_name_split = db_name.split("_")

                            print("\nVersion REP: {0}\n".format(db_name_split[1]))

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

                    for db_name, db_subset_list in set_dic[species][work].items():

                        # if multiple instruction like 2 or more schemas MLST
                        for subset in db_subset_list:

                            vir_db_path = db_dir + "/dbVIR/{0}_{1}".format(db_name, subset)

                            db_name_split = db_name.split("_")

                            print("\nVersion VIR: {0}\n".format(db_name_split[1]))

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
    parser.add_argument('-sb', '--subset', dest="subset", default='GN',
                        help="Comma separated value of database subset for ARM Database like (GN,Eff)")
    parser.add_argument('-V', '--version', action='version', version='prepare_mapping-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
