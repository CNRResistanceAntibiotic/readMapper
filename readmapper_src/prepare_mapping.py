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
                data_list = data.split(',')
                for data in data_list:
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
    outfile = os.path.join(out_dir, f'{sample_id}__calling__{dt_basename}.sh')
    out_dir = os.path.join(out_dir, sample_id, dt_basename)
    cmd = os.path.join(f'ariba run {os.path.join(db_dir, "ref_db")} {reads1} {reads2} {out_dir} ')

    with open(outfile, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write(cmd)


def gene_calling(db_dir, sample_id, reads1, reads2, out_dir):
    dt_basename = os.path.basename(db_dir)
    outfile = os.path.join(out_dir, f'{sample_id}__calling__{dt_basename}.sh')
    out_dir = os.path.join(out_dir, sample_id, dt_basename)
    cmd = os.path.join(f'ariba run {db_dir} {reads1} {reads2} {out_dir}\n')
    with open(outfile, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write(f'mkdir -p {os.path.dirname(out_dir)}\n')
        f.write(cmd)


def get_all_subsets(all_tsv_file):
    db_subset_list = []
    with open(all_tsv_file, 'r', encoding="utf8") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            subset_list = row['db_names'].split('|')
            for subset in subset_list:
                if subset not in db_subset_list:
                    db_subset_list.append(subset)
    return db_subset_list


def get_prepareref_seq_ariba(db_subset_list, db_arm_path, db_name, subsets_name):
    db_name_split = db_name.split("_")
    db_name = db_name_split[0] + "_ariba_" + db_name_split[1]
    seq_filename = os.path.join(db_arm_path, f"{db_name}_{subsets_name}.fa")
    if not os.path.exists(seq_filename):
        records_hash = {}
        for subset in db_subset_list:
            subset_seq_file = os.path.join(db_arm_path, f"{db_name}_{subset}.fa")
            for seq_record in SeqIO.parse(subset_seq_file, "fasta"):
                records_hash[seq_record.id] = seq_record
        # write output fasta
        with open(seq_filename, 'w') as handle:
            SeqIO.write(records_hash.values(), handle, 'fasta')
    return seq_filename


def get_prepareref_tsv_ariba(db_subset_list, db_arm_path, db_name, subsets_name):
    db_name_split = db_name.split("_")
    db_name = db_name_split[0] + "_ariba_" + db_name_split[1]
    tsv_filename = os.path.join(db_arm_path, f"{db_name}_{subsets_name}.tsv")
    if not os.path.exists(tsv_filename):
        tsv_hash = {}
        for subset in db_subset_list:
            subset_tsv_file = os.path.join(db_arm_path, f"{db_name}_{subset}.tsv")
            with open(subset_tsv_file) as tsv:
                reader = csv.reader(tsv, delimiter='\t')
                for row in reader:
                    id_prin = row[0] + '-' + row[3]
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
    sample_file = args.sampleFile
    subset_arm_list = args.subsetARM.split(',')
    subset_vir_list = args.subsetVIR.split(',')
    if args.force == "True":
        print("\nForce the preparation : \n", flush=True)
        force = args.force
    else:
        force = ""
    initial = args.initial
    # nucmer_min_id = args.nucmer_min_id
    # execute main
    main(setting_file, wk_dir, reads_dir, sample_file, force, initial, subset_arm_list, subset_vir_list)


def main(setting_file, wk_dir, reads_dir, sample_file, force, initial, subset_arm_list, subset_vir_list):
    db_dir = os.path.abspath(os.path.join(setting_file, os.pardir))
    if sample_file == '':
        sample_file = os.path.join(wk_dir, 'sample.csv')
    set_dic = read_setting_file(setting_file)
    sample_dic, sample_list = read_sample_file(sample_file)
    for sample_id in sample_list:
        out_dir = os.path.join(wk_dir, sample_id)
        species = sample_dic[sample_id].lower()
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        elif os.path.exists(out_dir) and force:
            cmd = f'rm {out_dir} -Rf; mkdir -p {out_dir}'
            os.system(cmd)
        with open(os.path.join(out_dir, 'sample.csv'), 'w') as f:
            f.write(f'{sample_id}\t{species}')
        reads1, reads2 = glob.glob(os.path.join(reads_dir, f'{sample_id}_*.fastq*'))
        for work in ['mlst', 'arm', 'rep', 'vir']:
            if work == 'mlst':
                if species in set_dic:
                    if 'mlst' in set_dic[species]:
                        print(f'Prepare ST detection for {sample_id}', flush=True)
                        mlst_db_hash = set_dic[species][work]
                        for key, value in mlst_db_hash.items():
                            # if multiple instruction like 2 or more schemas MLST
                            for subset in value:
                                mlst_db_path = os.path.join(db_dir, "dbMLST", f"{key}_{subset}")
                                if not os.path.isdir(mlst_db_path):
                                    print(f'Database directory {mlst_db_path} not found', flush=True)
                                else:
                                    with open(db_dir + '/info/mlst_trace.log', 'a') as f:
                                        f.write(f'{datetime.date.today()}\t{sample_id}\t{mlst_db_path}\t{initial}\n')
                                    # run the MLST calling
                                    mlst_calling(mlst_db_path, sample_id, reads1, reads2, wk_dir)
                else:
                    pass
            elif work == 'arm':
                subsets_name = db_name = ""
                db_subset_list = []
                if species in set_dic:
                    if 'arm' in set_dic[species]:
                        print(f'Prepare antibiotic resistance gene detection for {sample_id}', flush=True)
                        for db_name_tmp, db_subset_list in set_dic[species][work].items():
                            # get user selected subsets
                            if subset_arm_list:
                                db_subset_list = subset_arm_list
                            for subset in db_subset_list:
                                if subset == 'all':
                                    subsets_name = 'all'
                                    db_name = db_name_tmp
                                    db_subset_list = get_all_subsets(
                                        os.path.join(db_dir, 'dbARM', 'subsets', db_name_tmp + "_all.tsv"))
                                    break
                                if not subsets_name:
                                    subsets_name = subset
                                    continue
                                subsets_name = subsets_name + "-" + subset
                else:
                    subsets_name = "all"
                    for species, has_1 in set_dic.items():
                        for work, has_2 in has_1.items():
                            for db_name_tmp, db_subset_list in has_2.items():
                                if "arm" in db_name_tmp:
                                    db_name = db_name_tmp
                    db_subset_list = get_all_subsets(
                        os.path.join(db_dir, 'dbARM', 'subsets', db_name + "_all.tsv"))
                print(f"\nList ARM-DB Subset: {subsets_name}\n")
                db_name_split = db_name.split("_")
                db_name_ariba = f"{db_name_split[0]}_ariba_{db_name_split[1]}"
                print(f"\nVersion ARM-DB: {db_name_split[1]}\n")
                arm_db_ariba_path = os.path.join(db_dir, "dbARM", "ariba", f"{db_name_ariba}_{subsets_name}")

                arm_subset_tsv_global_path = os.path.join(db_dir, "dbARM", 'subsets', f"{db_name}_{subsets_name}.tsv")
                if not os.path.exists(arm_subset_tsv_global_path):
                    print(f"Create subset file in Database: {arm_subset_tsv_global_path}")
                    with open(arm_subset_tsv_global_path, 'w') as out:
                        writer = ""
                        pivot = 1
                        for subset in db_subset_list:
                            arm_subset_tsv_path = os.path.join(db_dir, "dbARM", 'subsets', f"{db_name}_{subset}.tsv")
                            with open(arm_subset_tsv_path, 'r') as input_tsv:
                                reader = csv.DictReader(input_tsv, delimiter='\t')
                                if pivot:
                                    writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter='\t')
                                    writer.writeheader()
                                    pivot = 0
                                for row in reader:
                                    writer.writerow(row)
                elif not os.path.isdir(arm_db_ariba_path):
                    print(f'Database directory {arm_db_ariba_path} not found', flush=True)
                    print('The ariba database construction started...', flush=True)
                    db_arm_path = os.path.abspath(os.path.join(db_dir, "dbARM", 'ariba'))
                    args_f = get_prepareref_seq_ariba(db_subset_list, db_arm_path, db_name, subsets_name)
                    args_m = get_prepareref_tsv_ariba(db_subset_list, db_arm_path, db_name, subsets_name)

                    cmd = f'ariba prepareref -f {args_f} -m {args_m} {arm_db_ariba_path}'
                    print(cmd)
                    os.system(cmd)
                    print('The ariba database construction finished.', flush=True)
                with open(os.path.join(db_dir, 'info', 'arm_trace.log'), 'a') as f:
                    f.write(f'{datetime.date.today()}\t{sample_id}\t{arm_db_ariba_path}\t{initial}\n')
                # run the ARM calling
                gene_calling(arm_db_ariba_path, sample_id, reads1, reads2, wk_dir)
            elif work == 'rep':
                if 'rep' in set_dic[species]:
                    print(f'Prepare replicon detection for {sample_id} ', flush=True)
                    for db_name, db_subset_list in set_dic[species][work].items():
                        # if multiple instruction like 2 or more schemas MLST
                        for subset in db_subset_list:
                            rep_db_path = db_dir + f"/dbREP/{db_name}_{subset}"
                            db_name_split = db_name.split("_")
                            print(f"\nVersion REP: {db_name_split[1]}\n")
                            if not os.path.isdir(rep_db_path):
                                print(f'Database directory {rep_db_path} not found', flush=True)
                            else:
                                with open(db_dir + '/info/rep_trace.log', 'a') as f:
                                    f.write(f'{datetime.date.today()}\t{sample_id}\t{rep_db_path}\t{initial}\n')
                                # run the REP calling
                                gene_calling(rep_db_path, sample_id, reads1, reads2, wk_dir)
                else:
                    pass
            elif work == 'vir':
                db_name = ""
                subsets_name = ""
                if 'vir' in set_dic[species]:
                    print(f'Prepare virulence detection for {sample_id}', flush=True)
                    for db_name, db_subset_list in set_dic[species][work].items():
                        # if multiple instruction like 2 or more schemas MLST
                        subsets_name = ""
                        # get user selected subsets
                        if subset_vir_list:
                            db_subset_list = subset_vir_list
                        for subset in db_subset_list:
                            if subset == 'all':
                                subsets_name = 'all'
                                db_subset_list = get_all_subsets(
                                    os.path.join(db_dir, 'dbVIR', 'subsets', db_name + "_all.tsv"))
                                break
                            if not subsets_name:
                                subsets_name = subset
                                continue
                            subsets_name = subsets_name + "-" + subset

                else:
                    subsets_name = 'all'
                    for species, has_1 in set_dic.items():
                        for work, has_2 in has_1.items():
                            for db_name_tmp, db_subset_list in has_2.items():
                                db_name = db_name_tmp

                print(f"\nList VIR-DB Subset: {subsets_name}\n")
                db_name_split = db_name.split("_")
                db_name_ariba = f"{db_name_split[0]}_ariba_{db_name_split[1]}"
                vir_db_path = os.path.join(db_dir, "dbVIR", "ariba", f"{db_name_ariba}_{subsets_name}")
                db_name_split = db_name.split("_")
                print(f"\nVersion VIR: {db_name_split[1]}\n")

                if not os.path.isdir(vir_db_path):
                    print(f'Database directory {vir_db_path} not found', flush=True)
                else:
                    with open(db_dir + '/info/vir_trace.log', 'a') as f:
                        f.write(f'{datetime.date.today()}\t{sample_id}\t{vir_db_path}\t{initial}\n')
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
    parser.add_argument('-set', '--setFile', dest="setFile", default='/usr/local/readmapper-v0.1/setting.txt',
                        help="setting file [setting.txt]")
    parser.add_argument('-F', '--force', dest="force", default=False, action='store_true',
                        help="Overwrite output directory, if it already exists [False]")
    parser.add_argument('-in', '--initial', dest="initial", default="RBO", help="Initial of user")
    parser.add_argument('-sbARM', '--subsetARM', dest="subsetARM", default='GN',
                        help="Comma separated value of database subset for ARM Database like (GN,Eff)")
    parser.add_argument('-sbVIR', '--subsetVIR', dest="subsetVIR", default='GN',
                        help="Comma separated value of database subset for VIR Database like (GN,Eff)")
    parser.add_argument('-V', '--version', action='version', version='prepare_mapping-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
