#!/usr/bin/python3
import argparse
import csv
import os
import subprocess
import multiprocessing
from readmapper_src.parser_readmapper import parse_mlst_detection, parse_arm_detection, parse_rep_detection, \
    parse_vir_detection
from readmapper_src import write_docx


def get_sample_id_list(sample_file):
    """
    Get the list of all id in the sample file
    :param sample_file: the path to the sample file
    :return: a list of id
    """

    sample_id_list = []
    with open(sample_file, 'r') as tsv_file:
        spam_reader = csv.reader(tsv_file, delimiter='\t')
        for row in spam_reader:
            if row:
                sample_id_list.append(row[0])

    return sample_id_list


def manage_ariba(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, db_type):
    """
    This function manage the run ariba for mlst database
    :param wk_dir: the working directory
    :param sample_id: the sample id treated
    :param sample_file: the sample file
    :param setting_file: the setting file
    :param db_path: the path of the database
    :param bash_file : the bash file to execute
    :param db_type : type of the database used
    :return: nothing
    """

    subgroup = os.path.basename(bash_file).split("_")[-1].split('.')[0]

    name = multiprocessing.current_process().name
    log_message = f"\n*********** {name} Starting ***************\n"

    log_message = log_message + f"{db_type.upper()} detection in process for {sample_id}... Take time (~1 to 5 min)\n"

    try:
        with open(bash_file, 'r') as f:
            read_lines = f.readlines()
            for line in read_lines:
                line.rstrip('\n')
                if "#!/bin/bash" != line:
                    process = subprocess.Popen(line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)\
                        .stdout.read()
                    log_message = log_message + process.decode("utf-8")
    except Exception as e:
        log_message = log_message + f"Exception in detection: {e}\n"

    finally:
        log_message = log_message + f"{db_type.upper()} detection done.\n"

    log_message = log_message + f"{db_type.upper()} parsing in process for {sample_id}.\n"

    try:
        if db_type == "mlst":
            log_message = log_message + parse_mlst_detection.main(sample_id, sample_file, setting_file, db_type,
                                                                  wk_dir, subgroup)
        if db_type == "arm":
            log_message = log_message + parse_arm_detection.main(sample_id, sample_file, setting_file, db_type,
                                                                 wk_dir, db_path, subgroup)
        if db_type == "rep":
            log_message = log_message + parse_rep_detection.main(sample_id, sample_file, setting_file, db_type,
                                                                 wk_dir, subgroup)

        if db_type == "vir":
            log_message = log_message + parse_vir_detection.main(sample_id, sample_file, setting_file, db_type,
                                                                 wk_dir, db_path, subgroup)

    except Exception as e:
        log_message = log_message + f"Exception in parse: {e}\n"

    finally:

        log_message = log_message + "{0} parsing done.\n".format(db_type.upper())
        os.remove(bash_file)

    log_message = log_message + "*********** {0} Exiting ***************\n".format(name)
    print(log_message, flush=True)


def pre_main(args):
    setting_file = args.setFile
    wk_dir = os.path.abspath(args.workDir)
    sample_file = args.sampleFile
    initial = args.initial
    db_path = args.databasePath

    # execute main
    main(setting_file, wk_dir, sample_file, initial, db_path)


def main(setting_file, wk_dir, sample_file, initial, db_path):
    sample_id_list = get_sample_id_list(sample_file)
    count = 0
    files = os.listdir(wk_dir + "/")
    call_ariba_files = []

    for f in files:
        if "__calling__" in f:
            call_ariba_files.append(f)

    # process ariba for each database selected
    for sample_id in sample_id_list:
        count += 1
        print(f"Sample {count}/{len(sample_id_list)}: {sample_id} \n", flush=True)
        jobs = []
        for call_ariba_file in call_ariba_files:
            if sample_id in call_ariba_file:
                bash_file = os.path.join(wk_dir, call_ariba_file)
                if "__calling__mlst" in call_ariba_file:
                    schema_mlst = call_ariba_file.split("__calling__mlst_")[1].split(".")[0]
                    manage_ariba(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "mlst")
                    """
                    p = multiprocessing.Process(
                        target=manage_ariba,
                        args=(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "mlst"),
                        name=f'mlst {sample_id} schema {schema_mlst}'
                    )
                    jobs.append(p)
                    p.start()
                    """
                elif "__calling__arm" in call_ariba_file:
                    manage_ariba(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "arm")
                    """
                    p = multiprocessing.Process(
                        target=manage_ariba,
                        args=(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "arm"),
                        name=f'arm {sample_id}'
                    )
                    jobs.append(p)
                    p.start()
                    """
                elif "__calling__rep" in call_ariba_file:
                    manage_ariba(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "rep")
                    """
                    p = multiprocessing.Process(
                        target=manage_ariba,
                        args=(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "rep"),
                        name=f'rep {sample_id}'
                    )
                    jobs.append(p)
                    p.start()
                    """
                elif "__calling__vir" in call_ariba_file:
                    manage_ariba(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "vir")
                    """
                    p = multiprocessing.Process(
                        target=manage_ariba,
                        args=(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "vir"),
                        name=f'vir {sample_id}'
                    )
                    jobs.append(p)
                    p.start()
                    """
        """
        # wait multiple jobs
        for job in jobs:
            job.join()
        """
        print("Write a report in docx file \n")
        write_docx.main(wk_dir, initial, sample_id)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='manager.py - Version ' + version())
    parser.add_argument('-sf', '--splFile', dest="sampleFile", default='/home/bacteriologie/ariba/test/sample.csv',
                        help='fasta file of database to use')
    parser.add_argument('-wd', '--wkDir', dest="workDir", default='/home/bacteriologie/ariba/test',
                        help='tsv file of database to use')
    parser.add_argument('-set', '--setFile', dest="setFile", default='/usr/local/readmapper-v0.1/setting.txt',
                        help="setting file [setting.txt]")
    parser.add_argument('-in', '--initial', dest="initial", default="RBO", help="Initial of user")
    parser.add_argument('-db', '--databasePath', dest="databasePath", default='',
                        help="Database directory path")
    parser.add_argument('-V', '--version', action='version', version='manager-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
