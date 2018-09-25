#!/usr/bin/python3
import argparse
import csv
import os
import subprocess
import multiprocessing
from readmapper.parser_readmapper import parse_mlst_detection, parse_arm_detection, parse_rep_detection, \
    parse_vir_detection
from readmapper import write_docx


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


def manage_ariba(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, db_type, call_ariba_file):
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

    subgroup = call_ariba_file.split("_")[-1].split('.')[0]

    name = multiprocessing.current_process().name
    print("\n*********** ", name, 'Starting ***************\n', flush=True)

    print("{0} detection in process for {1}... Take time (~1 to 5 min)".format(db_type.upper(), sample_id), flush=True)

    try:
        with open(bash_file, 'r') as f:
            read_lines = f.readlines()
            for line in read_lines:
                line.rstrip('\n')
                if "#!/bin/bash" != line:
                    subprocess.call(line, shell=True)
    except Exception as e:
        print(e, flush=True)

    finally:

        print("{0} detection done.".format(db_type.upper()), flush=True)

    print("{0} parsing in process for {1}.".format(db_type.upper(), sample_id), flush=True)

    try:

        parser = ""

        if db_type == "mlst":
            # parser = parse_mlst_detection.__file__
            parse_mlst_detection.main(sample_id, sample_file, setting_file, db_type, wk_dir, subgroup)

        if db_type == "arm":
            # parser = parse_arm_detection.__file__

            parse_arm_detection.main(sample_id, sample_file, setting_file, db_type, wk_dir, db_path, subgroup)

        if db_type == "rep":
            # parser = parse_rep_detection.__file__
            parse_rep_detection.main(sample_id, sample_file, setting_file, db_type, wk_dir, subgroup)

        if db_type == "vir":
            # parser = parse_vir_detection.__file__
            parse_vir_detection.main(sample_id, sample_file, setting_file, db_type, wk_dir, db_path, subgroup)

        # give execution permission
        #  os.chmod(parser, 0o775)
        """
        cmd = parser + " -s {0} -sf {1} -wd {2} -d {3} -st {4} -db {5}".format(sample_id,
                                                                               sample_file, wk_dir,
                                                                               db_type,
                                                                               setting_file,
                                                                               db_path)
        print("Command parsing used : ", cmd, flush=True)
        subprocess.call(cmd, shell=True)
        """
    except Exception as e:
        print(e, flush=True)

    finally:

        print("{0} parsing done.".format(db_type.upper()), flush=True)

        os.remove(bash_file)

    print("\n*********** ", name, 'Exiting ***************\n', flush=True)


def pre_main(args):
    # Put executable permission for all users in all file
    # subprocess.call(['chmod', '-R', 'a+x', os.path.realpath(os.path.dirname(sys.argv[0]))])

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
        print("Sample {0}/{1}: {2} \n".format(count, len(sample_id_list), sample_id), flush=True)

        jobs = []

        for call_ariba_file in call_ariba_files:

            if sample_id in call_ariba_file:

                if "__calling__mlst" in call_ariba_file:
                    bash_file = os.path.join(wk_dir, call_ariba_file)

                    p = multiprocessing.Process(
                        target=manage_ariba,
                        args=(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "mlst", call_ariba_file),
                        name='mlst {0}'.format(sample_id)
                    )
                    jobs.append(p)
                    p.start()

                if "__calling__arm" in call_ariba_file:
                    bash_file = os.path.join(wk_dir, call_ariba_file)

                    p = multiprocessing.Process(
                        target=manage_ariba,
                        args=(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "arm", call_ariba_file),
                        name='arm {0}'.format(sample_id)
                    )
                    jobs.append(p)
                    p.start()

                if "__calling__rep" in call_ariba_file:
                    bash_file = os.path.join(wk_dir, call_ariba_file)

                    p = multiprocessing.Process(
                        target=manage_ariba,
                        args=(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "rep", call_ariba_file),
                        name='rep {0}'.format(sample_id)
                    )
                    jobs.append(p)
                    p.start()

                if "__calling__vir" in call_ariba_file:
                    bash_file = os.path.join(wk_dir, call_ariba_file)

                    p = multiprocessing.Process(
                        target=manage_ariba,
                        args=(wk_dir, sample_id, sample_file, setting_file, db_path, bash_file, "vir", call_ariba_file),
                        name='vir {0}'.format(sample_id)
                    )
                    jobs.append(p)
                    p.start()

        # wait multiple jobs
        for job in jobs:
            job.join()

        print("Write a report in docx file \n", flush=True)

        write_docx.main(wk_dir, initial, sample_id)

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", flush=True)


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
