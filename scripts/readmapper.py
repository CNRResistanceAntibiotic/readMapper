#!/usr/bin/python3
# Aurélien Birer
# 05/2018


import argparse
import os

from readmapper_src import prepare_mapping, manager, write_merged_xlsx


def main(args):
    subset_list = []
    print("Version ReadMapper: ", version())
    sample_file = os.path.abspath(args.sampleFile)
    reads = os.path.abspath(args.reads)
    wk_dir = os.path.abspath(args.workDir)
    database = os.path.abspath(args.database)
    initial = args.initial
    force = args.force
    if args.subset:
        subset_list = args.subset.split(',')

    if not sample_file:
        print("Sample file is missing !\n", flush=True)
        print(usage, flush=True)
        exit(1)

    if not reads:
        print("Reads directory is missing !\n", flush=True)
        print(usage, flush=True)
        exit(1)

    if not wk_dir:
        print("Working directory is missing !\n", flush=True)
        print(usage, flush=True)
        exit(1)

    if not database:
        print("Database directory is missing !\n", flush=True)
        print(usage, flush=True)
        exit(1)

    if not initial:
        print("Initial is missing !\n", flush=True)
        print(usage, flush=True)
        exit(1)

    dir_path = os.path.dirname(os.path.realpath(__file__))

    set_file = os.path.join(database, "setting.txt")

    #

    # print folders/files path
    print("\nSample file: {0}".format(sample_file), flush=True)
    print("Reads directory: {0}".format(reads), flush=True)
    print("Work directory: {0}".format(wk_dir), flush=True)
    print("DataBase directory: {0}".format(database), flush=True)
    print("Setting File: {0}".format(set_file), flush=True)
    print("Initial user: {0}".format(initial), flush=True)
    print("Subset selected: {0}".format(subset_list), flush=True)
    print("Force: {0}".format(force), flush=True)
    print("Application run at : {0}\n".format(dir_path), flush=True)

    print("----------------", flush=True)
    print("PREPARE THE JOBS", flush=True)
    print("----------------", flush=True)

    print("\nRun the preparation : \n", flush=True)

    try:
        prepare_mapping.main(set_file, wk_dir, reads, sample_file, force, initial, subset_list)
    except Exception as e:
        print(e, flush=True)

    finally:

        print("Preparation mapping done.", flush=True)

    print("\n-----------------", flush=True)
    print("STARTING THE JOBS", flush=True)
    print("-----------------", flush=True)

    print("\nRun the manager : \n", flush=True)

    try:
        manager.main(set_file, wk_dir, sample_file, initial, database)

    except Exception as e:
        print(e, flush=True)

    finally:

        print("Jobs finished.", flush=True)

    print("\n---------------", flush=True)
    print("MERGING RESULTS", flush=True)
    print("---------------", flush=True)

    try:
        write_merged_xlsx.main(wk_dir, initial, sample_file)

    except Exception as e:
        print(e, flush=True)

    finally:

        print("Jobs finished.", flush=True)


def version():
    return "1.0"


def run():
    global usage

    usage = "readmapper [-sf sample file] [-rd reads directory] [-wd work directory] [-dd databases directory] [-i " \
            "initial of the user] [-sb subset to use] <-F Overwrite output directory (Default=False)> "

    parser = argparse.ArgumentParser(
        prog='readmapper',
        usage=usage,
        description='ReadMapper: pipeline CNR Resistance with Ariba tool - Version ' + version(),
    )

    parser.add_argument('-sf', '--sampleFile', dest="sampleFile", default='', help='Sample file')
    parser.add_argument('-rd', '--readsDir', dest="reads", default='', help="Reads directory")
    parser.add_argument('-wd', '--wkDir', dest="workDir", default='',
                        help="Working directory")
    parser.add_argument('-dd', '--databaseDir', dest="database", default='',
                        help="Setting file")
    parser.add_argument('-sb', '--subset', dest="subset", default='',
                        help="Comma separated value of database subset for ARM Database like (GN,Eff)")
    parser.add_argument('-i', '--initial', dest="initial", default='',
                        help="Initial of user")
    parser.add_argument('-f', '--force', dest="force", default='False',
                        help="Overwrite output directory")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0",
                        help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='parse_rep_detection-' + version(),
                        help="Prints version number")

    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
