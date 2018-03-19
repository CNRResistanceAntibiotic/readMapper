#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/06/2017

import argparse
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os, re
import time, datetime

#def os_cmd(cmd):
    #subprocess.Popen(['/bin/bash', '-c', cmd])
    #subprocess.Popen(cmd, shell=True, executable='/bin/bash')


def load_arm_db(inpfile):
    armDic = {}
    f = open(inpfile, 'r')
    for n, line in enumerate(f):
        line = line.strip()
        if n == 0:
            header = line.split('\t')
        elif line != '' and n > 0:
            line = line.split('\t')
            data = dict(zip(header,line))
            for key in header:
                if key not in data.keys():
                    data[key]=''
            key = data['key']
            armDic[key] = data
    f.close()
    print '\nLoading of %s done!' % inpfile
    print 'Number of records: %i' % len(armDic.keys())
    return armDic


def write_ariba_db(armDic, db_name, db_ver, outdir, db='all'):
    keys = armDic.keys()
    records = []
    txt = ''
    keys.sort()
    for key in keys:
        data = armDic[key]
        #if armVer in data['db_names'].split('|'):
        #print len(data.keys())
        if db in data['db_names'].split('|'):
            name = data['entry_name']
            seq  = Seq(data['dna_sequence'])
            ID   = '%s::%s' % (key, name.replace(' ','_').replace('\'','PRIME').replace('(','[').replace(')',']'))
            record = SeqRecord(id=ID, name=ID, description='', seq=seq)
            records.append(record)
            mol_type = data['mol_type']
            if mol_type == 'cds':
                mol_type = 1
            else:
                mol_type = 0
            comment = data['comment']
            if comment.strip() == '':
                comment = ID

            prot_snp = data['prot_snp']
            if data['prot_snp'] != '':
                prot_snp = sort_snp(prot_snp.strip())
            dna_snp = data['dna_snp']
            if data['dna_snp'] != '':
                dna_snp = sort_snp(data['dna_snp'].strip())

            if prot_snp != '':
                search_type = 1
                #print 'snp', prot_snp
                for snp in prot_snp.split(','):
                    #print snp
                    txt = txt + '%s\t%s\t%s\t%s\t%s\t%s\n' % (ID, mol_type, search_type, snp, '.', comment)

            elif prot_snp == '' and dna_snp != '':
                search_type = 1
                for snp in dna_snp.split(','):
                    txt = txt + '%s\t%s\t%s\t%s\t%s\t%s\n' % (ID, mol_type, search_type, snp, '.', comment)

            else:
                search_type = 0
                txt = txt + '%s\t%s\t%s\t%s\t%s\t%s\n' % (ID, mol_type, search_type, '.', '.', comment)

    outPrefix = '%s_%s_%s' % (os.path.join(outdir, db_name), db_ver, db)
    SeqIO.write(records, open(outPrefix + '.fa','w'), 'fasta')

    f = open(outPrefix + '.tsv','w')
    f.write(txt)
    f.close()

    #cmd = 'ariba prepareref -f %s.fa -m %s.tsv %s' % (outPrefix, outPrefix, outPrefix)
    #os.system(cmd)

    return outPrefix


def dna_translate(dna_sequence, fmt):
    if fmt != 'obj':
        dna_sequence = Seq(dna_sequence)
    try:
        prot_sequence = dna_sequence.translate(table='Bacterial', cds=True)
    except:
        try:
            prot_sequence = dna_sequence[:-1].translate(table='Bacterial', cds=True)
            dna_sequence = dna_sequence[:-1]
        except:
            try:
                prot_sequence = dna_sequence.reverse_complement().translate(table='Bacterial', cds=True)
                dna_sequence = dna_sequence.reverse_complement()
            except:
                try:
                    prot_sequence = dna_sequence[:-1].reverse_complement().translate(table='Bacterial', cds=True)
                    dna_sequence = dna_sequence[:-1].reverse_complement()
                except:
                    prot_sequence = ''
    return str(dna_sequence), str(prot_sequence)


def write_armDB(armDic, outdir, tag='nw', subset='all'):
    colnames = ['key','source','entry_name','alternative_names','nomenclature_trouble','mol_type','dna_accession',
                'blastn_evalue','dna_snp','dna_sequence','prot_accession','blastp_evalue','prot_snp','prot_sequence',
                'function_grp_names','mechanism_names','operon_grp_name','cluster90_grp_name','comment','taxonomy',
                'reference','curated_by','db_names','update']
    outfile = os.path.join(outdir, 'armDB_%s_%s.csv' % (tag,subset))
    f = open(outfile, 'w')
    f.write('\t'.join(colnames)+'\n')
    keys = armDic.keys()
    keys.sort()
    for nb, key in enumerate(keys):
        dtDic = armDic[key]
        line = []
        if subset in dtDic['db_names'].split('|'):
            if dtDic['mol_type'] == 'cds':
                dna_sequence, prot_sequence = dna_translate(dtDic['dna_sequence'], 'str')
                if prot_sequence == '':
                    print '\nRecord %i among %i' % (nb, len(keys))
                    print 'Error in cds translation during the write of file armDB_new.cvs: %s (%s)' % (
                    key, dtDic['entry_name'])
                    print 'cds sequence converted in dna sequence'
                    dtDic['mol_type'] = 'dna'
                    dtDic['update'] = 'cdsTodna'
            else:
                dna_sequence = dtDic['dna_sequence']
                prot_sequence = ''

            for colname in colnames:
                try:
                    if colname == 'dna_sequence':
                        line.append(dna_sequence)
                    elif colname == 'prot_sequence':
                        line.append(prot_sequence)
                    elif colname == 'dna_snp':
                        snps = ','.join(list(set(dtDic[colname].split(','))))
                        snps = sort_snp(snps)
                        line.append(snps)
                    elif colname == 'prot_snp':
                        snps = ','.join(list(set(dtDic[colname].split(','))))
                        snps = sort_snp(snps)
                        line.append(snps)
                    else:
                        line.append(dtDic[colname])
                except KeyError:
                    line.append('')

            f.write('\t'.join(line)+'\n')
    f.close()


def sort_snp(snps):
    snpDic = {}
    pattern = re.compile('[A-Za-z]([0-9]+)[A-Za-z]')
    snps = list(set(snps.split(',')))
    for snp in snps:
        if snp != '':
            match = pattern.match(snp)
            pos = int(match.group(1))
            try :
                snpDic[pos].append(snp)
                snpDic[pos].sort()
            except KeyError:
                snpDic[pos] = [snp]

    positions = snpDic.keys()
    positions.sort()
    snps = ''
    for pos in positions:
        for snp in snpDic[pos]:
            snps = snps + ',' + snp
    snps = snps[1:]

    return snps


def load_dbs(fasfile, tsvfile):
    dbname = os.path.splitext(os.path.basename(tsvfile))[0]
    dbDic = {}
    fas = open(fasfile, 'r')
    for rec in SeqIO.parse(fas, 'fasta'):
        dbDic[rec.id] = {'rec':rec, 'tsv':''}
    fas.close()
    print 'Number: %i' % len(dbDic.keys())

    tsv = open(tsvfile, 'r')
    for line in tsv:
        line = line.strip().split('\t')
        key = line[0]
        if key not in dbDic.keys():
            dbDic[key]['tsv'] = line
        else:
            if dbDic[key]['tsv'] == '':
                dbDic[key]['tsv'] = line
            else:
                if dbDic[key]['tsv'][3] == '.' or dbDic[key]['tsv'][3] == line[3]:
                    dbDic[key]['tsv'][3] = '%s' % line[3]
                else:
                    dbDic[key]['tsv'][3] = '%s,%s' % (dbDic[key]['tsv'][3], line[3])
                #print dbDic[key]['tsv']
    tsv.close()
    print '\nLoading of %s and %s done!' % (fasfile, tsvfile)
    print 'Database name: %s' % dbname
    print 'Number of records: %i' % len(dbDic.keys())
    return dbname, dbDic


def update_armDB(armDic, dbDic, dbname, outdir):
    today = str(datetime.date.today())
    nb_updated_record = 0
    nb_new_record = 0

    for nb, dbkey in enumerate(dbDic.keys()):
        nb += 1
        db_dnaSeq  = str(dbDic[dbkey]['rec'].seq)
        db_moltype    = dbDic[dbkey]['tsv'][1]

        if db_moltype == '1':
            try:
                db_protSeq = str(dbDic[dbkey]['rec'].seq.translate(table='Bacterial', cds=True))
            except:
                try :
                    db_protSeq = str(dbDic[dbkey]['rec'].seq[:-1].translate(table='Bacterial', cds=True))
                    dbDic[dbkey]['rec'].seq = dbDic[dbkey]['rec'].seq[:-1]
                except:
                    try:
                        db_protSeq = str(dbDic[dbkey]['rec'].seq.reverse_complement().translate(table='Bacterial', cds=True))
                        dbDic[dbkey]['rec'].seq = dbDic[dbkey]['rec'].seq.reverse_complement()
                    except:
                        try:
                            db_protSeq = str(dbDic[dbkey]['rec'].seq[:-1].reverse_complement().translate(table='Bacterial', cds=True))
                            dbDic[dbkey]['rec'].seq = dbDic[dbkey]['rec'].seq[:-1].reverse_complement()
                        except:
                            print 'Error in translation in lines 107-121'
                            print '>%s' % dbkey
                            print dbDic[dbkey]['rec'].seq
                            exit(1)
        else:
            db_protSeq = 'db_proSeq_none'
        #db_searchtype = dbDic[dbkey]['tsv'][2]
        db_snp        = dbDic[dbkey]['tsv'][3]
        #db_grp        = dbDic[dbkey]['tsv'][4]
        #db_comment    = dbDic[dbkey]['tsv'][5]

        found = False

        for armkey in armDic.keys():
            arm_dnaSeq = armDic[armkey]['dna_sequence']
            #NON CODING SEQUENCE
            if db_moltype == '0':
                arm_protSeq = 'arm_proSeq_none'
                if arm_dnaSeq == str(db_dnaSeq) or str(db_dnaSeq) in arm_dnaSeq :
                    found = True
                    print '\n%i Non coding sequence %s already in database as reference %s (%s)' % (nb, dbkey, armkey, armDic[armkey]['entry_name'])
                    #Alternative name update
                    db_altname = '%s:%s' % (dbname, dbkey)
                    if db_altname not in armDic[armkey]['alternative_names'].split(','):
                        nb_updated_record += 1
                        print 'Add alternative name %s from %s in %s (%s)' % (db_altname, dbkey, armkey, armDic[armkey]['entry_name'])
                        if armDic[armkey]['alternative_names'].strip() == '':
                            armDic[armkey]['alternative_names'] = db_altname
                        else:
                            armDic[armkey]['alternative_names'] = '%s,%s' % (armDic[armkey]['alternative_names'], db_altname)
                        try:
                            armDic[armkey]['update'] = armDic[armkey]['update'] + ',alternative_names:%s:%s' % (today,db_altname)
                        except KeyError:
                            armDic[armkey]['update'] = 'alternative_names:%s:%s' % (today, db_altname)
                    #SNP update
                    if db_snp != '.' and arm_dnaSeq == db_dnaSeq:
                        db_snp = db_snp.replace('U','T')
                        if db_snp not in armDic[armkey]['dna_snp']:
                            nb_updated_record += 1
                            print 'Add SNP %s from %s in %s' % (db_snp, dbkey, armkey)
                            if armDic[armkey]['dna_snp'].strip() == '':
                                armDic[armkey]['dna_snp'] = db_snp
                            else:
                                armDic[armkey]['dna_snp'] = '%s,%s' % (armDic[armkey]['dna_snp'], db_snp)
                                armDic[armkey]['dna_snp'] = ','.join(list(set(armDic[armkey]['dna_snp'].split(','))))
                            try:
                                armDic[armkey]['update'] = armDic[armkey]['update'] + ',dna_snp:%s:%s' % (today, db_snp)
                            except KeyError:
                                armDic[armkey]['update'] = 'dna_snp:%s:%s' % (today, db_snp)
                    break

            #CODING SEQUENCE
            else:
                try :
                    arm_protSeq = Seq(armDic[armkey]['dna_sequence']).translate(table='Bacterial', cds=True)
                    armDic[armkey]['prot_sequence'] = arm_protSeq
                except:
                    try:
                        arm_protSeq = Seq(armDic[armkey]['dna_sequence'][:-1]).translate(table='Bacterial', cds=True)
                        armDic[armkey]['dna_sequence'] = armDic[armkey]['dna_sequence'][:-1]
                        armDic[armkey]['prot_sequence'] = arm_protSeq
                    except:
                        try:
                            arm_protSeq = Seq(armDic[armkey]['dna_sequence']).reverse_complement().translate(table='Bacterial', cds=True)
                            armDic[armkey]['dna_sequence'] = str(Seq(armDic[armkey]['dna_sequence']).reverse_complement())
                            armDic[armkey]['prot_sequence'] = arm_protSeq
                        except:
                            try:
                                arm_protSeq = Seq(armDic[armkey]['dna_sequence'][:-1]).reverse_complement().translate(table='Bacterial', cds=True)
                                armDic[armkey]['dna_sequence'] = str(Seq(armDic[armkey]['dna_sequence'][:-1]).reverse_complement())
                                armDic[armkey]['prot_sequence'] = arm_protSeq
                            except:
                                armDic[armkey]['prot_sequence'] = ''
                                arm_protSeq = 'db_proSeq_none'

                if arm_protSeq == db_protSeq or db_protSeq in arm_protSeq :
                    found = True
                    print '\n%i Coding sequence %s already in database as reference %s (%s)' % (nb, dbkey, armkey, armDic[armkey]['entry_name'])
                    #Alternative name update
                    db_altname = '%s:%s' % (dbname, dbkey)
                    if db_altname not in armDic[armkey]['alternative_names'].split(','):
                        print 'Add alternative name %s from %s in %s (%s)' % (db_altname, dbkey, armkey, armDic[armkey]['entry_name'])
                        nb_updated_record += 1
                        if armDic[armkey]['alternative_names'].split(',') == []:
                            armDic[armkey]['alternative_names'] = db_altname
                        else:
                            armDic[armkey]['alternative_names'] = '%s,%s' % (armDic[armkey]['alternative_names'], db_altname)
                        try:
                            armDic[armkey]['update'] = armDic[armkey]['update'] + ',alternative_names:%s' % db_altname
                        except KeyError:
                            armDic[armkey]['update'] = 'alternative_names:%s' % db_altname
                    #SNP update
                    if db_snp != '.' and arm_protSeq == db_protSeq:
                        if db_snp not in armDic[armkey]['prot_snp']:
                            print 'Add SNP %s from %s in %s' % (db_snp, dbkey, armkey)
                            nb_updated_record += 1
                            if armDic[armkey]['prot_snp'].split(',') == []:
                                armDic[armkey]['prot_snp'] = db_snp
                            else:
                                armDic[armkey]['prot_snp'] = '%s,%s' % (armDic[armkey]['prot_snp'], db_snp)
                                armDic[armkey]['prot_snp'] = ','.join(list(set(armDic[armkey]['prot_snp'].split(','))))
                            try:
                                armDic[armkey]['update'] = armDic[armkey]['update'] + ',prot_snp:%s:%s' % (today, db_snp)
                            except KeyError:
                                armDic[armkey]['update'] = 'prot_snp:%s' % db_snp
                    break

        #time.sleep(1)
        if found == False :
            #CREATE NEW ITEM
            print '\n%i Putative new record: %s' % (nb, dbkey)
            nwkey = generate_DBkey()
            dtDic = parse_db_tsv(dbkey, dbDic, dbname, armDic)
            dtDic['key'] = nwkey
            answer = ''
            while answer != 'Y' and answer != 'N':
                answer = raw_input('Confirm the update of armDB [Y]:\n')
                if answer == '':
                    answer = 'Y'
            if answer == 'Y':
                dtDic['update'] = 'create:%s' % today
                armDic[nwkey] = dtDic
                nb_new_record += 1
                print 'New sequence %s (%s) created from %s' % (nwkey, armDic[nwkey], dbkey)
                print armDic[nwkey]
                write_armDB(armDic, outdir, 'up')
            else:
                print 'The armDB database has not been updated with %s' % dbkey


    print '\nFrom database %s' %  dbname
    print 'Number of updates: %i' %  nb_updated_record
    print 'Number of new records: %i' % nb_new_record

    return armDic


def generate_DBkey():
    time.sleep(2)
    y = time.localtime(time.time()).tm_year
    m = time.localtime(time.time()).tm_mon
    j = time.localtime(time.time()).tm_mday
    h = time.localtime(time.time()).tm_hour
    min = time.localtime(time.time()).tm_min
    s = time.localtime(time.time()).tm_sec
    key = '%i%i%i.%i%i%i' % (y,m,j,h,min,s)
    return key


def parse_db_tsv(dbkey, dbDic, dbname, armDic):
    dtDic = {}
    data = dbDic[dbkey]['tsv']
    if dbname == 'card':
        entry_name = data[0].split('.')[0]
        dtDic['entry_name'] = entry_name
        dtDic['dna_accession'] = data[0].split('.')[2]
        dtDic['dna_sequence'] = str(dbDic[dbkey]['rec'].seq)

    elif dbname == 'resfinder':
        entry_name = data[0].split('.')[0]
        dtDic['entry_name'] = entry_name
        dtDic['dna_accession'] = data[0].split('_')[-1]
        dtDic['dna_sequence'] = str(dbDic[dbkey]['rec'].seq)

    elif dbname == 'argannot':
        pattern = compile('Original name: [(A-Za-z)]+([A-Za-z0-9_-]+)')
        match = pattern.match(data[5])
        dtDic['entry_name'] =  match.group(1)
        try :
            dtDic['dna_accession'] = data[5].split(':')[2]
        except IndexError:
            dtDic['dna_accession'] = data[5].split(':')[2]
        dtDic['dna_sequence'] = str(dbDic[dbkey]['rec'].seq)

    if data[1] == '0':
        mol_type = 'dna'
        dtDic['mol_type'] = 'dna'
        dtDic['prot_sequence'] = 'db_proSeq_none'
    else:
        mol_type = 'cds'
        dtDic['mol_type'] = 'cds'
        try:
            dtDic['prot_sequence'] = str(dbDic[dbkey]['rec'].seq.translate(table='Bacterial', cds=True))
        except:
            try:
                dtDic['prot_sequence'] = str(dbDic[dbkey]['rec'].seq.reverse_complement().translate(table='Bacterial', cds=True))
                dtDic['dna_sequence'] = str(dbDic[dbkey]['rec'].seq.reverse_complement())
            except:
                print 'Error translation line 274'
                print '>%s_%s' % (dbname, dbkey)
                print  dbDic[dbkey]['rec'].seq
                print '>%s_%s' % (dbname, dbkey)
                print  dbDic[dbkey]['rec'].seq.translate(table='Bacterial', cds=False)
                exit(1)

    dtDic['alternative_names'] = ''
    dtDic['source'] = '%s_%s' % (dbname, entry_name)

    dtDic['blastn_evalue']  = '1E-50'
    dtDic['prot_accession'] = ''
    dtDic['dna_snp']  = ''
    dtDic['prot_snp'] = ''
    if data[3] != '.' and mol_type == 'dna':
        dtDic['dna_snp']  = data[3].replace('U','T')
        dtDic['prot_snp'] = ''
        dtDic['blastp_evalue'] = ''
    elif data[3] != '.' and mol_type == 'cds':
        dtDic['dna_snp']  = ''
        dtDic['prot_snp'] = data[3]
        dtDic['blastp_evalue'] = '1E-30'

    dtDic['function_grp_names']	= ''
    dtDic['mechanism_names']	= ''
    dtDic['cluster90_grp_name']	= ''
    dtDic['db_names']           = ''
    dtDic['operon_grp_name']    = ''
    dtDic['comment']            = dbDic[dbkey]['tsv'][5]
    dtDic['taxonomy']           = ''
    dtDic['reference']          = ''
    dtDic['curated_by']	        = ''

    dtDic = blast_db(armDic, dtDic, dbkey, dbDic, dbname)

    return dtDic


def blast_db(armDic, dtDic, dbkey, dbDic, dbname):

    blastdir = os.path.join(os.getcwd(), 'tmp')
    if os.path.exists(blastdir) == False:
        cmd = 'mkdir -p %s' % blastdir
        os.system(cmd)
    qfile = os.path.join(blastdir, 'query.fasta')
    tfile = os.path.join(blastdir, 'target.fasta')
    trecords = []
    if dtDic['mol_type'] == 'cds':
        qrec  = SeqRecord(dbDic[dbkey]['rec'].seq.translate(table='Bacterial',cds=True),
                      id=dbDic[dbkey]['rec'].id, name=dbDic[dbkey]['rec'].id, description='')

        for key in armDic.keys():
            if armDic[key]['dna_sequence'].strip() != '':
                seq = Seq(armDic[key]['dna_sequence'])
                try:
                    rec = SeqRecord(seq.translate(table='Bacterial', cds=True), id=key, name=key, description='')
                    trecords.append(rec)
                except:
                    try:
                        rec = SeqRecord(seq.reverse_complement().transcribe(table='Bacterial', cds=True), id=key, name=key, description='')
                        trecords.append(rec)
                    except:
                        continue
    else:
        qrec  = SeqRecord(dbDic[dbkey]['rec'].seq,
                      id=dbDic[dbkey]['rec'].id, name=dbDic[dbkey]['rec'].id, description='')
        for key in armDic.keys():
            seq = Seq(armDic[key]['dna_sequence'])
            rec = SeqRecord(seq, id=key, name=key, description='')
            if str(rec.seq) != '':
                trecords.append(rec)

    SeqIO.write([qrec], open(qfile,'w'), 'fasta')
    SeqIO.write(trecords, open(tfile,'w'), 'fasta')

    ofile = os.path.join(os.getcwd(), 'tmp', 'hits.csv')
    cmd = '/usr/local/bin/usearch61 -usearch_global %s -db %s -id 0.9 -blast6out %s -strand plus > %s' % \
          (qfile, tfile, ofile, os.path.join(blastdir, 'usearch61.log'))
    print '\n%s' % cmd
    os.system(cmd)

    id_pass = 90
    cv_pass = 80
    target_id = ''
    for line in open(ofile, 'r'):
        line   = line.strip().split('\t')
        query  = line[0]
        target = line[1]
        id_perc = float(line[2])
        alg_len = int(line[3]) - int(line[5])
        q_cov   = alg_len / (float(len(qrec)))*100
        if id_perc >= id_pass and q_cov >= cv_pass:
            id_pass = id_perc
            cv_pass = q_cov
            target_id = target
    if target_id != '':
        print '\nPutative ortholog found in database for %s in %s: %s (%s - id: %.2f, cv: %.2f)' \
              % (dbkey, dbname, armDic[target]['entry_name'], target, id_pass, cv_pass)
        dtDic['function_grp_names']= armDic[target]['function_grp_names']
        dtDic['mechanism_names']   = armDic[target]['mechanism_names']
        dtDic['cluster90_grp_name']= armDic[target]['cluster90_grp_name']
        dtDic['db_names']          = armDic[target]['db_names']
        dtDic['operon_grp_name']   = armDic[target]['operon_grp_name']
    else:
        print '\nNo putative ortholog found in database:'
        print dbkey
        print '\t'.join(dbDic[dbkey]['tsv'])
        print str(dbDic[dbkey]['rec'].seq)

        print '\nAminoglycoside Colistin Cycline Divers Efflux Elfamycin Ethambutol Ethionamide Fosfomycin'
        print 'Glycopeptide Isoniazid Lincosamide Lipopeptide Macrolide Nitro-imidazole Oxazolidinone '
        print 'Phenicol Permeability Pleuromutilin Pyrazinamide Quinolone Rifampin'
        print 'Streptogramin Streptogramin_A Streptogramin_B Sulfonamide Trimethoprime'
        enter_function_grp_names = ''
        while enter_function_grp_names == '':
            print 'Item concatenation possible via \"|\"'
            enter_function_grp_names = raw_input('\nEnter function_grp_names:\n')
            dtDic['function_grp_names'] = enter_function_grp_names.strip()

        print '\nantibiotic inactivation enzyme'
        print 'antibiotic resistant gene variant or mutant'
        print 'antibiotic target modifying enzyme'
        print 'antibiotic target protection protein'
        print 'antibiotic target replacement protein'
        print 'efflux pump conferring antibiotic resistance'
        print 'gene altering cell wall charge conferring antibiotic resistance'
        print 'gene conferring antibiotic resistance via molecular bypass'
        print 'gene conferring resistance via absence'
        print 'gene modulating antibiotic efflux'
        print 'gene modulating permeability to antibiotic'
        print 'gene involved in antibiotic sequestration'
        print 'gene involved in self resistance to antibiotic'
        print 'unknown'
        enter_mechanism_names = ''
        while enter_mechanism_names == '':
            print '\nItem concatenation possible via \"|\"\n'
            enter_mechanism_names = raw_input('Enter mechanism_names:\n')
            dtDic['mechanism_names'] = enter_mechanism_names.strip()

        enter_cluster90_grp_name = raw_input('Enter cluster90_grp_name [%s]:\n' % dbkey)
        if enter_cluster90_grp_name == '':
            dtDic['cluster90_grp_name'] = dbkey.strip()
        else:
            dtDic['cluster90_grp_name'] = enter_cluster90_grp_name.strip()
        enter_operon_grp_name = raw_input('Enter operon_grp_name:\n')
        dtDic['operon_grp_name'] = enter_operon_grp_name.strip()

        print '\nall GN GP AGly BAAR Bla Div Eff Eth Eti Fos Fq Fus '
        print   'Gly Lpp Plx Iso MLS Nit Pem Phe Qln Rif Sul Tet Tmt'
        enter_db_names = ''
        while enter_db_names == '':
            print '\nItem concatenation possible via \"|\"\n'
            enter_db_names = raw_input('Enter db names:\n')
            dtDic['db_names'] = enter_db_names

        comment = raw_input('Enter a comment about the record [%s]:\n' % dtDic['comment'])
        if comment == '':
            comment = dtDic['comment']
        dtDic['comment'] = comment

    return dtDic


def main(args):

    fasfile = args.fasDB
    tsvfile = args.tsvDB
    armFile = os.path.abspath(args.armFile)
    armVer  = os.path.splitext(os.path.basename(armFile))[0].split('_')[-1]
    subsetDB = args.subDB
    if subsetDB == '':
        subsetDB = 'all'
    outDir  = os.path.dirname(armFile)


    if os.path.exists(os.path.join(outDir, 'armDB_nw.csv')) == True:
        print '\nFound an putative update of %s (%s)' % (armFile, os.path.join(outDir, 'armDB_up.csv'))
        update = raw_input('Should I use this version of the database for the update [Y] ?\n')
        if update == 'Y' or update == '':
            armFile = os.path.join(outDir, 'armDB_nw.csv')

    armDic = load_arm_db(armFile)

    if fasfile != '' and tsvfile != '':
        print '\nDatabase file: %s %s' % (fasfile, tsvfile)
        fasfile = os.path.abspath(fasfile)
        tsvfile = os.path.abspath(tsvfile)
        dbname, dbDic = load_dbs(fasfile, tsvfile)
        armDic = update_armDB(armDic, dbDic, dbname, outDir)

    write_armDB(armDic, outDir)
    write_armDB(armDic, outDir, 'nw', subsetDB)
    write_ariba_db(armDic, 'armDB', 'nw', outDir, subsetDB)

def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='Update armDB from card, resfinder and argannot - Version ' + version())
    parser.add_argument('-fas', '--fasDB', dest="fasDB", default='', help='fasta file of database to use [card.fa]')
    parser.add_argument('-tsv', '--tsvDB', dest="tsvDB", default='', help='tsv file of database to use [card.tsv]')
    parser.add_argument('-arm', '--armFile', dest="armFile", default='/usr/local/readmapper-v0.1/dbARM/armDB_1_all.csv', help="ARM file to update")
    parser.add_argument('-sub', '--subDB',   dest="subDB", default='GN', help="Subset of armDB")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0", help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='update_armDB-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
