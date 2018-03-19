import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_virDB(virDB_file, sep='\t'):
    vfDic = {}
    geneDic = {}
    seqDic  = {}
    f = open(virDB_file,'r')
    for n,line in enumerate(f):
        if n == 0:
            header = line.strip().split(sep)
        else:
            data = dict(zip(header, line.strip().split(sep)))

            gene_id  = data['gene_id']
            gene_name = data['gene_name']
            geneDic[gene_id] = data

            dna_seq  = Seq(data['VF_sequence'])
            record = SeqRecord(dna_seq, id=gene_id, name=gene_name, description='')
            seqDic[gene_id] = record

            vf_id = data['VF_Accession']
            try :
                vfDic[vf_id].append(gene_id)
            except KeyError:
                vfDic[vf_id] = [gene_id]
    f.close()

    return geneDic, seqDic, vfDic


def clustering_sequence(seqDic, outprefix, force, idpass=0.90, covpass=0.90):
    print '\nStart sequence clustering...'
    cluster_file = outprefix + '.clt'
    if os.path.exists(cluster_file) == False or force == True:
        records = seqDic.values()
        SeqIO.write(records, open('tmp.fasta','w'), 'fasta')

        #if os.path.exists(cluster_file) == False:
        cmd = '/usr/local/bin/usearch61 -cluster_fast tmp.fasta -id %s -query_cov %s -target_cov %s -uc %s' % \
              (idpass, covpass, covpass, cluster_file)
        os.system(cmd)
        os.system('rm tmp.fasta')
    else:
        print "Warning: cluster file %s found" % cluster_file
        print "Please delete it to perform sequence clustering again"

    print '\nCollect the clusters in %s...' % cluster_file
    clusters = {}
    realclusters = {}
    for line in open(cluster_file, 'r'):
        line = line.strip().split('\t')
        gene_id = line[8]
        cluster_name = line[1]
        if gene_id in clusters.keys() and clusters[gene_id] != cluster_name:
            print '########'
            print 'Warning: %s in at least 2 clusters' % gene_id
            print '%s in %s' % (gene_id, clusters[gene_id])
            print '%s in %s' % (gene_id, cluster_name)
        clusters[gene_id] = cluster_name
        try:
            realclusters[cluster_name].append(gene_id)
        except KeyError:
            realclusters[cluster_name] = [gene_id]

    del_clusters = []
    for cluster_name in realclusters.keys():
        if len(list(set(realclusters[cluster_name]))) <= 1:
            del_clusters.append(cluster_name)

    for gene_id in clusters.keys():
        if clusters[gene_id] in del_clusters:
            clusters[gene_id] = ''
    print 'Number of genes: %i' % len(realclusters)
    print 'Number of genes with no homolog: %i' % len(del_clusters)
    print 'Number of clusters: %i' % (len(realclusters) - len(del_clusters))

    return clusters


def clustering_vf(geneDic, vfDic, clusDic):
    print "\nStart cluster parsing..."
    vf_list = vfDic.keys()
    vf_list.sort()
    vfCluster = {}
    for vf1 in vf_list:
        gene_list = vfDic[vf1]
        gene_list.sort()
        for gene_id in gene_list:
            q_id = '%s::%s::%s::%s::%s' % (geneDic[gene_id]['gene_id'], geneDic[gene_id]['gene_name'],
                                           geneDic[gene_id]['genus'], geneDic[gene_id]['species'],
                                           geneDic[gene_id]['strain'])
            q_vf_id = '%s::%s' % (geneDic[gene_id]['VF_Accession'], geneDic[gene_id]['VF_Name'])
            q_sequence = Seq(geneDic[gene_id]['VF_sequence'])

            if q_vf_id not in vfCluster.keys():
                vfCluster[q_vf_id] = {q_id:[]}
            elif q_id not in vfCluster[q_vf_id].keys():
                vfCluster[q_vf_id][q_id] = []

            try:
                q_cluster = clusDic[gene_id]
            except KeyError:
                #print 'No cluster for gene sequence %s' % gene_id
                q_cluster = ''

            if q_cluster != '':
                vf_list2 = vfDic.keys()
                vf_list2.sort()
                for vf2 in vf_list2:
                    gene_list2 = vfDic[vf2]
                    for gene_id2 in gene_list2:
                        t_id = '%s::%s::%s::%s::%s' % (geneDic[gene_id2]['gene_id'], geneDic[gene_id2]['gene_name'],
                                                       geneDic[gene_id2]['genus'], geneDic[gene_id2]['species'],
                                                       geneDic[gene_id2]['strain'])
                        t_vf_id = '%s::%s' % (geneDic[gene_id2]['VF_Accession'], geneDic[gene_id2]['VF_Name'])
                        t_sequence = Seq(geneDic[gene_id2]['VF_sequence'])

                        try:
                            t_cluster = clusDic[gene_id2]
                        except KeyError:
                            #print 'No cluster for gene sequence %s' % gene_id2
                            t_cluster = ''

                        if t_cluster == q_cluster:
                            if q_id != t_id:
                                vfCluster[q_vf_id][q_id].append(t_id)
                                vfCluster[q_vf_id][q_id] = list(set(vfCluster[q_vf_id][q_id]))

    return vfCluster


def write_cluster(vfCluster, outprefix):
    txt = ''
    for vf in vfCluster.keys():
        for gene in vfCluster[vf].keys():
            txt = txt + '%s\t%s\t%s\n' %  (vf, gene, ','.join(vfCluster[vf][gene]))
    f = open(outprefix+'.clu', 'w')
    f.write(txt)
    f.close()


def main(args):
    virDB_file = args.csvDB
    outprefix = os.path.splitext(virDB_file)[0]
    force = args.force

    geneDic, seqDic, vfDic  = load_virDB(virDB_file, sep='\t')
    clusDic = clustering_sequence(seqDic, outprefix, force, 0.80, 0.80)
    vfCluster = clustering_vf(geneDic, vfDic, clusDic)
    write_cluster(vfCluster, outprefix)

def version():
	return "1.0"


def run():
    parser = argparse.ArgumentParser(description='cluster_virDB - Version ' + version())
    parser.add_argument('-db', '--csvDB', dest="csvDB", default='/usr/local/readmapper-v0.1/dbVIR/virDB_B.csv', help='CSV database file [/usr/local/readmapper/dbVIR/virDB_B.csv]')
    parser.add_argument('-F', '--force', dest="force", action='store_true', default=False, help='Overwrite the clusters')
    parser.add_argument('-v', '--verbose',    dest="verbose",    default="0",help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version',    action='version', version='cluster_virDB-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
	run()