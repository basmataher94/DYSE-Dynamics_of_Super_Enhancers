'''
ALGORITHM TO COMPARE SUPER-ENHANCERS BETWEEN DIFFERENT STAGES IN A GIVEN CELL LINEAGE
OR BETWEEN DIFFERENT CELL CONDITIONS

VERSION 1.0  2018

CONTACT: tperkins@ohri.ca
'''


import os
import re
import subprocess
import argparse
import time
import DYSE_compare

from matplotlib_venn import venn2
from matplotlib import pyplot as plt
plt.switch_backend('agg')

from operator import itemgetter
from joblib import Parallel, delayed

'''
======================================================
CALLING ROSE SCRIPTS
======================================================
'''


def callRose(fn, rose_out, g, s=None):
    fileNames = ''
    [i,r,c] = fn.split(',')
    stitch = ''

    if s is not None:
        stitch += ' -s '+str(s)
    fname = formatFolder(i.split('/')[-1].split('.')[0], True)

    if c == 'None':
        command = ['-i '+i, ' -r '+r, ' -o '+rose_out+fname, ' -g '+g, stitch, ' -t 2500']
    else:
        command = ['-i '+i, ' -r '+r, ' -o '+rose_out+fname, ' -g '+g, ' -c '+c, stitch, ' -t 2500']
    print command

    command_1 = 'python ROSE_main.py '+''.join(command)
    print command_1
    print rose_out+fname

    subprocess.call(command_1, shell=True)
    fileNames += rose_out+fname+fname[:-1]+'_Gateway_SuperEnhancers.bed\n'
    return fileNames


def callGeneMapper(fnames, out_dir, g):
    for fn in fnames:
        command = ['-i ' + fn, '-o ' + out_dir, '-g ' + g]
        command_1 = 'python ROSE_geneMapper.py ' + ' '.join(command)
        print command_1
        subprocess.call(command_1, shell=True)


'''
========================================================
SUPER-ENHANCER FUNCTIONS
========================================================
'''


# Takes in SE files for all TFs corresponding to a certain stage
# and groups the together, combining overlaps
def groupSE(stageFile):
    filesSE = stageFile.rstrip('\n').split('\n')
    print filesSE
    allSE = []
    for fn in filesSE:
        f = open(fn).read().rstrip('\n').split('\n')
        for line in f:
            if not (line[0].startswith('#') or line[0] is 'REGION'):
                temp = line.rstrip('\n').split('\t')
                temp[1] = int(temp[1])
                temp[2] = int(temp[2])
                allSE.append(temp)

    sorted_allSE = sorted(allSE, key=itemgetter(0,1))
    if len(filesSE) > 1:
        return mergeOverlaps(sorted_allSE)
    else:
        return sorted_allSE


# Merge overlapping super-enhancers within the same stage
def mergeOverlaps(toMerge):
    for reg1,reg2 in zip(toMerge[:-1], toMerge[1:]):
        if reg1 == reg2:
            toMerge.remove(reg2)
        elif reg1[0] == reg2[0] and reg1[2] > reg2[1]:
            reg1[2] = reg2[2]
            for x in [3,4]:
                reg1[x] += ', '+reg2[x]
            toMerge.remove(reg2)
    return toMerge


# Compares the super-enhancers from different stages
def compareStages(stage1, stage2, outdir):
    s1 = stage1.split('/')[-1].split('.')[0]
    s2 = stage2.split('/')[-1].split('.')[0]

    outfile = outdir + 'Intersect_' + s1 + '_' + s2 + '.bed'
    command = 'bash -c \"bedtools intersect -a '+stage1+' -b '+stage2+' -wo > '+outfile+'\"'
    subprocess.call(command, shell=True)

    return [outfile, (s1,s2)]


# Visualizing results
def visual(ov_data, venn_out):
    for item in ov_data:
        data = item[0]
        names = item[1]

        plt.clf()
        venn2(subsets=data, set_labels=names)
        plt.title('Overlap between %s and %s' % (names[0],names[1]))
        plt.savefig(venn_out+names[0]+'_'+names[1]+'.png')


'''
========================================================
FILE FORMATTING FUNCTIONS
========================================================
'''


# Writes output to file
def formatOutFile(fname, data, plus=False):
    if plus:
        f = open(fname, 'w+')
    else:
        f = open(fname, 'w')
    for item in data:
        f.write(item + '\n')
    f.close()


# For formatting the overlap files
def formatIntersectFile(fname):
    result = [sub.replace('\r', '') for sub in open(fname).read().split('\n')[:-1]]
    formatOutFile(fname, result)


# Creates folders for output files
def formatFolder(dir_name, make=False):
    if dir_name[-1] != '/':
        dir_name += '/'
    if make:
        if not os.path.isdir(dir_name):
            subprocess.call(['mkdir', dir_name])
    return dir_name


# Extracts data from file
def fileToList(fname):
    return re.split('\n|\r\t', open(fname).read())


def main():
    '''
    main run call
    '''

    usage = '%(prog)s  [options] -i [INPUT_FILES] -o [OUTPUT_FOLDER] -t [ANALYSIS_TYPE] -g [GENOME] [OPTIONAL_FLAGS]'
    parser = argparse.ArgumentParser(prog='DYSE_main.py', usage=usage)
    # Required flags
    parser.add_argument("-i", "--i", dest="input", default=None,
                        help="Enter a comma separated list of stage/condition files containing TF SE file names")
    parser.add_argument("-o", "--out", dest="out", default=None,
                        help="Enter an output folder")
    parser.add_argument("-t", "--type", dest="type", default=1, choices=['1','2'],
                        help="Enter 1 for adjacent stage comparisons only, or 2 for comparing all stages together")
    parser.add_argument("-g", "--genome", dest="genome", default=None,
                        help="Enter the genome build (MM10,MM9,MM8,HG18,HG19)")

    # Optional flags
    parser.add_argument("-r", "--rose", dest="rose", default=None,
                        help="ROSE output directory")
    parser.add_argument("-s","--stitch", dest="stitch", default=12500,
                        help="Enter a max linking distance for stitching")

    # RETRIEVING FLAGS
    options = parser.parse_args()

    if not options.input or not options.out or not options.type or not options.genome:
        print("Hi there\nYour code seems to be missing some arguments")
        parser.print_help()
        exit()

    out_dir = formatFolder(options.out, True)

    bed_dir = formatFolder(out_dir+'bed_files/', True)

    sp_dir = formatFolder(out_dir+'stage_specific/', True)

    if options.rose:
        rose_out = formatFolder(out_dir+options.rose, True)
    else:
        rose_out = None

    inputFiles = options.input.split(',')

    bed_files = []
    sp_files = []

    for fn in inputFiles:
        stageFile = open(fn)
        if options.rose:
            rosefiles = stageFile.read().rstrip('\n').split('\n')
            if options.stitch:
                result = Parallel(n_jobs=12)(delayed(callRose)(f,rose_out,options.genome,options.stitch) for f in rosefiles)
            else:
                result = Parallel(n_jobs=12)(delayed(callRose)(f,rose_out,options.genome) for f in rosefiles)

            while len(result) < len(rosefiles):
                time.sleep(600)
            temp = ''.join(result)
            compileStageSE = groupSE(temp)
        else:
            compileStageSE = groupSE(stageFile.read())
        stageFile.close()

        compileStageSE = ['\t'.join(map(str, reg)) for reg in compileStageSE]
        fname = bed_dir + fn.split('.')[0] + '.bed'
        sp_name = sp_dir + fn.split('.')[0] + '_specific_SE.bed'

        bed_files.append(fname)
        sp_files.append(sp_name)

        subprocess.call(['touch', fname])
        subprocess.call(['touch', sp_name])

        formatOutFile(fname, compileStageSE)

    all_name = bed_dir + 'all_reg.bed'
    subprocess.call(['touch', all_name])
    for bd in bed_files:
        name = bd.split('/')[-1].split('.')[0]
        print name
        lookup = [x for x in bed_files if x is not bd]
        print lookup
        for f in lookup:
            open(all_name, 'a+').write(open(f).read())
        tmp_command = 'bash -c \"bedtools sort -i ' + all_name + ' > ' + bed_dir + 'all_reg_sorted.bed\"'
        subprocess.call(tmp_command, shell=True)
        command = 'bash -c \"bedtools intersect -a ' + bd + ' -b ' + all_name + ' -v > ' + sp_dir + name + '_specific_SE.bed\"'
        print command
        subprocess.call(command, shell=True)
        subprocess.call('rm ' + all_name, shell=True)

    # for bd in bed_files:
    #     name = bd.split('/')[-1].split('.')[0]
    venn_out = formatFolder(out_dir+'figures/', True)
    ov_dir = formatFolder(out_dir+'overlaps/', True)

    venn_ov_data = []
    ov_filenames = []

    if options.type == str(2):
        for i in range(0,len(bed_files)-1):
        # for stage1 in bed_files[:-1]:
            a = i+1
            stage1 = bed_files[i]
            for j in range(a,len(bed_files)):
            # for stage2 in bed_files[a:]:
                stage2 = bed_files[j]
                print stage1
                print stage2

                # ov_filenames.append(compareStages(stage1, stage2, ov_dir)[0])
                [outf, st_names] = compareStages(stage1, stage2, ov_dir)
                ov_filenames.append(outf)
                tmp = fileToList(outf)

                total_st1 = open(stage1).read().rstrip('\n').split('\n')
                total_st2 = open(stage2).read().rstrip('\n').split('\n')

                int_set = len(tmp) - 1

                venn_ov_data.append([(len(total_st1) - int_set, len(total_st2) - int_set, int_set), st_names])
                # a = a+1
    else:
        for stage1, stage2 in zip(bed_files[:-1], bed_files[1:]):
            [outf,st_names] = compareStages(stage1, stage2, ov_dir)
            ov_filenames.append(outf)
            tmp = fileToList(outf)

            total_st1 = open(stage1).read().rstrip('\n').split('\n')
            total_st2 = open(stage2).read().rstrip('\n').split('\n')

            int_set = len(tmp)-1

            venn_ov_data.append([(len(total_st1)-int_set, len(total_st2)-int_set, int_set), st_names])

    print venn_ov_data
    visual(venn_ov_data, venn_out)

    for item in ov_filenames:
        formatIntersectFile(item)

    genes = formatFolder(out_dir+'geneMapper/', True)
    callGeneMapper(sp_files+ov_filenames, genes, options.genome)


if __name__ == "__main__":
    main()
