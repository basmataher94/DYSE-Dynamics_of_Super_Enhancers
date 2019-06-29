# compare.py

# script to identify changes in super-enhancer activity between different stages in a given cell lineage
# or between different cell conditions


import subprocess
import DYSE_main
import argparse


def remodelling(stage, histone, ht_dir):
    rep = []
    act = []

    for i, j in zip(range(0,len(stage)-1), range(1,len(stage))):
        if not (histone[i]=='None' or histone[j]=='None'):
            st1name = stage[i].split('/')[-1].split('.')[0]
            st2name = stage[j].split('/')[-1].split('.')[0]
            rep_file = ht_dir+'repressed_SEs_'+st1name+'_to_'+st2name+'.bed'
            act_file = ht_dir+'activated_SEs_'+st1name+'_to_'+st2name+'.bed'

            rep.append(rep_file)
            act.append(act_file)

            command1 = 'bash -c \"bedtools intersect -a '+stage[i]+' -b '+histone[j]+' -wo > '+rep_file + '\"'
            subprocess.call(command1, shell=True)
            command2 = 'bash -c \"bedtools intersect -a '+stage[j]+' -b '+histone[i]+' -wo > '+act_file + '\"'
            subprocess.call(command2, shell=True)

    return [rep, act]


def main():
    usage = '%(prog)s  [options] -i [INPUT_FILES] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]'
    parser = argparse.ArgumentParser(prog='DYSE_compare.py', usage=usage)
    # required flags
    parser.add_argument("-i", "--i", dest="input", default=None,
                        help="Comma separated list of SE bed files from different stages/conditions")
    parser.add_argument("-hs", "--histone", dest="histone", default=None,
                        help="Comma separated list of histone H3K27me3 bed files")
    parser.add_argument("-o", "--out", dest="out", default=None,
                        help="Output folder")

    # RETRIEVING FLAGS
    options = parser.parse_args()

    if not options.input or not options.out or not options.histone:
        print("Hi there\nYour code seems to be missing some arguments")
        parser.print_help()
        exit()

    inputFiles = options.input.split(',')
    histoneFiles = options.histone.split(',')
    ht_dir = DYSE_main.formatFolder(options.out + 'SE_changes/', True)

    rep_files, act_files = remodelling(inputFiles, histoneFiles, ht_dir)


if __name__ == "__main__":
    main()