# diffexp.py

# This script is for identifying super-enhancer associated genes that are differentially expressed between two stages


from DYSE_main import formatFolder
import pandas as pd
import subprocess
import argparse


def diffexp(deFile, SEgenes):
    de_list = [x.rstrip().split('\t') for x in deFile]

    col = []
    for elem in de_list[0]:
        if elem != '':
            col.append(elem)

    diffexp_df = pd.DataFrame(de_list[1:], columns=col)
    geneList = set(diffexp_df['gene_id'].tolist())

    st_diffexp_genes = list(geneList & set(SEgenes))

    wanted = diffexp_df.loc[diffexp_df['gene'].isin(st_diffexp_genes)]

    return wanted

def main():

    '''
    main run call
    '''

    usage = '%(prog)s  [options] -i [INPUT_FILES] -d [RNA-SEQ_DIFF_EXP_FILE] -o [OUTPUT_FOLDER]'
    parser = argparse.ArgumentParser(prog='DYSE_diffexp.py', usage=usage)
    # Required flags
    parser.add_argument("-i", "--i", dest="input", default=None,
                        help="Comma separated list of SEgene files")
    parser.add_argument("-d", "--diffexp", dest="deFile", default=None,
                        help="RNA-seq differential expression file that includes stages of interest")
    parser.add_argument("-o", "--out", dest="out", default=None,
                        help="Output folder")

    # RETRIEVING FLAGS
    options = parser.parse_args()

    if not options.input or not options.deFile or not options.out:
        print("Hi there\nYour code seems to be missing some arguments")
        parser.print_help()
        exit()

    out_dir = formatFolder(options.out, True)

    inputFiles = options.input.split(',')

    deFile = open(options.deFile).read().rstrip('\n').split('\n')

    for stage in inputFiles:
        SEgenes = [item.split('\t')[0] for item in open(stage).read().rstrip('\n').split('\n')]
        diffexpSE = diffexp(deFile, SEgenes)
        temp = pd.DataFrame(columns=list(diffexpSE))
        last_col = []

        for index,row in diffexpSE.iterrows():
            if row['significant'] == 'yes' and float(row['log2(fold_change)']) > 0:
                last_col.append('upreg in ' + row['sample_2'])
                temp = temp.append(row, ignore_index=True)
            elif row['significant'] == 'yes' and float(row['log2(fold_change)']) < 0:
                last_col.append('downreg in ' + row['sample_2'])
                temp = temp.append(row, ignore_index=True)

        last_col_df = pd.DataFrame({'description': last_col})
        tofile = pd.concat([temp, last_col_df], axis=1, ignore_index=True)
        tofile.columns =  list(temp)+list(last_col_df)

        fname = out_dir+stage.split('/')[-1].split('.')[0]+'_SEgenes_diffexp.xls'
        subprocess.call(['touch', fname])

        pd.DataFrame.to_csv(tofile, path_or_buf=fname, sep='\t', header=True, index=False, line_terminator='\n')


if __name__ == "__main__":
    main()