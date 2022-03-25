'''
This assignment perfoms pairwise alignment to protein sequences contained in
seq1.txt and seq2.txt

For more information about the biopython package
https://biopython.org/docs/1.75/api/Bio.pairwise2.html

Please uncomment the area where it says PART1 or PART2 for the corresponding
alignment parameters
'''

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import os

def create_lists(): # test reads seq1.txt with a relative filepath
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, 'seq1.txt')

    with open(filename, "r") as f:
        next(f)             # skip first line in file
        list1 = [line.strip() for line in f if line.strip()]

    filename = os.path.join(dir, 'seq2.txt')

    with open(filename, "r") as f:
        next(f)             # skip first line in file
        list2 = [line.strip() for line in f if line.strip()]

    '''
    for i in list1:
        print(i)

    for i in list2:
        print(i)
    '''
    # return the lists as strings without newline
    return ''.join(list1).replace("\n", ""), ''.join(list2).replace("\n", "")

if __name__ == '__main__':
    seq1, seq2 = create_lists()
    matrix = matlist.blosum62

    # PART 1
    print("************* Global Alignment *************")
    for x in pairwise2.align.globaldx(seq1, seq2, matrix):
        print(format_alignment(*x))

    print("************* Local Alignment *************")
    for x in pairwise2.align.localdx(seq1, seq2, matrix):
        print(format_alignment(*x))
    
    '''
    # PART 2
    print("************* Global Alignment *************")
    for x in pairwise2.align.globalds(seq1, seq2, matrix, -2, -0.5):
        print(format_alignment(*x))

    print("************* Local Alignment *************")
    for x in pairwise2.align.localds(seq1, seq2, matrix, -2, -0.5):
        print(format_alignment(*x))
    '''
