""" Converts PDF graphics to bitmap png. """
import os
import sys


def main(prefix, N):
    fin_suffix = '.pdf'
    for i in range(N):
        fin_name = prefix + str(i) + fin_suffix
        fout_name = prefix + str(i)
        command = 'pdftoppm ' + fin_name + ' ' + fout_name + ' -png'
        print(command)


if __name__ == '__main__':
    if (len(sys.argv) != 3):
        print('Usage: python '+sys.argv[0]+' file_prefix NFILES')
        sys.exit()
    else:
        main(sys.argv[1], int(sys.argv[2]))
