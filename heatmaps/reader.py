'''
this file will produce a heatmap of the cluster scores
provided with a 2d array saved to a text file
DELINEATED BY COMMAS
'''

import numpy as np
import matplotlib.pyplot as plt


def genNumbersFromLine(line, tbins=35):
    output = np.zeros(tbins)
    line = line.split(',')
    for i in range(tbins):
        output[i] = float(line[i])
    return output[::-1]

def readFromFile(fName, rbins=35, tbins=35):
    rbins, tbins = int(rbins), int(tbins)
    output = np.zeros((rbins, tbins))
    r = 0
    with open(fName) as f:
        line = 1
        while line:
            line = f.readline()
            if len(line) < 3: break
            output[r,] = genNumbersFromLine(line)
            r+=1
    return output

def mapper(arr):
    plt.imshow(arr, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    fName = os.path.dirname(os.path.realpath(__file__))+ "/heatmap.png"
    print('saving heatmap as',fName)
    plt.savefig(fName)

if __name__ == '__main__':
    import sys
    import os
    args = [arg for arg in sys.argv[1:]]
    try:
        arr = readFromFile(*args)
    except Exception as e:
        print(e)
        print('supply scores file, rbins, tbins to the command line')
        print('eg python3 reader.py scores.txt 35 35')
        exit()
    maxscore = np.max(arr)
    print(f'maxscore is {maxscore}')
    arr = np.where(arr > maxscore * 1.5, arr, maxscore * 1.5)
    mapper(arr)