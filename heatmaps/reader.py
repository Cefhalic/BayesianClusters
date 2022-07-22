import numpy as np
import matplotlib.pyplot as plt

def gett(line):
    index1 = line.index(':') + 1
    index2 = line.index(',')
    return float(line[index1:index2]), index2+1

def reader(fName):
    rDict = {}
    with open(fName) as f:
        line = 1
        while line:
            line = f.readline()
            if len(line) < 4: continue
            if line[4] != 'R': continue

            #read the R value
            R, index = gett(line)
            line = line[index:]

            #read T value
            T, index = gett(line)
            line = line[index:]

            #read Score
            score, _ = gett(line)

            if R not in rDict:
                rDict[R] = {T:score}
            rDict[R][T] = score
    return rDict

def creator(rDict):
    #dict which has r values as a key
    #each r value contains a dict that looks like
    #t:score within

    rWidth = len(rDict.keys())
    #assume it's the same for T for now
    output = np.zeros((rWidth, rWidth))

    for rIndex, rKey in enumerate(sorted(rDict.keys())):
        #get the current T dict 
        tDict = rDict[rKey]
        for tIndex, tKey in enumerate(sorted(tDict.keys())):
            output[rIndex, tIndex] = tDict[tKey]
    return output


if __name__ == '__main__':
    import sys
    rd = reader(sys.argv[1])
    arr = creator(rd)

    plt.imshow(arr.transpose(), cmap='hot')
    fName = 'heatmap.png'    
    try:
        fName = sys.argv[2]
        assert(fName[-4:] == '.png')
    except:
        pass
    plt.savefig(fName)
