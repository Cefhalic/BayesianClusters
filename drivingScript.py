import fileinput
import sys, os, shutil
import numpy as np
import subprocess
import signal

# break_sig = False

# def handler(signum, frame):
#   print( "Ctrl-c was pressed. Terminating ASAP.", flush=True)
#   global break_sig
#   break_sig = True

# signal.signal(signal.SIGINT, handler)


class Box:
    def __init__(self, n, e, s, w, scaleFactor = 250):
        self.n = n
        self.e = e
        self.s = s
        self.w = w
        self.area = abs(self.n - self.s) * abs(self.e - self.w)
        self.scaleFactor = scaleFactor
    def __str__(self):
        return f'Box: n-{self.n}, e-{self.e}, s-{self.s}, w-{self.w}, area-{self.area}'
    def getArea(self):
        return self.area
    # def __eq__(self, other):
    #     return self.area == other.area
    def __lt__(self,other):
        return self.area < other.area

    def centre(self):
        return ((self.w + self.e - 2) / 8.0 , (self.s + self.n - 2) / 8 )

    def zoom(self):
        return max(abs(self.s - self.n), abs(self.w - self.e)) / 4.0

#let's write the function that sets the scan areas

def writeLine(file, box):
    for line in fileinput.input(file, inplace = 1):
        if 'zoom' in line:
            line = f'--zoom {box.zoom():.1f}um\n'
        if 'centre' in line:
            centre = box.centre()
            line = f'--centre {centre[0]:.1f}um {centre[1]:.1f}um\n'
        sys.stdout.write(line)


def reader(fName = '../1_un_red.csv', scaleFactor=250):
    scaleFactor = 250

    xShape, yShape = 520, 520
    print(f'setting default shape of {xShape} x {yShape}')
    print('.'*50)
    output = np.zeros((xShape, yShape))

    with open(fName) as f:
        line = f.readline()
        while line:
            line = f.readline()
            try:
                line = line[line.index(',') + 1:] #skip id
                line = line[line.index(',') + 1:] #skip frame
            except:
                print(line)
                break
            #read x value
            x = float(line[:line.index(',')])
            line = line[line.index(',') + 1:]

            #read y value
            y = float(line[:line.index(',')])
            
            # output[int(x//500), int(y//500)] += 1
            output[int(x//scaleFactor), int(y//scaleFactor)] += 1

    #normalise
    output /= output.max()
    return output

def validIndices(i, j, array): #generator that yields "good" indices
    for lI in range(-1,2):
        if (i + lI < 0) or (i + lI >= array.shape[0]): continue
        for lJ in range(-1,2):
            if (j + lJ < 0) or (j + lJ >= array.shape[1]): continue
            if (lI == 0) and (lJ == 0): continue
            yield i + lI, j + lJ

def NESWsetter(x, y, n, e, s, w):
    if x > e: e = x #just seen the easternmost point
    if x < w: w = x #just seen westernmost point
    if y < n: n = y
    if y > s: s = y
    return n, e, s, w

def bfs(i, j, array, visited, componentID): 
    westPoint, eastPoint = array.shape[1], 0
    northPoint, southPoint = array.shape[0], 0
    stacks = [[]]
    stacks[0].append((i,j))
    visited[i, j] = 1
    array[i,j] = componentID
    counter = 0
    while stacks[counter]:
        stacks.append([]) #new empty one at location (i + 1)
        #next, loop over the indices in stacks[i]
        for index in stacks[counter]:
            #loop over the neighbours of this index
            for x, y in validIndices(*index, array):
                if visited[x,y]: continue
                else:
                    #set everything
                    if array[x,y] == 1: 
                        array[x,y] = componentID
                        stacks[counter+1].append((x,y))
                        #set the maxvals
                        northPoint, eastPoint, southPoint, westPoint = NESWsetter(x,y, northPoint,eastPoint, southPoint,westPoint)
                    visited[x,y] = 1           
        counter += 1
    if len(stacks) == 2 and len(stacks[1]) == 0:
        print('reached early exit in bfs')
        return array, visited, Box(0, 0, 0, 0)
    return array, visited, Box(northPoint, eastPoint, southPoint, westPoint)

def blurAndThresh(arr, rounds=15, threshold = 0.2):
    if arr.max() > 1.0:
        #normalise
        arr /= arr.max()

    pad_width = 2 * rounds
    pad_arr = np.zeros((arr.shape[0] + pad_width, arr.shape[1] + pad_width))
    pad_arr[pad_width//2:-pad_width//2, pad_width//2:-pad_width//2] = arr
    arr = pad_arr
    arr_empty = np.zeros_like(arr)

    #blurring
    for i in range(-rounds, rounds + 1):
        for j in range(-rounds, rounds + 1):
            arr_temp = np.roll(arr, i, 0)
            arr_temp = np.roll(arr_temp, j, 1) 
            arr_empty += arr_temp
    #renormalise
    arr_empty /= arr_empty.max()

    #threshold
    arr_empty = np.where(arr_empty > threshold, 1, 0)

    #remove the padding, return it
    return arr_empty[pad_width//2:-pad_width//2, pad_width//2:-pad_width//2]
    
def roiID(arr):
    visited = np.zeros_like(arr)
    ComponentID = 2
    boxes = []
    #let's hit every pixel
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if not visited[i,j]:
                if arr[i,j] == 0.0:
                    visited[i,j] == 1
                    continue
                #retrieve the box
                arr, visited, box = bfs(i, j, arr, visited, ComponentID)
                ComponentID += 1
                if box.getArea() < 500: continue
                boxes.append(box)
    return sorted(boxes) #sorts them by area

if __name__ == '__main__':
    #read filename from command line
    try:
        args = sys.argv[1:]
        assert any(arg.endswith('txt') for arg in args)
        assert any(arg.endswith('csv') for arg in args)
        for arg in args:
            if arg.endswith('csv'):
                fName = arg
            if arg.endswith('txt'):
                configFile = arg
    except:
        print('supply a thunderstorm csv, and a config file to modify')
        exit()
    #read it
    image = reader(fName)

    #blur and thresh
    image = blurAndThresh(image)

    scanAreas = roiID(image) #returns the boxes sorted by area

    PATH = os.path.dirname(os.path.realpath(__file__)) #for now, just put it next to this python script
    lFName = PATH + '/tempconfig.txt'

    print('-'*50)
    print(f'we found {len(scanAreas)} boxes')
    print('-'*50)


    scanMore = True
    for scanArea in scanAreas[:10]:
        if not scanMore: break
        #make a copy of the config file we found
        #assume it doesn't have path attached
        print(f'centre: {scanArea.centre()}, zoom: {scanArea.zoom()}')

        outfile = PATH + f'/results_SA_{int(scanArea.centre()[0])}_{int(scanArea.centre()[1])}_Z_{int(scanArea.zoom())}.json'
        print('sending to', outfile)
        
        #make a new copy of the config file
        shutil.copyfile(configFile, lFName) 
        #write the centre, zoom to this file
        writeLine(lFName, scanArea)      

        print('running a scan')
        cmdLineArg = f'{PATH}/Cluster.exe -i {PATH}/1_un_red.csv --cfg {lFName} -o {outfile}'
        # print('running', cmdLineArg)
        # exit()
        # os.system(cmdLineArg)
        processes = []
        try:
            process = subprocess.Popen([cmdLineArg],  shell=True)
            processes.append(process)
            print('running with', process.pid)
            process.wait()
        except KeyboardInterrupt:
            print('\nctrl c pressed')
            for process in processes: 
                print('killing process', process.pid)
                process.kill()
            break
    print('exit')
    