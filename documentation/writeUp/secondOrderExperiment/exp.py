#let's do a small clustering experiment
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
np.random.seed(0)

#class for a clustered point

class Point:
	instances = []
	#xmin, xmax = 0, 1
	#ymin, ymax = 0, 1	
	ROIArea = 1
	def __init__(self, x, y):
		self.x = x
		self.y = y
		Point.instances.append(self) #adds to our list of instances
	def __str__(self):
		return f'p({self.x},{self.y})'
	def __repr__(self):
		return str(self)
	def dist2(self,other):
		return (self.x - other.x)**2 + (self.y - other.y)**2
	def dIJ(self, r):
		return sum(0 < self.dist2(other) <= r**2 for other in self.instances)
		
	def L(self, r):
		localPoints = self.dIJ(r)
		localPoints *= self.ROIArea / ( np.pi * ( len(self.instances) - 1 ) )
		return np.sqrt(localPoints)
	def scatter(self, fName=None):
		x,y=[],[]
		for p in self.instances:
			x.append(p.x)
			y.append(p.y)
			if id(p) == id(self): continue
		x.append(self.x)
		y.append(self.y)
		s, c = [0.5]*len(y), [1]*len(y)
		s[-1], c[-1] = 100, 44
		plt.figure(figsize=(15,7.5), dpi=80)
		plt.scatter(x,y, s = s)
		plt.savefig('bigscatter,png')

	def reset(self):
		self.instances = [self]
		
def printer(x,y,p1, *args):
	plt.clf()
	fNames = [arg for arg in args]
	p1.scatter(fNames[0])
	
	plt.figure(figsize=(15,7.5), dpi=80)

	plt.subplot(1,2, 1)
	plt.plot(x,y, label = 'L(r)')
	plt.xlabel('r')
	plt.ylabel('L(r)')
	plt.legend()
	plt.subplot(1,2,2)
	plt.plot(x,y, label = 'L(r)')
	plt.plot(x, np.linspace(0,0.5,500), label='linear best fit')
	plt.xlabel('r')
	plt.ylabel('L(r)')
	plt.legend()
	plt.show()
		
		
if __name__ == '__main__':
	#we don't need to keep track of them 
	points = []
	for x, y in zip((0.6 + 0.1*np.random.rand(100)), (0.6 + 0.1*np.random.rand(100))):
		#call constructor
		points.append(Point(x,y))
		
	#another cluster in the corner
	for x, y in zip(( 0.1*np.random.rand(100)), (0.3 + 0.1*np.random.rand(100))):
		#call constructor
		points.append(Point(x,y))
	for x, y in zip((0.5+ 0.1*np.random.rand(100)), (0.1 + 0.1*np.random.rand(100))):
		#call constructor
		points.append(Point(x,y))
	p1 = Point(0.5,0.5)
	p1.scatter()
	
	varmax = -99999.0
	rbest = 0
	x = np.linspace(0, 2,50)
	variances = []
	for r in tqdm(x):
		lscoresVar = np.var([p.L(r) for p in Point.instances])
		variances.append(lscoresVar)
		if lscoresVar > varmax :
			varmax = lscoresVar
			rbest = r
	print(rbest, varmax)

	



	
