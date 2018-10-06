import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats
import math

#outlier detection
from sklearn.cluster import DBSCAN
from sklearn import mixture

totalReads = 100000

#--------------------------------------------------------------------------------------------------------------------------------------
def KLdivergence( mu1, sigma1, mu2, sigma2 ):

	out = math.log(sigma2 / sigma1) + (sigma1**2 + (mu1 - mu2)**2)/(2.0*sigma2**2) - 0.5



#--------------------------------------------------------------------------------------------------------------------------------------
def displayProgress(current, total):

	barWidth = 70
	progress = float(current)/float(total)

	if progress <= 1.0:
		sys.stdout.write('[')
		pos = int(barWidth*progress)
		for i in range(barWidth):
			if i < pos:
				sys.stdout.write('=')
			elif i == pos:
				sys.stdout.write('>')
			else:
				sys.stdout.write(' ')
		sys.stdout.write('] '+str(int(progress*100))+' %\r')
		sys.stdout.flush()
#--------------------------------------------------------------------------------------------------------------------------------------
#eventalign output
sixmer2eventsBrdU = {}
f = open(sys.argv[1],'r')
currentRead = ''
readCounter = 0
for line in f:

	splitLine = line.rstrip().split('\t')

	#ignore the header line
	if splitLine[0] == 'contig':
		continue

	readIndex = splitLine[3]
	if readIndex != currentRead:
		currentRead = readIndex
		readCounter += 1
		displayProgress(readCounter, totalReads)

	if readCounter == totalReads:
		break

	eventTime = float(splitLine[8])
	if eventTime < 0.002:
		continue

	sixmer = splitLine[9]
	if sixmer not in sixmer2eventsBrdU:
		sixmer2eventsBrdU[sixmer] = []
	else:
		sixmer2eventsBrdU[sixmer].append( float(splitLine[6]) )
f.close()

#pore model
model = {}
f = open(sys.argv[2],'r')
for line in f:

	if line[0] == '#' or line[0] == 'k':
		continue
	
	splitLine = line.rstrip().split('\t')
	model[splitLine[0]] = [	float(splitLine[1]), float(splitLine[2]) ]
f.close()

for i, key in enumerate(sixmer2eventsBrdU):

	if key == 'NNNNNN':
		continue
	
	if len( sixmer2eventsBrdU[key] ) > 200:
		x = np.linspace( np.mean(sixmer2eventsBrdU[key])-15, np.mean(sixmer2eventsBrdU[key])+15, 1000 )
		plt.figure(i)

		#noise reduction
		ar = np.array(sixmer2eventsBrdU[key])
		db = DBSCAN( min_samples= (0.025*len( sixmer2eventsBrdU[key] )) ).fit(ar.reshape(-1,1))
		outliers_filtered = []
		for j, label in enumerate(db.labels_):
			if label == -1:
				continue
			else:
				outliers_filtered.append(sixmer2eventsBrdU[key][j])

		if len(outliers_filtered) > 50:
			densityFiltered = stats.kde.gaussian_kde( outliers_filtered )
			plt.plot(x, densityFiltered(x), label='Filtered BrdU Density, N='+str(len(outliers_filtered)))

			gmm = mixture.GMM(n_components=2, covariance_type='full')
			ar_filtered = np.array(outliers_filtered)
			out = gmm.fit(ar_filtered.reshape(-1,1))
			yMix = mlab.normpdf( x, out.means_[0], math.sqrt(out.covars_[0]) )
			plt.plot( x, yMix, label='Fit Distribution (1)')
			yMix = mlab.normpdf( x, out.means_[1], math.sqrt(out.covars_[1]) )
			plt.plot( x, yMix, label='Fit Distribution (2)')

			if abs(model[key][0] - out.means_[0]) < abs(model[key][0] - out.means_[1]):
				print key + '\t' + str(out.means_[1][0]) + '\t' + str(math.sqrt(out.covars_[1])) + '\t' + str(out.means_[0][0]) + '\t' + str(math.sqrt(out.covars_[0])) + '\t' + str(KLdivergence( out.means_[1][0], math.sqrt(out.covars_[1]), out.means_[0][0], math.sqrt(out.covars_[0]) ))
			else:
				print key + '\t' + str(out.means_[0][0]) + '\t' + str(math.sqrt(out.covars_[0])) + '\t' + str(out.means_[1][0]) + '\t' + str(math.sqrt(out.covars_[1])) + '\t' + str(KLdivergence( out.means_[0][0], math.sqrt(out.covars_[0]), out.means_[1][0], math.sqrt(out.covars_[1]) ))

		#density
		densityBrdU = stats.kde.gaussian_kde( sixmer2eventsBrdU[key] )
		plt.plot(x, densityBrdU(x), label='BrdU Density')
		
		#pore model
		yModel = mlab.normpdf( x, model[key][0], model[key][1] )
		plt.plot(x, yModel, label='6mer Pore Model')

		#plotting stuff
		plt.xlabel('pA')
		plt.ylabel('Count')
		plt.title( key + '  N=' + str(len(sixmer2eventsBrdU[key])) )
		plt.legend(loc='upper right', framealpha=0.5)
		plt.savefig( key + '.pdf' )
		plt.close()
