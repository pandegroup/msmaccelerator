#===============================================================================
# ssaCalculator.py
#
# Please reference
# N Singhal Hinrichs, and VS Pande. J. Chem. Phys. 2007. Calculation of the 
# distribution of eigenvalues and eigenvectors in Markovian State Models for
# molecular dynamics.
#
# Written 10/1/10 by
# Dan Ensign <densign@mail.utexas.edu>
# Gregory Bowman <gregoryrbowman@gmail.com>
# Sergio Bacallado <sergiobacallado@gmail.com>
# Stanford University
# Pande group
#
# Copyright (C) 2008  Stanford University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


from random import gammavariate, random
from copy import copy

from math import log, fabs
from numpy import real, round, float64, common_type
from os.path import exists
from os import mkdir
from pickle import Pickler
from scipy.linalg import eigvals, lu, solve, LinAlgError
from scipy import argmax, diag, identity, int32, matrix, zeros
import scipy.sparse.linalg.eigen
from scipy.io import mmread, mmwrite


class ssaCalculator( object ):
        """ Class for objects to govern sensitivity calculations. """
        def __init__( self, lagTime, tcMatrix, priorCounts=.5, evalList=[1], nNewSamples=1, timeUnits="ns" ):
                self.lagTime = lagTime          # a number in appropriate units 

                self.tcMatrix = tcMatrix        # a regular list of lists of numbers, square

                self.priorCounts = priorCounts  # this integer is the uniform alpha parameter of the Dirichlet prior
						# which must be greater than 0. By default, it is 0.5

                self.evalList = evalList	# the list of eigenvalues (as arguments of a sorted list from largest to
						# smallest) to be used in sensitivity computations, by default
						# only the second largest eigenvalue is used, [1].

                self.nNewSamples = nNewSamples  # the number of new samples that we plan to simulate.

		self.timeUnits = timeUnits	

                self.evalList.sort()

		# derived properties
                self.dim = len( self.tcMatrix )
		self.weights = self.computeWeights()	# the total number of transitions observed yet from each state
		self.varianceContributions = []		# a list of sundry results
                self.compute()

	def computeWeights( self ):
		# sum the rows of the tcMatrix and return as a matrix
		weights = [] 
		for i in range( self.dim ):
			row = self.tcMatrix[ i ]
			rowsum = float( matrix( row ).sum() )
			weights.append( rowsum )
		return matrix( weights )

        def compute( self ):
                # 1. generate the Dirichlet matrix 
                self.dirmat = self.computeDirichletMatrix()
		del self.tcMatrix
		

                # 2. compute the eigenvalues
                self.avgEvals = list( eigvals( self.dirmat ) )
		# Sergio: take only the real part of the eigenvalues, such that the program
		# doesnt give weird results in computing the LU decomposition of A
                self.avgEvals = list(real(self.avgEvals))
		#self.avgEvals = list( arpack.eigen( self.dirmat, 1+max(self.evalList) )[0] )
		self.avgEvals.sort()
		self.avgEvals.reverse() # longest time scales first
		f = open("eigVals.dat", 'a')
		for ind in range(len(self.avgEvals)):
		  f.write(str(float(self.avgEvals[ind]))+" ")
		f.write("\n")
		f.close()

                for evalIndex in self.evalList:
                        eval = self.avgEvals[ evalIndex ]
			
                        # 3. compute the characteristic polynomial-generating matrix for the average matrix
#                        A = self.computeCharPolGenerator( eval )
			# to save memory will modify dirmat and then fix later
			for i in range(self.dim):
				self.dirmat[i,i] -= eval

                        # 4. Decompose the characteristic polynomial generator
                        ( P, L, U, Q ) = self.computeDecomposition( self.dirmat )

			# fix dirmat
			for i in range(self.dim):
				self.dirmat[i,i] += eval
			
			# 5. solve for the projection vectors x1 and x2
			( x1, x2 ) = self.computeProjectionVectors( P, L, U )
			del P
			del L
			del U

			# 6. compute the matrix of first derivatives
			normalization = x2.T * x1
			myNumerator = Q * x2 * x1.T * Q 
			del x1
			del x2
			firstDerivatives = myNumerator / normalization
			print "myNumerator", myNumerator
			print "normalization", normalization
			if normalization==0: 
				print "ALERT"
			del normalization
			del myNumerator
			del Q

			# 7. compute the contributions of each state to the eigenvalue uncertainty
			# and the reductions given 1 additional sample and m additional samples
			qlist = self.computeVarianceContributions( firstDerivatives )
			wlist = self.weights

			# Added by Sergio June/10
			# Store raw qlist in object
			self.qlist = qlist

			self.computeReductions( qlist, wlist )

        def computeDirichletMatrix( self ):
		dirmat = zeros([self.dim,self.dim], float64)
		i = 0
		while i < self.dim :
			
			# big note = if we didn't observe the state, then we assume that the transition
        		# probability of the state to itself is 1.0 -- therefore it won't show up as
        		# influencing the uncertainty. We also do not include the prior counts for such
			# a state.
			"""
			if self.weights[0,i] == 0 :
				dirmat[ i, i ] = 1
				
			else:
				row = arr2lst( self.tcMatrix[i] + self.priorCounts )
				frow = DirichletDistribution( row )

				for g in range( len( row ) ) :
					dirmat[i,g] = frow.mean( g )
			"""
			# Sergio: I am removing Dan's modification to the algorithm for states that 
			# haven't been observed.
			row = arr2lst( self.tcMatrix[i] + self.priorCounts )
			
			# Sergio: I'm trying a new thing here: adding 2 prior counts to self-transition
			# probabilities to try to avoid the possible effect of "traps"
			# row[i] += 2

			frow = DirichletDistribution( row )

			for g in range( len( row ) ) :
				dirmat[i,g] = frow.mean( g )

			i += 1


		return dirmat

        def computeCharPolGenerator( self, eval ):
		A = self.dirmat.copy()
		for i in range(self.dim):
			A[i,i] -= eval
                return A

        def computeDecomposition( self, A ):
		# first decomposition
		( P, L, U ) = lu( A )
	
		smallestIndex = findsmallestdiag( U )
		smallest = U[ smallestIndex, smallestIndex ]

		Q = identity( self.dim, dtype=float64 )
		if smallestIndex+1 != self.dim :
			del P
			del L
			del U

			# exchange smallestIndex row with dim-1 row
			# multiplying A by this matrix on both sides (Q.A.Q) will ensure that the
			# smallest element of U will be in the lower right corner.
			swaprow( Q, smallestIndex, self.dim-1 )
			
			# recompute the decomposition
			Q = matrix( Q )
			A = matrix( A )
			( P, L, U ) = lu( Q*A*Q )
	
		return ( P, L, U, Q )

	def computeProjectionVectors( self, P, L, U ) :	
		eK = matrix( identity( self.dim, float64 )[ 0: ,( self.dim - 1 ) ] ).T
		U = matrix(U, float64)
		U[ self.dim - 1, self.dim - 1 ] = 1.0
		# Sergio: I added this exception because in rare cases, the matrix
		# U is singular, which gives rise to a LinAlgError.
		try: 
			x1 = matrix( solve( U, eK ), float64 )
		except LinAlgError:
			print "Matrix U was singular, so we input a fake x1\n"
			print "U: ", U
			x1 = matrix(ones(self.dim))

		#print "x1", x1
		del U

		LT = matrix( L, float64, copy=False ).T
		PT = matrix( P, float64, copy=False ).T

		x2 = matrix( solve( LT*PT, eK ), float64 )
		del L
		del P
		del LT
		del PT
		del eK

		return ( x1, x2 )

	def computeVarianceContributions( self, firstDerivatives ) :

		qlist = []
		i = 0
		while i < self.dim :
			#print "Compute var %d" % i
			# derivatives of the eigenvalue for this state
			s = matrix( firstDerivatives[ i,0: ], float64 ).T
			##print s.shape

			# cross probability matrix for this state
			Pi = matrix( self.dirmat[ i,0: ], float64 ).T
			##print Pi.shape
			part1 = diag( arr2lst( Pi.T ) )
			part1 = matrix(part1, float64)
			##print part1.shape
			##print common_type(part1)
			Cp = matrix( part1 - Pi * Pi.T )
			##print common_type(Cp)
			##print Cp.shape
			del part1

			# degree of sensitivity for this state
			q = float( abs( s.T * Cp * s ) )
			del s
			del Pi
			del Cp

			qlist.append( q )

			i += 1

		return matrix( qlist )	

	def computeReductions( self, qlist, wlist ) :

		deltalist = qlist / ( 1. + self.weights )

		deltalist_1add = qlist / ( 2. + self.weights )
		reductionlist = deltalist - deltalist_1add
		fractionlist = reductionlist/reductionlist.sum()
		recclist = matrix( identity( self.dim )[ argmax( fractionlist ),0: ], "int" )

		deltalist_1add = arr2lst( deltalist_1add )
		reductionlist = arr2lst( reductionlist )
		fractionlist = arr2lst( fractionlist )
		recclist = arr2lst( recclist )

		currentEvalVariance = [ deltalist_1add, reductionlist, fractionlist, recclist ]

 		if self.nNewSamples > 1 :
			deltalist_madd = qlist / ( 1. + self.nNewSamples + self.weights ) 
			reductionlist_m = deltalist - deltalist_madd
			self.stateToSampleMore = argmax(reductionlist_m)
			fractionlist_m = reductionlist_m / reductionlist_m.sum()
			recclist_m = matrix( round( self.nNewSamples* fractionlist_m ), "int" )
			
			deltalist_madd = arr2lst( deltalist_madd )
			reductionlist_m = arr2lst( reductionlist_m )
			fractionlist_m  = arr2lst( fractionlist_m )
			recclist_m = arr2lst( recclist_m )

			currentEvalVariance.extend( [ deltalist_madd, reductionlist_m, fractionlist_m, recclist_m] ) 
	
		currentEvalVariance.insert( 0, arr2lst( deltalist ) )
		currentEvalVariance.insert( 0, range( self.dim ) )

		currentEvalVariance = matrix( currentEvalVariance, "float64" ).T   
	
		self.varianceContributions.append( currentEvalVariance )

	def displayContributions( self, bling = False ):
		if bling:
			sep1 = ":"
			sep2 = "|" 
		else:
			sep1 = " "
			sep2 = " "

		i = 0
		for evalIndex in self.evalList :

			eval = self.avgEvals[ evalIndex ] 

			if bling:
				s1 = "eigenvalue %6.6f" % eval

				try:
					timescale = -self.lagTime/log( eval )
					s1 += "   time scale %7.4f %s" % ( timescale, self.timeUnits )

				except ZeroDivisionError :
					s1 += "   time scale Inf %s" % self.timeUnits 

				except ValueError :
					pass
				
				#print s1
				#print "-"*len( s1 ) 

			# here are the results for this eigenvalue
			s2 = ""
			contributions = self.varianceContributions[ i ] 

			frm = "%6.6e " * 3
			v = "%d" % ( int( log( self.nNewSamples ) + 1 ) ) 
			frm += "%" + v + "d" 
	
			for j in range( contributions.shape[0] ):
				row = tuple( arr2lst( contributions[ j,0: ] ) )

				rowfrm = "%3d " + sep1 + " %6.6e " + sep2 + " \t" + frm
				if self.nNewSamples  > 1 :
					rowfrm += "\t "+ sep2 + " " + frm
				#print rowfrm % row

			if bling :	
				print
			i += 1

	def resamplingDistribution( self, eigIndex ):
		# just return the list of suggested new samples, to reduce the variance of that eigenvalue
		evalVarianceData = self.varianceContributions[ eigIndex ].T
		
		# take the last row
		dist = int32( evalVarianceData[ -1,0: ] )
		return tuple( arr2lst( dist ) )		

	def save( self, filename ):
		FILE = open( "%s/%s" % ( self.outputDir, filename ), "w" )
		p = Pickler( FILE )
		p.dump( self )
		FILE.close()	
		return True	


def openmat( filename ):
	FILE = open( filename )
	text = FILE.readlines()
	FILE.close()
	elems = []
	for line in text :
		row = []
		vals = line.split()
		for val in vals :
			row.append( int( val ) )

		elems.append( tuple( row ) )

	return tuple( elems )

def show( arr, label = "array", format = "% .6e " ):
        #format = "% 3.8f "
        nindex = arr.shape[0]
        s = "%s =\t" % label
        for i in range( nindex ):
                for j in range( nindex ):
                        elem = arr[i,j]
                        s += format % elem
                s += "\n\t"
        print s

def mathform( arr, format="%10.10f" ):
	s = "{ "
	arr = arr2lst( arr )

	try:
		dim = len( arr )
	except:
		dim = 0 

	if dim == 0 : 
		# just a number
		return format % arr

	elif dim == 1 :
		# either a 1-d array
		for i in range( dim ):
			s += format % arr[ i ]

			if i == dim-1 :
				s += " }"
			else:
				s += ", "
	else:
		for i in range( dim ):
			s += mathform( arr[ i ] ) 
			if i == dim-1 :
				s += " }"
			else:
				s += ", "

        return s

def makemat( dim, r=10, offset=0 ):
        a = zeros( (dim, dim), "float64" )
        for i in range( dim ):
                for j in range( dim ):
                        n = float( int( r*random() - offset ) )
                        a[i][j] = n
        return a

def sumNormalize( u ):
	v = copy( u )
	tot = float( u.sum() )
	return v/tot

def normalizeRows( mat ):
        # divide each row by the sum of the entries.
        dim = mat.shape[0]
        newmat = []
        i = 0
        while i < dim :
                row = mat[ i,0: ]
                newmat.append( arr2lst( sumNormalize( row ) ) )
                i += 1
	newmat = matrix( newmat )
        return newmat

def swapcol( mat, x, y ):
        # this will actually modify mat
        newcolx = copy( mat[ 0:,y ] )
        newcoly = copy( mat[ 0:,x ] )
        mat[0:,x] = newcolx
        mat[0:,y] = newcoly

def swaprow( mat, x, y ):
        # this will actually modify mat 
        newrowx = copy( mat[ y,0: ] )
        newrowy = copy( mat[ x,0: ] )
        mat[ x,0: ] = newrowx
        mat[ y,0: ] = newrowy

def findsmallestdiag( mat ):
        # which is the smallest (absolute) value on the diagonal of mat?
	diag = abs( matrix( mat ).diagonal() )
	minIndex = diag.argmin()

        return minIndex

def arr2lst( arr ):
        try:
		if arr.shape[0] == 1 :
			lst = arr.tolist()[0]
		else:
        		lst = arr.tolist()

        except AttributeError:
                lst = arr
        return lst
                  
def charpolmatrix( M, evalIndex = 0 ):
	# return a matrix A = M - l*I (l is the eigenvalue with index evalIndex)
	dim = M.shape[ 0 ]
	evals = eigvals( M )
	myeval = evals[ evalIndex ]
	return M - identity( dim )*myeval 

def decompose( matrix ):
	# Returns the decomposition of a matrix A where
	#
	# Q.A.Q = P.L.U
	#
	# P.L.U is the factoring of Q.A.Q such that L is a lower triangular matrix with 1's
	# on the diagonal and U is an upper triangular matrix; P is the permutation (row-swapping
	# operations) required for this procedure. The permutation matrix Q is chosen such that 
	# the last element of U is its smallest diagnoal element. If A has a zero eigenvalue, 
	# then U's last element will be zero.
	
	dim = matrix.shape[ 0 ]

	# first decomposition
	( P, L, U ) = lu( matrix )
	
 	# detect the smallest element of U
	smallestIndex = findsmallestdiag( U )
	smallest = U[ smallestIndex, smallestIndex ]

	#show( matrix, "M" )
	#show( U, "U" )
	#print "Smallest element is %f at %d" % ( smallest, smallestIndex )

	# is the permutation Q not just the identity matrix?
	Q = identity( dim )
	if smallestIndex+1 != dim :
		# trick: exchange row 'smallestIndex' with row 'dim-1' of the identity matrix
		swaprow( Q, smallestIndex, dim-1 )

	return ( P, L, U, Q )

class DirichletDistribution( object ):
	# Dirichlet distribution
	# alphas are a list of reals > 0 (counts)
	#
	# methods:
	# 	sample() - return a random sample from the distribution

	def __init__( self, counts ):
		
		self.alphas = matrix(counts) # + 1 # Sergio: commented out the +1 in order to make
						   # the input of ssaCalculator (priorCounts), the alpha 
						   # parameters instead of "pseudo counts".
		
		self.nparams = self.alphas.shape[ 1 ]
		self.alphasSum = self.alphas.sum()

	def sample( self ):
		# sample vector X = (x1, x2, ..., xK )
		# generated from xi = yi/Ytot
		# with each yi sampled from gamma distribution
		# 	p( yi ) = gammavariate( alphai, 1 )
		# and Ytot the sum of the yi's.

		ylist = []
		ysum = 0
		n = 0
		while n < self.nparams :
			alpha = self.alphas[ 0, n ]
			yi = gammavariate( alpha, 1 )
			ylist.append( yi )
			ysum += yi 
			n += 1
		
		xlist = []
		n = 0
		while n < self.nparams :
			xlist.append( ylist[ n ]/ysum )
			n += 1

		return tuple( xlist )

	def mean( self, index ):
		# return the expectation of parameter i: ai/aSum
		return float( self.alphas[ 0,index ] ) / float( self.alphasSum )

	def var( self, index ):
		# return the variance of parameter i:  
		ai = float( self.alphas[ index ] )
		a0 = float( self.alphasSum )
		num = ai*(a0 - ai )
		den = a0**2 * ( a0+1 )
		return num/den

