'''
QAP Analysis
Example data: http://anjos.mgi.polymtl.ca/qaplib/inst.html

Time to local optima: timetolocal()

Random problem generator: randomgenerator()


'''
import os
import numpy as np
import scipy as sp
import itertools
import math
import matplotlib as mpl
import glob
import pandas as pd
import time
import re

class QAP():
    '''
    QAP Object
    Example:
        q = qap.QAP('chr12a')
    '''
    def __init__(self, filename = None, silent = False):
        filename = re.sub('\.dat$', '', filename)
        self.filename = filename
        if self.filename is not None:
            self.outputs = {}
            self.read()
            self.dominance()
            try:
                self.sln()
            except:
                self.globalsol = None
                self.globalfit = None
                self.bestfit = 1.e100

            try:
                self.measure()
            except:
                self.variance = None
                self.rugdness = None


    def read(self):
        filename = self.filename + '.dat'
        if not os.path.isfile(filename):
            path = os.path.dirname(__file__)
            _filename = os.path.join(path, filename)
            _filename = os.path.abspath(_filename)
            if os.path.isfile(_filename):
                filename = _filename
            else:
                raise IOError('file {} not found'.format(filename))

        with open(filename, 'rt') as f:
            # self.logger_file_info(f)
            for i, line in enumerate(f):
                if i == 0:
                    #Matrix size
                    self.size = int(line.split()[0])
                    self.weight = np.ndarray(
                        (self.size, self.size),
                        dtype = 'int',
                    )
                    self.distance = np.ndarray(
                        (self.size, self.size),
                        dtype = 'int',
                    )
                if i == 1:
                    if len(line) > 1:
                        startweight = i
                    else:
                        startweight = i + 1
                if i >= 1:
                    if (startweight <= i < self.size + startweight):
                        ind = i - startweight
                        self.weight[ind,:] = np.fromstring(
                            line.strip(),
                            dtype = int,
                            sep = " ",
                        )
                    if i >= self.size + startweight:
                        if i == self.size + startweight:
                            if len(line) > 1:
                                startdistance = i
                            else:
                                startdistance = i + 1
                        elif (startdistance <= i < self.size + startdistance):
                            ind = i - startdistance
                            self.distance[ind,:] = np.fromstring(
                                line.strip(),
                                dtype = int,
                                sep = " ",
                            )

    def sln(self):
        '''Get the solution data if present'''
        filename = self.filename + '.sln'
        if not os.path.isfile(filename):
            path = os.path.dirname(__file__)
            _filename = os.path.join(path, "../qap", filename)
            _filename = os.path.abspath(_filename)
            if os.path.isfile(_filename):
                filename = _filename
            else:
                raise IOError('file {} not found'.format(filename))

        with open(filename, 'rt') as f:
            # self.logger_file_info(f)
            for i, line in enumerate(f):
                if i == 0:
                    #Matrix size
                    self.globalfit = float(line.split()[0])
                elif i == 1:
                    ind = i - 2
                    self.globalsol = np.fromstring(
                        line.strip(),
                        dtype = int,
                        sep = " ",
                    )

    def measure(self):
        '''Read the measure file'''
        filename = self.filename + '.measure'
        if not os.path.isfile(filename):
            path = os.path.dirname(__file__)
            _filename = os.path.join(path, "../qap", filename)
            _filename = os.path.abspath(_filename)
            if os.path.isfile(_filename):
                filename = _filename
            else:
                raise IOError('file {} not found'.format(filename))

        with open(filename, 'rt') as f:
            # self.logger_file_info(f)
            for i, line in enumerate(f):
                if i == 3:
                    self.variance = float(line.split()[0])
                    self.rugdness = float(line.split()[1])
                elif i == 5:
                    self.globalfit = int(line.split()[-1])


    def cost(self,sol):
        '''
        Solution is in the form [a, b, c, d, ...]
        where facility a is at location 1, b is at location 2...

        For the evaluation of cost, use a b c as the indices of the
        non-zero values of a NxN matrix. For example, the solution

        [1, 0, 2]

        would become

        s = [
            [0, 1, 0],
            [1, 0, 0],
            [0, 0, 1],
        ]

        d is the matrix that gives the distance between any two
        locations i and j.

        ==> s.d.sT gives the distance between any two
        facilities for a particular solution s where sT is s transposed

        ==> s.d.sT o w gives the cost associated with each
        link between facilities, where o is the entrywise product

        ==> The sum of all elements in s.d.sT o w is the total cost of
        the solution s
        '''

        s = np.zeros((self.size, self.size))
        
        for i in range(self.size):
            s[i,int(sol[i])] = 1
        return np.sum(np.dot(np.dot(s,self.distance),s.T) * self.weight)

    def dominance(self):
        matrices = [self.distance, self.weight]

        dom = []
        for m in matrices:
            dom += [100 * np.var(m) / np.mean(m)]

        self.dom_d = dom[0]
        self.dom_w = dom[1]
        self.dom_dw = (np.min(dom), np.max(dom))

        self.outputs['dom_d'] = (
            self.dom_d,
            'Distance dominance',
        )
        self.outputs['dom_w'] = (
            self.dom_d,
            'Weight dominance',
        )
        # self.outputs['dom_dw'] = (
        #     self.dom_dw,
        #     'Distance/Weight dominance',
        # )


    def initialize(self, n):
        '''Initialize random population'''
        sols = np.ndarray((n,self.size), dtype = 'int')
        for i in range(n):
            sols[i] = np.random.permutation(self.size)

        return sols

    def initall(self):
        '''Initialize every possible solution'''
        sols = np.ndarray((sp.misc.factorial(self.size),self.size), dtype = 'int')
        for i, row in enumerate(itertools.permutations(range(self.size))):
            sols[i,:] = row

        return sols

    def neighbors(self, sol):
        '''Find 2-opt neighbors of a solution'''
        neighbors = np.ndarray((int(0.5*self.size*(self.size - 1)), self.size),dtype = 'int')
        neighbors[:,:] = sol
    
        row = 0
        for i in range(self.size):
            for j in range(i+1, self.size):
                neighbors[row,i] = sol[j]
                neighbors[row,j] = sol[i]
                row += 1

        return neighbors

    def randomwalk(
        self,
        steps = 10,
    ):
        sol = np.ndarray((steps, self.size))
        sol[0] = self.initialize(1)

        for i in range(1, steps):
            sol[i] = sol[i - 1]
            swaps = np.random.randint(self.size, size = 2)
            sol[i, swaps[0]] = sol[i-1, swaps[1]]
            sol[i, swaps[1]] = sol[i-1, swaps[0]]
        return sol

    def autocorrlen(
        self,
        n = 100000,
        walk = None,
    ):
        if walk is None:
            sols = self.randomwalk(n)
        else:
            sols = walk
        solcosts=[]
        for sol in sols:
            solcosts.append(self.cost(sol))
        costs = pd.Series(solcosts)
        walks = pd.DataFrame({'original': costs[1:-1], 'shifted': costs.shift()[1:-1]})
        rho1 = walks.corr()['original']['shifted']
        length = -1 / np.log(np.abs(rho1))

        return length

    def coefficientVariation(
        self,
        steps = 1000,
        nw=100
        ):
        lw = int(steps/nw)
        sol = self.randomwalk(steps)
        costs = pd.Series(self.cost(sol))
#        i = 0
#
#        while i<steps-lw:
#            mean= costs[i:i+lw].mean()
#            var= costs[i:i+lw].var()
#            i+=lw;
#            print(var/mean)

        cv = np.ndarray((nw))

        for i in range(nw):
            start = i * lw
            end = start + lw
            mean = costs[start:end].mean()
            var = costs[start:end].var()
            cv[i] = var/mean
        cv = pd.Series(cv)
        ax = cv.plot()
        fig = ax.get_figure()
        fig.savefig(self.filename+'.cv.plot.pdf')
        return cv

    def infostability(
        self,
        steps = 100,
        walk = None,
        guess = 100,
    ):
        '''Calculate the critical epsilon for information stability'''
        epsilon = guess
        S, H, M, walk, costs = self.infocontent(steps = steps, walk = walk, epsilon = epsilon, returnwalk = True)
        diff = costs.diff()
        #Find an approximation
        for i in range(30):
            stable = self.stabilitycheck(diff, epsilon)
            if not stable:
                epsilonold = epsilon
                epsilon *= 2
            else:
                break
        #Bisection search
        upper = epsilon
        lower = epsilonold
        for i in range(20):
            mid = 0.5 * (upper + lower)
            stable = self.stabilitycheck(diff, mid)
            if stable:
                upper = mid
            else:
                lower = mid
        if stable:
            epsilon = upper
        else:
            epsilon = 2 * upper - lower
        return epsilon

    def stabilitycheck(
        self,
        diff,
        epsilon,
    ):
        combed = diff.copy(True)
        combed[abs(diff) < epsilon] = 0
        S = np.sign(combed).shift(-1)[:-1].astype(int)
        stable = (S == 0).all(0)
        return stable


    def walktest(
        self,
        steps = 100000,
    ):
        walk = self.randomwalk(steps)
        acl = self.autocorrlen(walk = walk)

        self.outputs['autocorrlen'] = (acl,'Autocorrelation length')

    def localsearch(
        self,
        n = 2000,
        plot = False,
        ):
        start = time.time()
        #Initialize n random solutions
        sols = self.initialize(n)
        #Initialize paths as a dictionary
        paths = {}
        #Array for locally optimal solutions
        self.localsols = np.empty_like(sols)
        for i in range(n):
            #Initialize cost path as an array
            path = []
            #Set the old cost so that it never converges initially
            oldcost = 1.0e100
            #Grab the current solution
            sol = sols[i]
            #Calculate the intial cost
            cost = self.cost(sol)
            #Write the initial cost to the path
            path += [cost]
            #Record all of the solutions that are at the current best fitness
            sublist = np.ndarray((0, self.size))
            x=0
            while (cost < oldcost):
                #Find neighbors
                neighbors = self.neighbors(sol)
                bestNeighCost = self.cost(neighbors[0])
                ind = 0
                bestind = 0
                for neigh in neighbors:
                    #Calculate cost of neighbors
                    neighborcost = self.cost(neigh)
                    #Find best neighbor
                    if(neighborcost<bestNeighCost):
                        bestNeighCost = neighborcost
                        bestind = ind
                    ind = ind+1
                oldcost = cost
                cost = bestNeighCost
                #Prevent infinite loops by checking if a previous solution of
                #given fitness has already been visited
                if cost < oldcost:
                    sublist = np.ndarray((0, self.size))
                else:
                    if sol.tolist() in sublist.tolist():
                        break
                    else:
                        sublist = np.concatenate((sublist, sol[np.newaxis,:]))

                path += [cost]
                sol = neighbors[bestind]
            #Write final solution into array of the locally optimal solutions
            self.localsols[i] = sol
            #Write each path taken as a Series
            paths[i] = pd.Series(np.array(path)[:-1])
        #Convert paths to a DataFrame
        paths = pd.DataFrame(paths)

        #Calculate path lengths by counting the length of paths
        pathlengths = paths.count(numeric_only = True)

        #Fitness of the local minima
        self.localmin = paths.min()

        #Find the index of the best local optima found in the localsearch
        bestind = self.localmin.idxmin()
        if self.globalsol is None:
            if (self.localmin[bestind] < self.bestfit):
                #If the global optima are not given, then record
                #the best solution found
                self.bestsol = self.localsols[bestind]
                self.bestfit = self.localmin[bestind]
            if self.globalfit is not None:
                #If the solution is not given, but the fitness is, then
                #if a solution of corresponding fitness is found, record this
                #as the global optimum
                if float(self.bestfit) == float(self.globalfit):
                    self.globalsol = self.bestsol



        #Find all paths that lead to the global optimum
        if self.globalfit is not None:
            ind = paths.min() == self.globalfit
        else:
            ind = paths.min() == self.bestfit
        gpaths = paths[ind[ind == True].index]
        gpathlengths = gpaths.count(numeric_only = True)

#print(time.time() - start)


        self.outputs['variance'] = (
            self.variance,
            'Variance',
        )
        self.outputs['rugdness'] = (
            self.rugdness,
            'Rugdness',
        )

    def locprop(self, n = None, plot = False):
        start = time.time()
        if n is not None:
            sols = self.initialize(n)
        else:
            sols = self.initall()
            n = sols.shape[0]
        isloc = np.ndarray((n), dtype = 'bool')
        isloc[:] = False
        costs = np.ndarray((n), dtype = 'int')
        for i in range(n):
            costs[i] = self.cost(sols[i, np.newaxis, :])[0]
            #Find neighbors
            neighbors = self.neighbors(sols[i,:])
            #Calculate cost of neighbors
            neighborcosts = self.cost(neighbors)
            #Determine if local minimum
            if np.min(neighborcosts) >= costs[i]:
                isloc[i] = True

        lower = np.min(costs) - 1
        # upper = np.max(costs) + 1
        upper = np.min(costs) + 0.25*(np.max(costs) - np.min(costs))
        bins = np.linspace(lower, upper, 20)

        loc = costs[isloc]
        nonloc = costs[np.invert(isloc)]

        labels = []
        for i in bins[:-1]:
            labels += [int(i)]
        loc = pd.cut(loc, bins, labels = labels)
        nonloc = pd.cut(nonloc, bins, labels = labels)

        loc = pd.Series(loc)
        nonloc = pd.Series(nonloc)

        loc = loc[loc != 0]
        nonloc = nonloc[nonloc != 0]

        count_loc = loc.value_counts()
        count_nonloc = nonloc.value_counts()
        # print(loc)

        proportion = count_loc.div(
            count_loc.add(
                count_nonloc,
                fill_value = 0,
            ),
            fill_value = 0,
        )

        # print(count_loc, count_nonloc)
        if plot:
            fig, axes = plt.subplots(
                nrows = 1,
                ncols = 1,
            )
            proportion.plot(kind = 'bar', ax = axes)
            axes.set_ylabel('cost')
            axes.set_xlabel('proportion')
