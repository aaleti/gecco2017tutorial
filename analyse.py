import pandas as pd
from qap import QAP
import numpy as np
class Analyze():
    '''
    Compare the properties of multiple problems
    Parameters:
        path: directory containing problem files
        string: string matching the filenames
        n: number of randomly started local searches to perform
        steps: number of steps taken on the random walk

    Example:
        To compare all QAP problems with filenames beginning with 'chr':
        a = qap.Analyze(
            path = '/path/to/problems',
            string = 'chr*.dat',
            n = 100,
            steps = 1000,
        )

        a.data is a Pandas Dataframe containing the properties of each
        problem run. The usual Pandas tricks can be applied, such as

        a.data.mean() will find the means
        a.data.corr() will find the correlations
        a.data.to_csv() will convert the data to CSV format

        MPL plots can be made, for example:
            a.plot('variance', 'fdc')
        will plot fdc against variance

    '''
    def __init__(self,
        problem = QAP,
        datfile = 'overlap/ov0_adqap_16.0.0.5.0.5.2.dat',
        n = 5000, #Local search
        steps = 500000, #Random walk
    ):
        # properties = [
        #     'overlap',
        #     'rugdness',
        #     'variance',
        #     'meanpathlength',
        #     'pathlengthratioglobal',
        #     'pathlengthratioall',
        #     'flipdistlocal',
        #     'flipdistglobal',
        #     'locflipdist',
        #     'fdc',
        #     'infocontent',
        #     'partialinfocontent',
        #     'infostability',
        #     'autocorrlen',
        # ]
        #Create a test object to know all the attributes
        testobj = problem(filename = datfile)
        testobj.localsearch(n = 10)
        testobj.walktest(steps = 10)
        properties = testobj.outputs.keys()

        self.data = pd.DataFrame()
        self.labels = {}
        for p in properties:
            self.data[p] = np.ndarray(1)


        q = problem(filename = datfile)
        q.localsearch(n = n)
        q.walktest(steps = steps)
            # S, H, M = q.infocontent(steps = n)
            # epsilon = q.infostability()
            # acl = q.autocorrlen()
            #Get data
        for p in properties:
                # try:
            self.data[p][i] = q.outputs[p][0]
                # except:
                #     if p == 'infocontent':
                #         self.data[p][i] = H
                #     elif p == 'partialinfocontent':
                #         self.data[p][i] = M
                #     elif p == 'infostability':
                #         self.data[p][i] = epsilon
                #     elif p == 'autocorrlen':
                #         self.data[p][i] = acl
        #Get labels
        for p in properties:
            # try:
            self.labels[p] = q.outputs[p][1]
            # except:
            #     if p == 'infocontent':
            #         self.labels[p] = 'Information content'
            #     elif p == 'partialinfocontent':
            #         self.labels[p] = 'Partial information content'
            #     elif p == 'infostability':
            #         self.labels[p] = 'Information stability'
            #     elif p == 'autocorrlen':
            #         self.labels[p] = 'Autocorrelation length'

#print(time2human(time.time() - start))

    def plot(self, x, y):
        '''
        Plot x against y
        '''
        ax = self.data.plot(
            x = x,
            y = y,
            linewidth = 0,
            marker = 'x',
            legend = False,
        )
        ax.set_xlabel(self.labels[x])
        ax.set_ylabel(self.labels[y])
        ax.text(
            0.95,
            0.95,
            '{:4.3f}'.format(self.data.corr()[x][y]),
            horizontalalignment='right',
            verticalalignment='top',
            transform = ax.transAxes,
        )
        ax.set_xlim(self.data[x].min(),self.data[x].max())
        ax.set_ylim(self.data[y].min(),self.data[y].max())
