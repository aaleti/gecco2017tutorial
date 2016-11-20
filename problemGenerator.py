
def randomgenerator(size = 10, number = 3, prefix = 'cc'):
    for i in range(number):
        distance = np.random.randint(100, size = (size, size))
        distance[np.tril_indices(size)] = 0
        distance = distance + distance.T

        weight = np.random.randint(100, size = (size, size))
        weight[np.tril_indices(size)] = 0
        weight = distance + distance.T

        path = os.path.dirname(__file__)
        filename = '{}_{}_{}.dat'.format(prefix, size, i)
        _filename = os.path.join(path, "../qap", filename)
        _filename = os.path.abspath(_filename)

        with open(_filename, 'wt') as f:
            f.write('{}\n'.format(size))
            f.write('\n')

            for i in range(size):
                f.write(' '.join([str(n) for n in distance[i,:]]) + '\n')
            f.write('\n')

            for i in range(size):
                f.write(' '.join([str(n) for n in weight[i,:]]) + '\n')
