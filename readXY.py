
def readFASTA(fileName):
    with open(fileName, 'r') as file:
        v = []
        genome = ''
        for line in file:
            if line[0] != '>':
                genome += line.strip()
            else:
                v.append(genome)
                genome = ''
        v.append(genome)
        del v[0]
        return v


def fetchXY():
    X_hot = readFASTA('hotSpot.fasta')
    X_cold = readFASTA('coldSpot.fasta')
    X = X_hot + X_cold

    Y_hot = [+1 for _ in range(1, 478 + 1, 1)]  # 478
    Y_cold = [-1 for _ in range(1, 572 + 1, 1)]  # 572
    Y = Y_hot + Y_cold

    return X, Y




