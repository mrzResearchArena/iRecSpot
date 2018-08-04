nameFASTA = []

def readFASTA(fileName):
    with open(fileName, 'r') as file:
        v = []
        genome = ''
        for line in file:
            if line[0] != '>':
                genome += line.strip()
            else:
                nameFASTA.append(line.replace('\n', ''))
                v.append(genome)
                genome = ''
        v.append(genome)
        del v[0]
        return v


def fetchX():
    X = readFASTA('testFASTA.fasta')
    return X


def isValid(x):
    for base in x:
        if not base in {'A','C','G','T'}:
            return 'Invalid'
    return 'Valid'


##################################################
