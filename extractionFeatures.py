import numpy as np
import itertools
import readXY

############################################################

m1 = list(itertools.product('ACGT', repeat=1))
m2 = list(itertools.product('ACGT', repeat=2))
m3 = list(itertools.product('ACGT', repeat=3))
m4 = list(itertools.product('ACGT', repeat=4))
m5 = list(itertools.product('ACGT', repeat=5))


def isValid(x):
    for base in x:
        if not base in {'A','C','G','T'}:
            return 'Invalid'
    return 'Valid'

def kmers(seq, k):
    v = []
    for i in range(len(seq) - k + 1):
        v.append(seq[i:i + k])
    return v

def pseudoKNC(x):
    ### k-mer ###
    ### A, AA, AAA
    for i in range(1, k + 1, 1):
        v = list(itertools.product('ACGT', repeat=i))
        # seqLength = len(x) - i + 1
        for i in v:
            # print(x.count(''.join(i)), end=',')
            t.append(x.count(''.join(i)))
    ### --- ###

def zCurve(x):
    ### Z-Curve ### total = 3
    A = x.count('A'); C = x.count('C'); G = x.count('G'); T = x.count('T')
    x_ = (A + G) - (C + T); y_ = (A + C) - (G + T); z_ = (A + T) - (C + G)
    # print(x_, end=','); print(y_, end=','); print(z_, end=',')
    t.append(x_); t.append(y_); t.append(z_)
    ### print('{},{},{}'.format(x_, y_, z_), end=',')
    ### --- ###

def gapping1_1(x, g):
    ### g-gap
    '''
    AA      0-gap (2-mer)
    A_A     1-gap
    A__A    2-gap
    A___A   3-gap
    A____A  4-gap
    '''

    m = m2

    for i in range(0, g + 1, 1):
        V = kmers(x, i + 2)
        # seqLength = len(x) - (i+2) + 1
        #
        for gGap in m:
            # print(gGap[0], end='')
            # print('-'*i, end='')
            # print(gGap[1], end=' ')

            C = 0
            for v in V:
                if v[0] == gGap[0] and v[-1] == gGap[1]:
                    C += 1
            # print(C, end=',')
            t.append(C)

    ### --- ###

def gapping1_2(x, g):

    m = m3

    for i in range(1, g + 1, 1):
        V = kmers(x, i + 3)
        # seqLength = len(x) - (i+2) + 1

        # print(V)
        for gGap in m:
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print('-'*i, end='')
            # print(gGap[2], end=' ')
            C = 0
            for v in V:
                if v[0] == gGap[0] and v[1] == gGap[1] and v[-1] == gGap[2]:
                    C += 1
            # print(C, end=',')
            t.append(C)

            # print(gGap[0], end='')
            # print('-' * i, end='')
            # print(gGap[1], end='')
            # print(gGap[2], end=' ')

            C = 0
            for v in V:
                if v[0] == gGap[0] and v[-2] == gGap[1] and v[-1] == gGap[2]:
                    C += 1
            # print(C, end=',')
            t.append(C)

def gapping1_3(x, g):

    m = m4

    for i in range(1, g + 1, 1):
        V = kmers(x, i + 4)
        # seqLength = len(x) - (i+2) + 1

        # print(V)
        for gGap in m:
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print('-'*i, end='')
            # print(gGap[2], end=' ')
            C = 0
            for v in V:
                if v[0] == gGap[0] and v[1] == gGap[1] and v[2] == gGap[2] and v[-1] == gGap[3]:
                    C += 1
            # print(C, end=',')
            t.append(C)
            # print(gGap[0], end='')
            # print('-' * i, end='')
            # print(gGap[1], end='')
            # print(gGap[2], end=' ')

            C = 0
            for v in V:
                if v[0] == gGap[0] and v[-3] == gGap[1] and v[-2] == gGap[2] and v[-1] == gGap[3]:
                    C += 1
            # print(C, end=',')
            t.append(C)

def gapping2_2(x, g):
    ### gapping ### total = [(64xg)] = 2,304 [g=9]
    # AA_AA       1-gap
    # AA__AA      2-gap
    # AA___AA     3-gap
    # AA____AA    4-gap
    # AA_____AA   5-gap upto g

    m = m4

    for i in range(1, g + 1, 1):
        V = kmers(x, i + 4)
        # seqLength = len(x) - (i+2) + 1
        # print(V)
        for gGap in m:
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print(gGap[2], end='')
            # print('-'*i, end=' ')
            C = 0
            for v in V:
                if v[0] == gGap[0] and v[1] == gGap[1] and v[-2] == gGap[2] and v[-1] == gGap[3]:
                    C += 1
            print(C, end=',')

    ### --- ###

def gapping2_3(x, g):
    ### gapping ### total = [(64xg)] = 2,304 [g=9]
    # AAA_AA       1-gap
    # AAA__AA      2-gap
    # AAA___AA     3-gap
    # AAA____AA    4-gap
    # AAA_____AA   5-gap upto g

    # AA_AAA       1-gap
    # AA__AAA      2-gap
    # AA___AAA     3-gap
    # AA____AAA    4-gap
    # AA_____AAA   5-gap upto g

    m = m5

    for i in range(1, g + 1, 1):
        V = kmers(x, i + 5)
        # seqLength = len(x) - (i+2) + 1
        # print(V)
        for gGap in m:
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print(gGap[2], end='')
            # print('-'*i, end=' ')
            C = 0
            for v in V:
                if v[0] == gGap[0] and v[1] == gGap[1] and v[2] == gGap[2] and v[-2] == gGap[3] and v[-1] == gGap[4]:
                    C += 1
            # print(C, end=',')
            t.append(C)
            # print('-' * i, end='')
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print(gGap[2], end=' ')
            C = 0
            for v in V:
                if v[0] == gGap[0] and v[1] == gGap[1] and v[-3] == gGap[2] and v[-2] == gGap[3] and v[-1] == gGap[4]:
                    C += 1
            # print(C, end=',')
            t.append(C)

def gapping1_(x, g):
    ### gapping ### total = [(16xg)x2] = 288 [g=9]
    # A_       1-gap
    # A__      2-gap
    # A___     3-gap
    # A____    4-gap
    # A_____   5-gap upto g

    # _A       1-gap
    # __A      2-gap
    # ___A     3-gap
    # ____A    4-gap
    # _____A   5-gap upto g

    m = m1

    for i in range(1, g + 1, 1):
        V = kmers(x, i + 1)
        # seqLength = len(x) - (i+2) + 1
        # print(V)
        for gGap in m:
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print('-'*i, end='')
            C = 0
            for v in V:
                if v[0] == gGap[0]:
                    C += 1
            print(C, end=',')
            # print('-' * i, end='')
            # print(gGap[0], end='')
            # print(gGap[1], end=' ')
            C = 0
            for v in V:
                if v[-1] == gGap[0]:
                    C += 1
            print(C, end=',')

def gapping2_(x, g):
    ### gapping ### total = [(16xg)x2] = 288 [g=9]
    # AA_       1-gap
    # AA__      2-gap
    # AA___     3-gap
    # AA____    4-gap
    # AA_____   5-gap upto g

    # _AA       1-gap
    # __AA      2-gap
    # ___AA     3-gap
    # ____AA    4-gap
    # _____AA   5-gap upto g

    m = m2
    for i in range(1, g + 1, 1):
        V = kmers(x, i + 2)
        # seqLength = len(x) - (i+2) + 1
        # print(V)
        for gGap in m:
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print('-'*i, end='')
            C = 0
            for v in V:
                if v[0] == gGap[0] and v[1] == gGap[1]:
                    C += 1
            print(C, end=',')
            # print('-' * i, end='')
            # print(gGap[0], end='')
            # print(gGap[1], end=' ')
            C = 0
            for v in V:
                if v[-2] == gGap[0] and v[-1] == gGap[1]:
                    C += 1
            print(C, end=',')

def gapping3_(x, g):
    ### gapping ### total = [(64xg)x2] = 1,152 [g=9]
    # AAA_       1-gap
    # AAA__      2-gap
    # AAA___     3-gap
    # AAA____    4-gap
    # AAA_____   5-gap upto g

    # _AAA       1-gap
    # __AAA      2-gap
    # ___AAA     3-gap
    # ____AAA    4-gap
    # _____AAA   5-gap upto g

    m = list(itertools.product('ACGT', repeat=3))
    for i in range(1, g + 1, 1):
        V = kmers(x, i + 3)
        # seqLength = len(x) - (i+2) + 1
        # print(V)
        for gGap in m:
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print(gGap[2], end='')
            # print('-'*i, end=' ')
            C = 0
            for v in V:
                if v[0] == gGap[0] and v[1] == gGap[1] and v[2] == gGap[2]:
                    C += 1
            print(C, end=',')
            # print('-' * i, end='')
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print(gGap[2], end=' ')
            C = 0
            for v in V:
                if v[-3] == gGap[0] and v[-2] == gGap[1] and v[-1] == gGap[2]:
                    C += 1
            print(C, end=',')

def gapping_1_(x, g):
    ### g-gap ### total = (16x(g+1)) = 160 [g=9]
    # _A_       1-gap
    # _A__      2-gap
    # _A___     3-gap
    # _A____    4-gap
    # _A_____   5-gap upto g

    # _A_       1-gap
    # __A_      2-gap
    # ___A_     3-gap
    # ____A_    4-gap
    # _____A_   5-gap upto g

    m = list(itertools.product('ACGT', repeat=1))
    for i in range(1, g + 1, 1):
        V = kmers(x, i + 2)
        # seqLength = len(x) - (i+2) + 1
        for gGap in m:
            # print('-', end='')
            # print(gGap[0], end='')
            # print('-'*i, end='')

            C = 0
            for v in V:
                if v[1] == gGap[0]:
                    C += 1
            print(C, end=',')

            if i > 1:
                # print('-'*i, end='')
                # print(gGap[0], end='')
                # print('-', end='')
                C = 0
                for v in V:
                    if v[-2] == gGap[0]:
                        C += 1
                print(C, end=',')

    ### --- ###

def gapping_2_(x, g):
    ### g-gap ### total = (16x(g+1)) = 160 [g=9]
    # _AA_       1-gap
    # _AA__      2-gap
    # _AA___     3-gap
    # _AA____    4-gap
    # _AA_____   5-gap upto g

    # _AA_       1-gap
    # __AA_      2-gap
    # ___AA_     3-gap
    # ____AA_    4-gap
    # _____AA_   5-gap upto g

    m = list(itertools.product('ACGT', repeat=2))
    for i in range(1, g + 1, 1):
        V = kmers(x, i + 2)
        # seqLength = len(x) - (i+2) + 1
        #
        for gGap in m:
            # print('-', end='')
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print('-'*i, end='')

            C = 0
            for v in V:
                if v[1] == gGap[0] and v[2] == gGap[1]:
                    C += 1
            print(C, end=',')

            if i > 1:
                # print('-'*i, end='')
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print('-', end='')

                C = 0
                for v in V:
                    if v[-3] == gGap[0] and v[-2] == gGap[1]:
                        C += 1
                print(C, end=',')

def gapping_3_(x, g):
    ### g-gap ### total = (16x(g+1)) = 160 [g=9]
    # _AAA_       1-gap
    # _AAA__      2-gap
    # _AAA___     3-gap
    # _AAA____    4-gap
    # _AAA_____   5-gap upto g

    # _AAA_       1-gap
    # __AAA_      2-gap
    # ___AAA_     3-gap
    # ____AAA_    4-gap
    # _____AAA_   5-gap upto g

    m = list(itertools.product('ACGT', repeat=3))
    for i in range(1, g + 1, 1):
        V = kmers(x, i + 3)
        # seqLength = len(x) - (i+2) + 1
        #
        for gGap in m:
            # print('-', end='')
            # print(gGap[0], end='')
            # print(gGap[1], end='')
            # print(gGap[2], end='')
            # print('-'*i, end='')

            C = 0
            for v in V:
                if v[1] == gGap[0] and v[2] == gGap[1] and v[3] == gGap[2]:
                    C += 1
            print(C, end=',')

            if i > 1:
                # print('-'*i, end='')
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print('-', end='')

                C = 0
                for v in V:
                    if v[-4] == gGap[0] and v[-3] == gGap[1] and v[-2] == gGap[2]:
                        C += 1
                print(C, end=',')

def none():
    # 5-l with 1-dance (256*5 = 1280 features)
    '''
    _AAAA, A_AAA, AA_AA AAA_A, AAAA_
    '''
    l = -5
    m = list(itertools.product('ACGT', repeat=4))
    for i in range(1, (l + 1), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('-', end='')
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print(gGap[3], end=' ')
                C = 0
                for v in V:
                    if v[1] == gGap[0] and v[2] == gGap[1] and v[3] == gGap[2] and v[4] == gGap[3]:
                        C += 1
                print(C, end=',')

        elif i == 2:
            for gGap in m:
                # print(gGap[0], end='')
                # print('-', end='')
                # print(gGap[1], end='')
                # print(gGap[2], end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[2] == gGap[1] and v[3] == gGap[2] and v[4] == gGap[3]:
                        C += 1
                print(C, end=',')

        elif i == 3:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print('-', end='')
                # print(gGap[2], end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1] and v[3] == gGap[2] and v[4] == gGap[3]:
                        C += 1
                print(C, end=',')

        elif i == 4:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print('-', end='')
                # print(gGap[2], end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1] and v[2] == gGap[2] and v[4] == gGap[3]:
                        C += 1
                print(C, end=',')



        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print('-', end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1] and v[2] == gGap[2] and v[3] == gGap[3]:
                        C += 1
                print(C, end=',')

    # 5-l with 2-dance (64*4 = 256 features)
    '''
    __AAA, A__AA, AA__A AAA__
    '''
    l = -5
    m = list(itertools.product('ACGT', repeat=3))
    for i in range(1, (l), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('-', end='')
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print(gGap[3], end=' ')
                C = 0
                for v in V:
                    if v[2] == gGap[0] and v[3] == gGap[1] and v[4] == gGap[2]:
                        C += 1
                print(C, end=',')

        elif i == 2:
            for gGap in m:
                # print(gGap[0], end='')
                # print('-', end='')
                # print(gGap[1], end='')
                # print(gGap[2], end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[3] == gGap[1] and v[4] == gGap[2]:
                        C += 1
                print(C, end=',')

        elif i == 3:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print('-', end='')
                # print(gGap[2], end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1] and v[4] == gGap[2]:
                        C += 1
                print(C, end=',')


        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print('-', end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1] and v[2] == gGap[2]:
                        C += 1
                print(C, end=',')

    # 4-l with 1-dance (64*4 = 256 features)
    '''
    _AAA, A_AA, AA_A, AAA_
    '''
    l = -4
    m = list(itertools.product('ACGT', repeat=3))
    for i in range(1, (l + 1), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('-' * i, end='')
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end=' ')
                C = 0
                for v in V:
                    if v[1] == gGap[0] and v[2] == gGap[1] and v[3] == gGap[2]:
                        C += 1
                print(C, end=',')
        elif i == 2:
            for gGap in m:
                # print(gGap[0], end='')
                # print('-', end='')
                # print(gGap[1], end='')
                # print(gGap[2], end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[2] == gGap[1] and v[3] == gGap[2]:
                        C += 1
                print(C, end=',')

        elif i == 3:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print('-', end='')
                # print(gGap[2], end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1] and v[3] == gGap[2]:
                        C += 1
                print(C, end=',')


        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print('-', end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1] and v[2] == gGap[2]:
                        C += 1
                print(C, end=',')

    # 3-l with 1-dance (16*2 = 32 features)
    '''
    _AA, AA_,
    '''
    l = -3
    m = list(itertools.product('ACGT', repeat=2))
    for i in range(1, (l), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('-', end='')
                # print(gGap[0], end='')
                # print(gGap[1], end=' ')
                C = 0
                for v in V:
                    if v[1] == gGap[0] and v[2] == gGap[1]:
                        C += 1
                print(C, end=',')


        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print('-', end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1]:
                        C += 1
                print(C, end=',')

def moving5_34(x):
    # 5-l with 3-dance (16*2 = 32 features)
    '''
    ___AA, AA___
    '''
    l = 5
    m = m2
    for i in range(1, (l - 2), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('-', end='')
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print(gGap[3], end=' ')
                C = 0
                for v in V:
                    if v[3] == gGap[0] and v[4] == gGap[1]:
                        C += 1
                # print(C, end=',')
                t.append(C)

        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print('-', end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1]:
                        C += 1
                # print(C, end=',')
                t.append(C)

    # 5-l with 4-dance (4*2 = 8 features)
    '''
    ____A, A____
    '''
    l = 5
    m = list(itertools.product('ACGT', repeat=1))
    for i in range(1, (l - 2), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('-', end='')
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print(gGap[3], end=' ')
                C = 0
                for v in V:
                    if v[4] == gGap[0]:
                        C += 1
                # print(C, end=',')
                t.append(C)

        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print(gGap[2], end='')
                # print('-', end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0]:
                        C += 1
                # print(C, end=',')
                t.append(C)

def moving4_23(x):
    # 4-l with 2-dance (16*2 = 32 features)
    '''
    __AA, AA__
    '''
    l = 4
    m = m2
    for i in range(1, (l - 1), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('--', end='')
                # print(gGap[0], end='')
                # print(gGap[1], end=' ')
                C = 0
                for v in V:
                    if v[2] == gGap[0] and v[3] == gGap[1]:
                        C += 1
                # print(C, end=',')
                t.append(C)

        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print(gGap[1], end='')
                # print('--', end=' ')

                C = 0
                for v in V:
                    if v[0] == gGap[0] and v[1] == gGap[1]:
                        C += 1
                # print(C, end=',')
                t.append(C)


    # 4-l with 3-dance (4*2 = 8 features)
    '''
    ___A, A___
    '''
    l = 4
    m = list(itertools.product('ACGT', repeat=1))
    for i in range(1, (l - 1), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('___', end='')
                # print(gGap[0], end=' ')
                C = 0
                for v in V:
                    if v[3] == gGap[0]:
                        C += 1
                # print(C, end=',')
                t.append(C)

        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print('___', end=' ')
                C = 0

                for v in V:
                    if v[1] == gGap[0]:
                        C += 1
                # print(C, end=',')
                t.append(C)

def moving3_2(x):
    # 3-l with 2-dance (4*2 = 8 features)
    '''
    __A, A__,
    '''
    l = 3
    m = m1
    for i in range(1, (l), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('--', end='')
                # print(gGap[0], end=' ')
                C = 0
                for v in V:
                    if v[2] == gGap[0]:
                        C += 1
                # print(C, end=',')
                t.append(C)

        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print('--', end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0]:
                        C += 1
                # print(C, end=',')
                t.append(C)

def moving2_1(x):
    # 2-l with 1-dance (4*2 = 8 features)
    '''
    _A, A_,
    '''
    l = 2
    m = m1
    for i in range(1, (l + 1), 1):
        V = kmers(x, l)
        # print(v)
        # seqLength = len(x) - (i+2) + 1

        # print(V)

        if i == 1:
            for gGap in m:
                # print('-', end='')
                # print(gGap[0], end=' ')
                C = 0
                for v in V:
                    if v[1] == gGap[0]:
                        C += 1
                # print(C, end=',')
                t.append(C)

        else:
            for gGap in m:
                # print(gGap[0], end='')
                # print('-', end=' ')
                C = 0
                for v in V:
                    if v[0] == gGap[0]:
                        C += 1
                # print(C, end=',')
                t.append(C)

def generateFeatures(k, g, x, y):
    ## begin okay ##
    gapping1_1(x, g)
    gapping1_2(x, g)
    gapping1_3(x, g)

    moving5_34(x)
    moving4_23(x)
    moving3_2(x)
    moving2_1(x)
    ## end okay ##

    ## begin newFeatures ##
    # pseudoKNC()
    zCurve(x)
    # gapping1_(x, g)
    # gapping2_(x, g)
    # gapping3_(x, g)

    # gapping_1_(x, g)
    # gapping_2_(x, g)
    # gapping_3_(x, g)

    # gapping2_2(x, g)
    gapping2_3(x, g)

    ## end newFeatures ##

    ######## class ########
    # print(y)
    t.append(y)
    #######################

    # T = [] # All instance ...

def saveCSV(X, Y):

    F = open('fullDataset.csv', 'w')

    for x, y in zip(X, Y):
        for each in x:
            F.write(str(each) + ',')
        F.write(str(y) + '\n')

    F.close()


if __name__ == '__main__':

    k = 0   # if k=0 then k-mer won't work
    g = 15  # if g=-1 then g-gap won't work # [0 to 8] are the best

    X, Y = readXY.fetchXY()

    T = []
    for x, y in zip(X, Y):
        t = []
        generateFeatures(k, g, x, y)
        t = np.array(t)
        T.append(t)

    T = np.array(T)

    X = T[:,:-1]
    Y = T[:,-1]

    saveCSV(X,Y)





