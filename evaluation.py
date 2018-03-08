import numpy as np
import itertools
import readiRecX

############################################################
picklePath = 'iRecSpotModel.pkl'

############################################################

############################################################
X = readiRecX.fetchX()
Y = [-100 for _ in range(len(X))] # assign garbage values

print('Please enter the stepSize:')
stepSize = int(input())

############################################################


############################################################
def FASTAtoPrediction(stepSize, X, Y):

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

    def evaluateModel(X_test):

        storeMeanSD = np.array([
      [235.52857143, 235.84764582],
      [1.01142857, 1.23282983],
      [0.97714286, 1.18057751],
      [1.37428571, 1.75116443],
      [1.02952381, 1.27319798],
      [1.21333333, 1.48907556],
      [1.05333333, 1.36312433],
      [2.04, 1.98382027],
      [11.44380952, 8.45656605],
      [2.28666667, 2.13379013],
      [1.97142857, 2.31786907],
      [9.33142857, 7.16929052],
      [0.34571429, 0.64550144],
      [1.05619048, 1.26930231],
      [1.47333333, 1.78559136],
      [0.99714286, 1.2440301],
      [0.48285714, 0.75981917],
      [0.60857143, 1.00480885],
      [9.10857143, 8.92367759],
      [1.00761905, 1.20630844],
      [6.30571429, 5.52353111],
      [0.81809524, 1.08419652],
      [3.1047619, 2.94140818],
      [1.52666667, 1.70705911],
      [4.33428571, 3.61157738],
      [6.76285714, 5.88346767],
      [4.9752381, 4.33318932],
      [1.53047619, 1.73466746],
      [2.58857143, 2.85694856],
      [0.42666667, 0.74248831],
      [3.18571429, 3.72261886],
      [5.6, 5.14050211],
      [4.19619048, 4.17584719],
      [1.26666667, 1.44595769],
      [0.6952381, 0.95442637],
      [1.09809524, 1.34479462],
      [1.38, 1.61226844],
      [2.40095238, 2.44290987],
      [0.68285714, 0.9697282],
      [1.83904762, 2.02503358],
      [1.02285714, 1.23421517],
      [13.62857143, 9.93690298],
      [1.43428571, 1.70600221],
      [1.20952381, 1.40200698],
      [18.10761905, 12.74597375],
      [0.68, 0.96237801],
      [0.86095238, 1.18309483],
      [1.24285714, 1.48743033],
      [3.92380952, 3.50997283],
      [0.70952381, 1.06833416],
      [0.3752381, 0.6776403],
      [3.08380952, 3.56953042],
      [1.90666667, 2.16219759],
      [0.96, 1.17526897],
      [1.3352381, 1.62480656],
      [1.11333333, 1.40729844],
      [1.16571429, 1.40582655],
      [2.02285714, 2.04134676],
      [0.99714286, 1.23635079],
      [4.03333333, 4.53009338],
      [5.13809524, 4.13382979],
      [1.04, 1.23532143],
      [1.37904762, 1.81557784],
      [7.87238095, 7.21934564],
      [1.19904762, 1.46170124],
      [0.78, 1.05206373],
      [1.0447619, 1.26675295],
      [0.59047619, 0.86734667],
      [0.93904762, 1.12972776],
      [1.36095238, 1.6682466],
      [1.80857143, 2.01104637],
      [1.60380952, 1.83385744],
      [2.30666667, 2.23420452],
      [2.17142857, 2.38270412],
      [2.66952381, 2.42336342],
      [13.28095238, 9.81418405],
      [0.61333333, 0.89390174],
      [1.05619048, 1.30116606],
      [1.57142857, 1.97584047],
      [4.09619048, 3.50853175],
      [10.6447619, 9.99973791],
      [2.14952381, 2.17639386],
      [3.84761905, 3.53453965],
      [3.2847619, 2.8327769],
      [1.06095238, 1.38636182],
      [16.26, 12.02323227],
      [0.36095238, 0.65914739],
      [2.13904762, 2.15090215],
      [3.66666667, 2.87805049],
      [0.95238095, 1.2685271],
      [2.15142857, 2.26209403],
      [1.28380952, 1.50662814],
      [0.64, 0.96309026],
      [0.9352381, 1.26174339],
      [0.72285714, 0.98867464],
      [1.24571429, 1.55490184],
      [0.51238095, 0.78700097],
      [1.63904762, 1.93173967],
      [4.32285714, 3.68853454],
      [1.13428571, 1.3005734],
      [5.86380952, 4.9222067],
      [0.75428571, 1.0648522],
      [1.88761905, 2.15157507],
      [1.31047619, 1.45431534],
      [15.59333333, 10.65705013],
      [2.20666667, 2.11213125],
      [0.79142857, 1.14281309],
      [1.22761905, 1.41104455],
      [0.60761905, 0.89567163],
      [1.74857143, 1.99943393],
      [2.0752381, 2.08145009],
      [0.62761905, 0.89517173],
      [0.85904762, 1.19089364],
      [3.3752381, 3.54214842],
      [0.96095238, 1.13615854],
      [1.72, 1.77648127],
      [1.21809524, 1.47391807],
      [4.69142857, 4.5069283],
      [2.37142857, 2.22948699],
      [3.4047619, 2.87103369],
      [1.58190476, 1.68054872],
      [1.52285714, 1.5374535],
      [1.35238095, 1.49590522],
      [1.76095238, 2.08280924],
      [1.04857143, 1.64592509],
      [10.55047619, 7.51581347],
      [0.5352381, 0.7995214],
      [0.93333333, 1.17162213],
      [0.45809524, 0.7538196],
      [3.36761905, 3.51298259],
      [0.83238095, 1.21518137],
      [0.48095238, 0.77589828],
      [1.22285714, 1.45597339],
      [0.68857143, 0.98423906],
      [2.09047619, 1.97951528],
      [0.42571429, 0.73075676],
      [1.22, 1.39324219],
      [12.25333333, 11.3481953],
      [1.48380952, 1.76317483],
      [1.59238095, 1.74997661],
      [1.7047619, 1.85587258],
      [1.81904762, 2.12907739],
      [18.19238095, 14.83134194],
      [0.99047619, 1.18879882],
      [0.52380952, 0.79276776],
      [0.68095238, 0.96516421],
      [1.3752381, 1.61721577],
      [1.07904762, 1.28078809],
      [0.87904762, 1.13709699],
      [3.72857143, 3.56512778],
      [0.65333333, 0.97091677],
      [0.64190476, 0.95983937],
      [3.70666667, 3.1874237],
      [1.84, 2.13833936],
      [1.31619048, 1.52971356],
      [4.89619048, 4.07894533],
      [1.2647619, 1.48529176],
      [4.30095238, 3.45844708],
      [1.36761905, 1.51566394],
      [1.21904762, 1.44934108],
      [9.02190476, 8.60755461],
      [3.69142857, 3.17280451],
      [2.45714286, 2.44337227],
      [2.51047619, 2.74787966],
      [1.94857143, 2.04781483],
      [1.35142857, 1.5446753],
      [2.94761905, 2.61572517],
      [0.71142857, 0.95496839],
      [1.62666667, 1.84330276],
      [1.50666667, 1.87018037],
      [0.99809524, 1.26566233],
      [5.89428571, 4.94553822],
      [1.03619048, 1.28011337],
      [0.72571429, 0.99666382],
      [0.7447619, 1.10047444],
      [21.35142857, 16.6236585],
      [5.85047619, 5.03317799],
      [1.45904762, 1.64830044],
      [1.78952381, 1.94849998],
      [2.96380952, 2.68019629],
      [2.70857143, 2.81608016],
      [1.38380952, 1.6739815],
      [1.41809524, 1.70026782],
      [18.66761905, 13.31254328],
      [1.52952381, 1.96967117],
      [1.45238095, 1.63652333],
      [7.5552381, 5.69572072],
      [6.54761905, 5.54188545],
      [1.53619048, 1.78168587],
      [12.49047619, 9.40508499],
      [1.4552381, 1.78736203],
      [5.83714286, 4.90938377],
      [3.77047619, 3.22426585],
      [0.69714286, 0.9490303],
      [0.94666667, 1.12250287],
      [8.13714286, 6.07714308],
      [0.76666667, 0.97971489],
      [1.13428571, 1.37808207],
      [2.11619048, 2.64031469],
      [1.05047619, 1.25823942],
      [3.29047619, 2.9078555],
      [1.47714286, 1.69310969],
      [1.29619048, 1.50075688],
      [3.55904762, 4.00462103],
      [6.88380952, 6.22329149],
      [1.55809524, 1.68323843],
      [0.88190476, 1.2578728],
      [0.84952381, 1.07889925],
      [0.42095238, 0.72107],
      [1.11333333, 1.46043481],
      [5.10285714, 4.3557178],
      [5.00190476, 4.22532007],
      [0.50380952, 0.77489462],
      [0.92, 1.20367691],
      [0.56285714, 0.86208575],
      [1.3247619, 1.61016302],
      [1., 1.23056317],
      [6.63619048, 6.00785397],
      [1.25714286, 1.46858255],
      [0.59238095, 1.09871745],
      [3.13619048, 2.74460939],
      [1.08666667, 1.28326716],
      [1.47142857, 1.72809911],
      [10.55142857, 10.53093958],
      [1.16285714, 1.32706315],
      [0.55238095, 0.83015721],
      [1.48, 1.70324061],
      [2.66190476, 2.96759903],
      [1.06, 1.33567176],
      [0.87238095, 1.22781224],
      [2.17238095, 2.38802549],
      [4.2952381, 3.44928926],
      [0.97333333, 1.37156613],
      [5.28761905, 4.55768728],
      [0.99142857, 1.28838136],
      [5.53714286, 4.93028072],
      [3.85238095, 3.29474604],
      [0.54857143, 0.81183089],
      [0.89238095, 1.1104864],
      [1.84857143, 2.04355499],
      [2.61619048, 2.32180651],
      [1.25714286, 1.45031078],
      [1.3152381, 1.44639328],
      [1.89333333, 1.97894954],
      [1.38761905, 1.62223803],
      [3.90095238, 3.9887579],
      [2.13238095, 2.39455713],
      [0.7152381, 0.99897862],
      [1.4847619, 1.81587082],
      [0.80857143, 1.16722008],
      [1.11904762, 1.38361812],
      [1.43619048, 1.59290832],
      [1.36190476, 1.5977592],
      [17.95238095, 13.56085469],
      [7.24666667, 6.67569246],
      [1.58857143, 1.87521754],
      [2.42666667, 2.43522823],
      [1.02095238, 1.21166515],
      [0.94190476, 1.20216321],
      [3.63619048, 3.26292664],
      [1.81809524, 1.9746624],
      [0.80571429, 0.99922147],
      [1.4647619, 1.65477554],
      [7.5447619, 7.51699004],
      [0.61809524, 0.85180827],
      [2.22285714, 2.3278544],
      [5.67619048, 4.83892006],
      [7.17904762, 6.85011915],
      [6.58571429, 6.56068998],
      [1.39904762, 1.57776184],
      [0.54380952, 0.88181888],
      [1.4047619, 1.60563774],
      [2.57047619, 2.83604096],
      [1.09809524, 1.3728303],
      [0.61619048, 0.94934803],
      [1.03619048, 1.22459897],
      [0.53047619, 0.83239457],
      [1.5352381, 1.66738923],
      [1.24666667, 1.52162611],
      [0.82571429, 1.02590527],
      [0.89238095, 1.12919653],
      [1.25047619, 1.41324509],
      [2.91142857, 2.45158119],
      [4.29333333, 4.26610711],
      [0.57047619, 0.80904692],
      [0.36285714, 0.70288618],
      [0.99333333, 1.28617311],
      [1.70095238, 1.87923466],
      [1.60571429, 1.93532714],
      [0.96571429, 1.22931096],
      [4.48380952, 3.76376518],
      [5.33047619, 4.92537886],
      [2.85142857, 2.84401173],
      [1.56285714, 1.72056232],
      [1.18857143, 1.44375927],
      [1.45142857, 1.6382403],
      [1.51428571, 1.66884212],
      [0.70761905, 0.93158291],
      [1.10285714, 1.31326254],
      [1.92095238, 2.18052834],
      [5.39809524, 4.51743671],
      [1.73904762, 1.93968555],
      [1.32380952, 1.566891],
      [1.0752381, 1.28250357],
      [0.82666667, 1.04309673],
      [21.30952381, 18.71740735],
      [0.94285714, 1.32685501],
      [1.81714286, 1.93732546],
      [4.5847619, 3.79811434],
      [1.50380952, 1.65169026],
      [2.10095238, 2.03430087],
      [1.75904762, 2.09671831],
      [2.0552381, 2.29176501],
      [0.87428571, 1.18699529],
      [16.54285714, 13.26982152],
      [2.39142857, 2.61718185],
      [0.50095238, 0.75498284],
      [1.57047619, 1.79344966],
      [2.01238095, 2.30585983],
      [0.82190476, 1.09095685],
      [1.17428571, 1.43065277],
      [0.59714286, 0.86111303],
      [3.30571429, 3.80461696],
      [3.6952381, 3.16926853],
      [7.29714286, 6.46008193],
      [3.4047619, 3.014722],
      [0.70666667, 1.0418177],
      [2.70666667, 3.13955613],
      [2.57428571, 2.83661531],
      [1.91619048, 2.23108446],
      [1.61238095, 1.79396848],
      [1.51714286, 1.85014579],
      [1.2647619, 1.5639986],
      [2.08190476, 2.24307531],
      [4.62, 4.34745465],
      [0.97142857, 1.26834833],
      [1.02952381, 1.317315],
      [0.74952381, 1.11403375],
      [1.08952381, 1.23272241],
      [28.45047619, 22.02790275],
      [1.34571429, 1.77911313],
      [1.4352381, 1.66077899],
      [1.36, 1.50045073],
      [0.73619048, 1.00566243],
      [2.32095238, 2.53148536],
      [5.32857143, 4.99558308],
      [5.4247619, 4.89280874],
      [4.54571429, 3.97256133],
      [2.06, 2.18636903],
      [2.45047619, 2.24200692],
      [1.92190476, 2.19427745],
      [1.95047619, 2.38276046],
      [1.20190476, 1.76451698],
      [0.8352381, 1.12397941],
      [4.96666667, 4.04166789],
      [0.66761905, 0.9793093],
      [4.76952381, 3.92819173],
      [0.59333333, 0.84478982],
      [6.24666667, 5.54145531],
      [2.02857143, 1.95159101],
      [2.16666667, 2.08208534],
      [1.27238095, 1.45736195],
      [1.29428571, 1.50809426],
      [2.35714286, 2.64126078],
      [0.63238095, 0.91448609],
      [4.91047619, 4.72133806],
      [1.88285714, 2.01862552],
      [0.54, 0.8158081],
      [1.49904762, 1.71450369],
      [4.12285714, 3.64771963],
      [2.67619048, 2.39879032],
      [1.52095238, 1.7495684],
      [1.03904762, 1.24724422],
      [0.75142857, 1.12889447],
      [3.32, 3.54944395],
      [6.58285714, 6.51235393],
      [2.18380952, 2.03996376],
      [5.94285714, 5.10653843],
      [5.71809524, 5.08317076],
      [1.68380952, 1.92898314],
      [3.40571429, 2.94734194],
      [2.08190476, 2.46285862],
      [1.25142857, 1.48247633],
      [1.31428571, 1.55664037],
      [1.67619048, 2.10530051],
      [0.61428571, 0.87220886],
      [2.81142857, 3.13120437],
      [7.90761905, 7.33147416],
      [0.64, 1.01696093],
      [0.42095238, 0.72370676],
      [9.4847619, 8.51319738],
      [12.81714286, 9.37488496],
      [1.64666667, 1.89080818],
      [0.47238095, 0.76575456],
      [1.42380952, 1.64357472],
      [2.09238095, 2.26443287],
      [1.58571429, 1.6604032],
      [6.52952381, 5.48188226],
      [2.9047619, 3.5604275],
      [2.2047619, 2.25789118],
      [6.59809524, 5.67101202],
      [0.53428571, 0.87137861],
      [2.76666667, 2.97862757],
      [1.83904762, 1.81522635],
      [1.47333333, 1.80521819],
      [2.53714286, 2.43291936],
      [4.52761905, 3.96817087],
      [0.43809524, 0.72537425],
      [1.66857143, 1.98198142],
      [5.65142857, 4.69492359],
      [1.24, 1.44536039],
      [5.57714286, 5.15994199],
      [2.27809524, 2.4457536],
      [4.15047619, 3.5846531],
      [1.22666667, 1.47165386],
      [3.69142857, 3.08236503],
      [0.44761905, 0.76196428],
      [1.45809524, 1.61205465],
      [1.20380952, 1.52764517],
      [1.14380952, 1.33783221],
      [11.38380952, 9.14811064],
      [0.67142857, 0.94397481],
      [1.22571429, 1.47793476],
      [0.8447619, 1.17501619],
      [2.54285714, 2.96625007]])

        scalingX_test = []
        for i in range(X_test.shape[1]):
            eachFeature = X_test[:, i]
            v = (eachFeature - storeMeanSD[i][0]) / (storeMeanSD[i][1])
            scalingX_test.append(v)

        X_test = np.array(scalingX_test).T

        # from sklearn.externals import joblib
        # model = joblib.load('/home/mrz/MyDrive/Education/Bioinformatics/Hotspots/MLApproach/Optimal/iRecModel.pkl')
        import pickle
        with open(picklePath, 'rb') as pickleFile:
            model = pickle.load(pickleFile)

        # y_artificial = model.predict(X_test)
        #
        # for each in y_artificial:
        #     if each == 1:
        #         print('hotspot')
        #     else:
        #         print('coldspot')

        HC = model.predict_proba(X_test)

        Types = []
        Probabilities = []

        for zero, one in HC:
            if zero>one:
                if zero>=0.70:
                    # print('Type: {}, Probability: {:0.4f}'.format('coldspot', zero))
                    # Types.append('coldspot')
                    # Probabilities.append(zero)
                    return 'coldspot', zero
                else:
                    # print('Type: {}, Probability: {:0.4f}'.format('None', zero))
                    # Types.append('None')
                    # Probabilities.append(zero)
                    return 'None', zero
            else:
                if one>=0.70:
                    # print('Type: {}, Probability: {:0.4f}'.format('hotspot', one))
                    # Types.append('hotspot')
                    # Probabilities.append(one)
                    return 'hotspot', one
                else:
                    # print('Type: {}, Probability: {:0.4f}'.format('None', one))
                    # Types.append('None')
                    # Probabilities.append(one)
                    return 'None', one

        return Types, Probabilities

    if __name__ == '__main__':

        # print('Enter the stepSize:')
        # stepSize = int(input())  # 50 to 200
        # stepSize = 1000  # 50 to 200

        k = 0   # if k=0 then k-mer won't work
        g = 15  # if g=-1 then g-gap won't work # [0 to 8] are the best

        SequencesNumbers = []
        Types = []
        Probabilities = []

        for x, y in zip(X, Y):

            sequencesNumber = []
            types = []
            probabilities = []

            begin = 1
            end = stepSize

            for i in range(0, len(x), stepSize):
                t = []


                # print('{} to {}'.format(begin, end))
                sequencesNumber.append((begin, end))

                begin = end+1
                end = begin + stepSize - 1


                generateFeatures(k, g, x[i:i + stepSize], y)

                t = np.array(t)
                t = np.reshape(t, (1, -1))
                X_test = t[:, :-1]

                X_test = X_test[:, [9954, 30769, 30710, 12413, 13005, 40130, 35526, 4949, 1610, 18664, 20716, 813, 33333, 38464, 22625,
                            20755, 40214, 35206, 7711, 32559, 2535, 38322, 9325, 21984, 2985, 7903, 6863, 32069, 32588, 18839, 30459,
                            5848, 8421, 13508, 37409, 34394, 24020, 33987, 27838, 11912, 37378, 1486, 21097, 12169, 1993, 31502,
                            32006, 32638, 4697, 23060, 23479, 20219, 20483, 11249, 24933, 22411, 37332, 4021, 18452, 20211, 6387,
                            25548, 13650, 7367, 11380, 16047, 18128, 13069, 23433, 29403, 40324, 10603, 9261, 14066, 7131, 1712,
                            17305, 35091, 31996, 7139, 3231, 13060, 2414, 3659, 25992, 1701, 17295, 2405, 8612, 23847, 18410, 22611,
                            33263, 29993, 15424, 36977, 25151, 19953, 7868, 10503, 3273, 31944, 27362, 24247, 1266, 4053, 24520,
                            24757, 10901, 38219, 9012, 15383, 37173, 34539, 12218, 32396, 32760, 8414, 4042, 2498, 29431, 2477,
                            36947, 36948, 12550, 1228, 15439, 32703, 19031, 14059, 27281, 27279, 32711, 33722, 5586, 27544, 40354,
                            6337, 14282, 21480, 14527, 26020, 1252, 29682, 11792, 26262, 21203, 30633, 31328, 2360, 21414, 14606,
                            9672, 38379, 28653, 9671, 19731, 8484, 11910, 22009, 9223, 9666, 3123, 13154, 16136, 22442, 6219, 40070,
                            40167, 13547, 25846, 5673, 14700, 26959, 31432, 1039, 9767, 38867, 12362, 6199, 37097, 28228, 33411,
                            1813, 14636, 33791, 1069, 4239, 24897, 1072, 24900, 4696, 8076, 10732, 12349, 1076, 10745, 23764, 18025,
                            33785, 5741, 18473, 38859, 26343, 4733, 22528, 23798, 29736, 10802, 40233, 8566, 3067, 30767, 23420,
                            19098, 13095, 11039, 3214, 31478, 10630, 2930, 14582, 20657, 6270, 33860, 34838, 31486, 20649, 12280,
                            34833, 26443, 7964, 32645, 8584, 20680, 6255, 9125, 24848, 36377, 18435, 4788, 38338, 14719, 20709,
                            19777, 16163, 34985, 11207, 25819, 19122, 11966, 9138, 25945, 452, 4166, 18428, 8557, 32135, 14235, 7998,
                            14130, 19628, 14508, 4704, 25365, 18690, 5269, 2584, 2585, 12905, 37277, 19426, 13955, 31089, 10184,
                            24607, 25497, 28480, 15042, 20403, 12652, 16393, 6051, 3425, 19366, 35374, 26674, 22175, 28421, 27090,
                            7664, 5134, 22785, 22155, 24710, 27721, 31719, 23959, 23076, 11498, 6512, 30304, 30386, 21575, 13340,
                            2022, 17703, 18213, 5104, 16504, 6581, 14332, 30407, 30597, 778, 22177, 31679, 32236, 36690, 35516,
                            15191, 11694, 20308, 5331, 4385, 4490, 18897, 16372, 24162, 27904, 39092, 20296, 17698, 10090, 5177,
                            12661, 28442, 35473, 30487, 1568, 21795, 15966, 31202, 27671, 32258, 7786, 8708, 7487, 12196, 5442,
                            22213, 40435, 14301, 30121, 8681, 14930, 7701, 16313, 3727, 6997, 4427, 23555, 33156, 40660, 34745, 3298,
                            22737, 23669, 20229, 6675, 7437, 11622, 13406, 36304, 36607, 9313, 3924, 3856, 7831, 27979, 8886, 22730,
                            26695, 14005, 21851, 32781, 32969, 6552, 20976, 31159, 8769, 2010, 32955, 23097, 36106, 24003, 4533,
                            7796, 11050, 33012, 6648, 21559, 22209, 38571, 24011, 3888, 6596, 31631, 20325, 8233, 14395, 3849, 18661,
                            2543, 39207, 5210, 16933, 39209, 27834, 31828, 1451, 33693, 25587, 12616, 34628]]

                type, probability = evaluateModel(X_test)

                types.append(type)
                probabilities.append(probability)

            SequencesNumbers.append(sequencesNumber)
            Types.append(types)
            Probabilities.append(probabilities)

        return SequencesNumbers, Types, Probabilities

###############################################################################################


###############################################################################################

SequencesNumbers, Types, Probabilities = FASTAtoPrediction(stepSize, X, Y)

C=1
for sequencesNumber, types, probabilities in zip(SequencesNumbers, Types, Probabilities):
    print('Test FASTA #{}:'.format(C))
    C += 1

    for x, y, z in zip(sequencesNumber, types, probabilities):
        # print(x)
        # print(y)
        # print(z)
        print('{} to {}: predicted Type: {}, probality {:0.4f}.'.format(x[0], x[1], y, z))

    print()

print('[Note: If probability < 0.70 then predicted type is assign to \'None\'.]')

###############################################################################################


