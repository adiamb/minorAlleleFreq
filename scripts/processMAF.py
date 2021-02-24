from collections import defaultdict, Counter

def parsePed(pedFile):
    '''process ped file and make 2 dicts or hashed lists.
    the dxdict stores indexes for diagnosis keyed by diagnois and the rsdict stores the alleles keyed by snp id
    '''
    DXDict = defaultdict(list)
    RsDict = defaultdict(list)
    with open(pedFile, 'rt') as pedIn:
        for n, line in enumerate(pedIn):
            lineParse = line.strip().split(' ')
            Dx = lineParse[5]
            DXDict[Dx].append(n)
            rsIterator = 0
            for rs in range(6, len(lineParse), 2): # since every 2 alleles are one person make an instantiate an iterator 
                rsIterator += 1
                RsDict['snp_{}'.format(rsIterator-1)].append(lineParse[rs]) # a1
                RsDict['snp_{}'.format(rsIterator-1)].append(lineParse[rs+1]) #a2
        return DXDict, RsDict

def getMaf(DXDict, RsDict):
    '''use these dicts to count alleles and write out the frequencies
    '''
    outFile= open('outs/sim1Text.freqparsed', 'w')
    outFile.write('SNP,CASE.MAF,CASE.COUNT,CONTROL.MAF,CONTROL.COUNT,TOTAL.MAF,TOTAL.COUNT'+'\n')
    caseAlleles = len(DXDict.get('2'))*2
    contAlleles = len(DXDict.get('1'))*2
    totalAlleles = caseAlleles+contAlleles
    Alleles = [caseAlleles, contAlleles, totalAlleles]
    for snp, val in RsDict.items():
        print(snp, len(val))
        caseVal = Counter([val[i] for i in DXDict.get('2')])
        controlVal = Counter([val[i] for i in DXDict.get('1')])
        TotalVal = Counter(val)
        sortAlleles = TotalVal.most_common()
        if len(caseVal) == 2 and len(controlVal) ==2 and len(TotalVal) ==2: # discard monomorphic 
            getMinor = sortAlleles[1][0]
            outStr = []
            outStr.append(snp)
            for n, item in enumerate([caseVal, controlVal, TotalVal]):
                count = item.get(getMinor)
                maf = item.get(getMinor)/Alleles[n]
                outStr.append(str(maf))
                outStr.append(str(count))
            outFile.write(','.join(outStr)+'\n')
    outFile.close()

def main():
    pedFile = 'data/sim1Text.ped'
    DXDict, RsDict=parsePed(pedFile)
    getMaf(DXDict, RsDict)

if __name__ == "__main__":
    main()










    
