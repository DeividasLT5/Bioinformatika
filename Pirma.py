from Bio import SeqIO
from Bio.Seq import Seq

def readFasta(filename):
    for seq_record in SeqIO.parse(filename, "fasta"):
        sequence = seq_record.seq
    return sequence

def buildFrames(dnr):
    frame1Segments = []
    frame2Segments = []
    frame3Segments = []
    frames = []
    for readingFrame in [0, 1, 2]:
        position = 0
        for x in dnr:
            segment = dnr[readingFrame:][position:position+3]
            position = position + 3
            if segment == "":
                break
            if  readingFrame == 0:
                frame1Segments.append(segment)
            elif readingFrame == 1:
                frame2Segments.append(segment)
            else:
                frame3Segments.append(segment)
    frames.append(frame1Segments)
    frames.append(frame2Segments)
    frames.append(frame3Segments)
    return frames

def parseSegmentsFromFrames(frameList):
    parsedSegments = []
    for frame in frameList:
        afterStartCodon = False
        temp = []
        for codon in frame:
            if codon == "ATG":
                afterStartCodon = True
            if afterStartCodon:
                temp.append(codon)
            if (codon == "TAA" or codon == "TAG" or codon == "TGA") and afterStartCodon:
                afterStartCodon = False
                parsedSegments.append(temp)
                temp = []
    return parsedSegments

def filterBySize(parsedSegments):
    filteredSegments = []
    for segment in parsedSegments:
        size = 0
        for codon in segment:
            size += len(codon)
        if size >= 100:
           filteredSegments.append(segment)
    return filteredSegments

def buildCodonList():
    defaultCodonList = []
    symbols = ["T", "A", "C", "G"]
    for firstSymbol in symbols:
        for secondSymbol in symbols:
            for thirdSymbol in symbols:
                defaultCodonList.append(firstSymbol + secondSymbol + thirdSymbol)
    # defaultCodonList = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
    #                     "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
    #                     "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
    #                     "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"]
    return defaultCodonList

def codonTotals(defaultCodonList, segmentList):
    total = 0
    freqList = [0] * len(defaultCodonList)
    for segment in segmentList:
        for codon in segment:
            i = 0
            while i < len(defaultCodonList):
                if codon == defaultCodonList[i]:
                    freqList[i] += 1
                    total += 1
                    break
                i += 1
    return freqList, total

def codonFrequencies(defaultCodonList, totalFilteredList, total):
    frequenciesList = [0] * len(defaultCodonList)
    j = 0
    while j < len(frequenciesList):
        frequenciesList[j] = round(totalFilteredList[j] / total, 5)
        j += 1
    return frequenciesList

def buildAcidList():
    defaultCodonList = buildCodonList()
    AcidList = []
    for codon in defaultCodonList:
        Acid = Seq(codon).translate()
        if Acid not in AcidList:
            AcidList.append(Acid) 
    AcidList.remove("*")
    return AcidList

def buildDicodonAcidList():
    codonAcidList = buildAcidList()
    dicodonAcidList = []
    for firstAcid in codonAcidList:
        for secondAcid in codonAcidList:
            dicodonAcidList.append(firstAcid + secondAcid)
    return dicodonAcidList

def dicodonTotals(defaultDicodonList, segmentList):
    freqList = [0] * len(defaultDicodonList)
    previousCodon = ""
    for segment in segmentList:
        for codon in segment:
            i = 0
            codonPair = Seq(previousCodon + codon).translate()
            previousCodon = codon
            while i < len(defaultDicodonList):
                if codonPair == defaultDicodonList[i]:
                    freqList[i] += 1
                    break
                i += 1
    return freqList, sum(freqList)

def dicodonFrequency(defaultDicodonList, totalFilteredList, total):
    frequenciesList = [0] * len(defaultDicodonList)
    j = 0
    while j < len(defaultDicodonList):
        frequenciesList[j] = round(totalFilteredList[j] / total, 5)
        j += 1
    return frequenciesList

def frequenciesFromFile(filename):
    sequence = readFasta(filename)

    frame = buildFrames(sequence)
    reverseFrame = buildFrames(sequence.reverse_complement())
    genomeFrames = frame + reverseFrame

    codingSegments = parseSegmentsFromFrames(genomeFrames)
    filteredSegments = filterBySize(codingSegments)

    codonTotalList, codonTotal = codonTotals(buildCodonList(), filteredSegments)
    codonFrequencyList = codonFrequencies(buildCodonList(), codonTotalList, codonTotal)

    dicodonTotalList, dicodonTotal = dicodonTotals(buildDicodonAcidList(), filteredSegments)
    dicodonFrequencyList = dicodonFrequency(buildDicodonAcidList(), dicodonTotalList, dicodonTotal)
    return codonFrequencyList, dicodonFrequencyList

def compareFrequencies(file1, file2):
    codonFrequencies1, dicodonFrequencies1 = frequenciesFromFile(file1)
    codonFrequencies2, dicodonFrequencies2 = frequenciesFromFile(file2)

    codonNumber = len(buildCodonList())
    codonDistances = [0] * codonNumber
    i = 0
    while i < codonNumber:
        frequency1 = codonFrequencies1[i]
        frequency2 = codonFrequencies2[i]
        if frequency1 == 0 and frequency2 == 0:
                codonDistances[i] = 1
                i += 1
                continue
        codonDistances[i] = abs(frequency1 - frequency2) / max(frequency1, frequency2)
        i += 1

    dicodonNumber = len(buildDicodonAcidList())
    dicodonDistances = [0] * dicodonNumber
    i = 0
    while i < dicodonNumber:
        frequency1 = dicodonFrequencies1[i]
        frequency2 = dicodonFrequencies2[i]
        if frequency1 == 0 and frequency2 == 0:
            dicodonDistances[i] = 1
            i += 1
            continue
        dicodonDistances[i] = abs(frequency1 - frequency2) / max(frequency1, frequency2)
        i += 1    

    return codonDistances, dicodonDistances    

def writeResultsToFiles(codonDistances, dicodonDistances, name1, name2):
    file = open("distances/" + name1[:10] + name2[:10] + ".txt", "w")
    codonNames = buildCodonList()
    acidNames = buildAcidList()

    file.write(str(len(codonNames)) + "\n")
    i = 0
    while i < len(codonNames):
        file.write(str(codonNames[i]) + "   " + str("%.3f" % codonDistances[i]) + "\n")
        i += 1

    file.write("-------------------------------------------------------------------\n         ")
    i = 0
    while i < len(acidNames):
        file.write(str(acidNames[i]) + "         ")
        i += 1

    file.write("\n" + str(len(acidNames)))

    i = 0
    while i < len(acidNames):
        file.write("\n" + str(acidNames[i]))
        j = 0
        while j < len(acidNames):
            file.write("   " + str("%.5f" % dicodonDistances[j + (i * 20)]))
            j += 1
        i += 1        


def main():
    filenames = ["bacterial1.fasta", "bacterial2.fasta", "bacterial3.fasta", "bacterial4.fasta",
                 "mamalian1.fasta", "mamalian2.fasta", "mamalian3.fasta", "mamalian4.fasta"]

    offset = 1
    for firstFilename in filenames:
        for secondFilename in filenames[offset:]:
            print("Comparing " + firstFilename + " and " + secondFilename + "...")
            codonFrq, dicodonFrq = compareFrequencies(firstFilename, secondFilename)
            writeResultsToFiles(codonFrq, dicodonFrq, firstFilename, secondFilename)
        offset += 1

if __name__ == "__main__":
    main()

                                                                        
