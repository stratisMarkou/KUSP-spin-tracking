class multi:

    def __init__(self, bincodes = None, quickScales = None, scalesWithZeroes = None):
        self.scale = [0]*32
        self.bincode = bincodes[:]
        
        if not (quickScales == None):
            for i in range(len(bincodes)):
                self.scale[bincodes[i]] = quickScales[i]

        elif not (scalesWithZeroes == None):
            for b in bincodes:
                self.scale[b] = scalesWithZeroes[b]

                    
    def __str__(self):
        result = ""
        for b in self.bincode:
            result += str((self.scale[b], b)) + " "
        return result

    def __add__(self, other):
        resultBincodes = self.bincode[:]
        resultScales = self.scale[:]
        if type(other) == multi:
            for i in other.bincode:
                if not (i in resultBincodes):
                    resultBincodes.append(i)
                resultScales[i] += other.scale[i]
            return multi(bincodes = resultBincodes, scalesWithZeroes = resultScales)
        
        else:
            if not (0 in resultBincodes):
                resultBincodes.append(0)
            resultScales[0] += other
            return multi(bincodes = resultBincodes, scalesWithZeroes = resultScales)

    def __radd__(self, other):
        resultBincodes = self.bincode[:]
        resultScales = self.scale[:]
        if not (0 in resultBincodes):
            resultBincodes.append(0)
        resultScales[0] += other
        return multi(bincodes = resultBincodes, scalesWithZeroes = resultScales)

        
    def __sub__(self, other):
        resultBincodes = self.bincode[:]
        resultScales = self.scale[:]
        for i in other.bincode:
            if not (i in resultBincodes):
                resultBincodes.append(i)
            resultScales[i] -= other.scale[i]
        return multi(bincodes = resultBincodes, scalesWithZeroes = resultScales)

    def __mul__(self, other):
        resultScales = [0]*32
        if type(other) == multi:
            resultBincodes = []
            for b1 in self.bincode:
                for b2 in other.bincode:
                    newBincode, sign = self.multiplyBincodes(b1,b2)
                    resultScales[newBincode] += sign*self.scale[b1]*other.scale[b2]
                    if not (newBincode in resultBincodes):
                        resultBincodes.append(newBincode)
            return multi(bincodes = resultBincodes, scalesWithZeroes = resultScales)
        
        elif type(float(other)) == float:
            for b in self.bincode:
                resultScales[b] = self.scale[b] * other
            return multi(bincodes = self.bincode[:], scalesWithZeroes = resultScales)
    
    def multiplyBincodes(self, bincode1, bincode2):
        count = 0
        blade  = bincode1^bincode2
        while(bincode1 > 0):
            bincode1 = bincode1 >> 1
            count += str(bin(bincode1 & bincode2)).count("1")
        return (blade,(-1)**count)

    def __rmul__(self, other):
        if type(float(other)) == float:
            resultScales = [0]*32
            for b in self.bincode:
                resultScales[b] = self.scale[b] * other
            return (multi(bincodes = self.bincode[:], scalesWithZeroes = resultScales))

    def applyRotorPair(self, other):
        return other*self*other.tilda()
    
    def getRotorNorm(self):
        return (self*self.tilda()).getBladeScale(0)
    
    def tilda(self):
        resultScale = [0]*32
        for b in self.bincode:
            if b == 0:
                resultScale[b] = self.scale[b]
            else:
                resultScale[b] = -1*self.scale[b]
        return multi(bincodes = self.bincode[:], scalesWithZeroes = resultScale)

    def getMagnitudeOfOrder(self, order):
        result = 0
        for bincode in self.bincode:
            if bin(bincode).count("1") == order:
                result += (self.scale[bincode])**2
        return result**0.5
    
    def getBladeScale(self, bincode):
        return self.scale[bincode]

