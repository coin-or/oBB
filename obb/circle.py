# Circle class
class circle:

    lbound = 0
    lReduced = None
    uReduced = None

    def __init__(self, xcc, rcc):
        self.xc = xcc
        self.r = rcc
 
    # Replicate python2 object comparison in python3
    # (only used for breaking ties in priority queue)  
    def __lt__(self, other):
        return hex(id(self)) < hex(id(other))
