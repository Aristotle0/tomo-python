import math

def nextpow2(i):
    """
    Find the next power of two
    """
    buf = math.ceil(math.log(i)/math.log(2))
    return int(math.pow(2, buf))
