
from ngsolve import BSpline

def GetBHCurve():
    B = [0, 1.413, 1.594, 1.751, 1.839, 1.896, 1.936, 1.967, 2.008, 2.042, 2.073, 2.101, 2.127, 2.151, 2.197, 2.240, 2.281, 2.321, 2.361, 2.400, 2.472]
    h = [0, 00.400, 00.801, 01.601, 02.402, 03.203, 04.003, 04.804, 06.405, 08.007, 09.608, 11.210, 12.811, 14.412, 17.615, 20.818, 24.020, 27.223, 30.426, 33.629, 39.634]
    H = [1e4*hi for hi in h]
    return BSpline(2, [0] + B, H)



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    bh = GetBHCurve()
    b = np.linspace(0, 2.471, 1000)
    h = [bh(bi) for bi in b]
    plt.plot(h, b)
    plt.show()
