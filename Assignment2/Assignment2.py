import numpy as np

Theta = np.linspace(0, np.pi, 10)


def coordlength(theta):
    coordlength = 1 * np.sin(theta)
    return coordlength


def ai(theta: npdarray, A: npdarray) -> float:
    n = np.arrange(1, len(A))+1
    ai = np.sum(n * A * np.sin(n * theta) / np.sin(theta))
    
    return ai

def b(AR: npdarray) -> npdarray:
    b = AR * Pi / 4
    return b

def Cl(theta: npdarray, A: npdarray, AR) -> float:
    c0 = coordlength(theta)
    b = b(AR)

    cl = 4*b/c0 * np.sum(A * np.sin(n * theta))
    
    return Cl

def Cdi(Cl, alpha):
    cdi = Cl * ai