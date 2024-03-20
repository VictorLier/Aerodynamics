import numpy as np
import matplotlib.pyplot as plt

# From xfoil
alphaL0 = -4.0
Cd_friction = 0.004


def Cl(alpha, AR):
    Cl = 2 * np.pi / (1 + 2 / AR) * (alpha - alphaL0)

    return Cl

def Cdi(Cl, AR, alpha):
    Cl = 2 * np.pi / (1 + 2 / AR) * (alpha - alphaL0)
    Cdi = 1/(np.pi * AR) * Cl**2

    return Cdi


AR = np.array([4, 6, 8, 10, 100000])

alpha = np.linspace(-6, 10, 20)

plt.figure()
for i in range(len(AR)):
    plt.plot(alpha, Cl(alpha, AR[i]), label='AR = {0}'.format(AR[i]))
plt.xlabel('Angle of Attack (alpha)')
plt.ylabel('Coefficient of Lift (Cl)')
plt.title('Coefficient of Lift')
plt.legend()

plt.figure()
for i in range(len(AR)):
    plt.plot(alpha, Cdi(Cl(alpha, AR[i]), AR[i], alpha), label='AR = {0}'.format(AR[i]))
plt.xlabel('Angle of Attack (alpha)')
plt.ylabel('Coefficient of Drag (Cdi)')
plt.title('Coefficient of Drag')
plt.legend()

plt.show()



def
