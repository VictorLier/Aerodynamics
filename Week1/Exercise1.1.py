import numpy as np
import matplotlib.pyplot as plt

Alpha = np.array([-2, 0, 2, 3, 6, 8, 10, 12, 14])
Cl = np.array([0.05, 0.25, 0.44, 0.64, 0.85, 1.08, 1.26, 1.43, 1.56])
Cd = np.array([0.006, 0.006, 0.006, 0.007, 0.0075, 0.0092, 0.0115, 0.0150, 0.0186])
Cm = np.array([-0.042, -0.040, -0.038, -0.036, -0.036, -0.036, -0.034, -0.030, -0.025])

Alpha = Alpha / 180*np.pi

#print(Alpha)


xcp  = (Cd  * np.sin(Alpha ) + Cl  * np.cos(Alpha ) - 4*Cm )/(4*(Cl *np.cos(Alpha ) + Cd  * np.sin(Alpha )))

print(xcp)

plt.plot(xcp, Alpha)
plt.show()
