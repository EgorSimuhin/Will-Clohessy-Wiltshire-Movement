import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def func(y,x):
    dydx=np.sin(x)+y 
    return dydx

x2 = np.linspace(0,3,301) # интервал поиска решения
y0 = [2] # Значения функции в точке х2[0]
y = odeint(func, y0, x2);
print(y)
yt=-np.sin(x2)/2-np.cos(x2)/2+5/2*np.exp(x2)

plt.plot(x2,y,':r*',x2,yt,'--g',lw=3,ms=4)

plt.legend(['Численное решение','Точное решение'])

p=np.max(abs((yt-y.T)))

print("Максимальная погрешность метода odeint равна %.8f" % p)
plt.show()
