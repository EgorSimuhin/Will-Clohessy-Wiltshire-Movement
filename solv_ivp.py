import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
def func(t,y):
    dydt=np.sin(t)+y;
    return dydt
x2 = np.linspace(0,3,301) # интервал поиска решения
y1 = [2] # Значения функции в точке х2[0]
y = solve_ivp(func,t_span =(0,3),y0=y1,t_eval=x2)

yt=-np.sin(x2)/2-np.cos(x2)/2+5/2*np.exp(x2)

plt.plot(x2,y.y.T,':r*',x2,yt,'--g',lw=3,ms=4)
plt.legend(['Численное решение','Точное решение'])
p=np.max(abs((yt-y.y)))
print("Максимальная погрешность метода solve_ivp равна %.8f" % p)
plt.show()
