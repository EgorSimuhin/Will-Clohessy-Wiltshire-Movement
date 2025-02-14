import numpy as np                                                                            
import matplotlib.pyplot as plt

R1 = 7 * 10**6
gamma = 398600.4415 * 1000**3
omega_per = (gamma / R1**3) ** 0.5


#Функция для правой части системы ДУ. Решение в ИСО, в базисе I
def F_iso(t, y):
    position = y[0:3]
    velocity = y[3:6]
    radius = np.linalg.norm(position)
    result = np.zeros(6)
    result[0:3] = velocity
    result[3:6] = -gamma / (radius**3) * position
    return result


#Функция для правой части системы ДУ. Решение в НИСО, в базисе J
def F_niso(t, y):
    position = y[0:3]
    velocity = y[3:6]
    ax = 3 * omega_per**2 * position[0] + 2 * omega_per * velocity[1]
    ay = -2 * omega_per * velocity[0]
    az = -omega_per**2 * position[2]
    result = np.zeros(6)
    result[0:3] = velocity
    result[3:6] = [ax, ay, az] 
    return result
    
#Шаг метода Рунге_Кутта
def calc_step(t, y, delta_t, func):
    k_1 = func(t, y)
    k_2 = func(t + delta_t / 2, y + delta_t / 2 * k_1)
    k_3 = func(t + delta_t / 2, y + delta_t / 2 * k_2)
    k_4 = func(t + delta_t, y + delta_t * k_3)
    return y + delta_t / 6 * (k_1 + 2 * (k_2 + k_3) + k_4)


#Функция для расчета в НИСО относительного вектора
def solve_ballistics_niso(y_0, delta_t, final_time, func):
    current_y = y_0
    current_t = 0
    solutions = [current_y]
    times = [current_t]
    while current_t < final_time:
        current_y = calc_step(current_t, current_y, delta_t, func)
        solutions.append(current_y)
        current_t += delta_t
        times.append(current_t)
    return solutions, times


#Функция для расчета относительного вектора в ИСО
def solve_ballistics_iso(y_01, y_02, delta_t, final_time, func):
    current_y1 = y_01
    current_y2 = y_02
    current_t = 0
    solutions_y1 = [current_y1]
    solutions_y2 = [current_y2]
    times = [current_t]
    r_otn_i_0 = y_02[0:3] - y_01[0:3]
    r_otn_j = [r_otn_i_0]
    while current_t < final_time:
        current_y1 = calc_step(current_t, current_y1, delta_t, func)
        current_y2 = calc_step(current_t, current_y2, delta_t, func)
        r_otn_i = current_y2[0:3] - current_y1[0:3]
        j1 = current_y1[0:3] / np.linalg.norm(current_y1[0:3])
        j3 = np.cross(current_y1[0:3], current_y1[3:6]) / np.linalg.norm(np.cross(current_y1[0:3], current_y1[3:6]))
        j2 = np.cross(j3, j1)
        A = np.array([j1, j2, j3])
        #A = A.transpose()
        r_otn_j.append(np.dot(A, r_otn_i))
        solutions_y1.append(current_y1)
        solutions_y2.append(current_y2)
        current_t += delta_t
        times.append(current_t)
    return r_otn_j, solutions_y1, solutions_y2, times


r_otn_j = np.array([10, 10, 10])
v_otn_j = np.array([0.1, 0.1, 0.1])
delta_t = 5 
final_time = 2*24*3600


start_data_2_i = np.array([r_otn_j[0] + R1, r_otn_j[1], r_otn_j[2], v_otn_j[0] - omega_per * r_otn_j[1], v_otn_j[1] + omega_per * (r_otn_j[0] + R1), v_otn_j[2]])
start_data_1_i = np.array([R1, 0, 0, 0, omega_per * R1, 0])

r_otn_iso, solutions_y1, solutions_y2, t = solve_ballistics_iso(start_data_1_i, start_data_2_i, delta_t, final_time, F_iso)

y_0 = np.zeros(6)
y_0[0:3] = r_otn_j[:]
y_0[3:6] = v_otn_j[:]
solutions, times = solve_ballistics_niso(y_0, delta_t, final_time, F_niso)

solutions = np.array(solutions)
rx_niso = solutions[:, 0]
ry_niso = solutions[:, 1]
rz_niso = solutions[:, 2]

r_otn_iso = np.array(r_otn_iso)
rx_iso = r_otn_iso[:, 0]
ry_iso = r_otn_iso[:, 1]
rz_iso = r_otn_iso[:, 2]

print("max error ry", max(ry_niso - ry_iso))
print("max error rx",max(rx_niso - rx_iso))
print("max error rz" ,max(rz_niso - rz_iso))

fig1 = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
ax.plot(rx_niso, ry_niso, rz_niso, label='НИСО')
ax.plot(rx_iso, ry_iso, rz_iso, label='ИСО')
ax.legend()

ax.set_xlabel('ρ_x (m)')
ax.set_ylabel('ρ_y (m)')
ax.set_zlabel('ρ_z (m)')
ax.set_title('Траектория второго спутника относительно первого в базисе j')
plt.show()
