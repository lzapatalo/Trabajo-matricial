import math
import numpy as np
import matplotlib.pyplot as plt

dof=2 #En "x" y en "y"
E = 210e9 #[Pa]
W = 450*7.9*0.5*9.8 #[N/m] carga aferente
#Àreas #m^2
A7 = 9.44*10**(-4)
A8 = 5.3*10**(-4)
A9 = 2*9.27*10**(-4)
A10 = 2*7.66*10**(-4)
A11 = 2*3.4*10**(-4)


A = [
  A7,#1
  A7,#2
  A7,#3
  A8,#4
  A7,#5
  A8,#6
  A7,#7
  A9,#8
  A9,#9
  A9,#10
  A9,#11
  A9,#12
  A9,#13
  A9,#14
  A9,#15
  A9,#16
  A9,#17
  A9,#18
  A9,#19
  A9,#20
  A11,#21
  A11,#22
  A11,#23
  A11,#24
  A11,#25
  A11,#26
  A11,#27
  A11,#28
  A11,#29
  A11,#30
  A11,#31
  A11,#32
  A11,#33
  A11,#34
  A11,#35
  A11,#36
  A11,#37
  A11,#38
  A11,#39
  A11,#40
  A11,#41
  A11,#42
  A11,#43
  A9,#44
  A10,#45
  A10,#46
  A10,#47
  A10,#48
  A10,#49
  A10,#50
  A10,#51
  A10,#52
  A10,#53
  A10,#54
  A10,#55
  A10 #56
  ]

nodos = np.array([
  [-1.470,2.2653],#1
  [0,2.4347],#2
  [1.470,2.2653],#3
  [-8.640,0.500],#4
  [-7.240,0.6613],#5
  [-5.770,0.8306],#6
  [-4.355,0.9936],#7
  [-2.885,1.163],#8
  [-1.470,1.326],#9
  [0,1.4954],#10
  [1.470,1.326],#11
  [2.885,1.163],#12
  [4.355,0.9936],#13
  [5.770,0.8306],#14
  [7.240,0.6613],#15
  [8.640,0.500],#16
  [-8.640,0],#17
  [-7.240,0],#18
  [-5.770,0],#19
  [-4.355,0],#20
  [-2.885,0],#21
  [-1.470,0],#22
  [0,0],#23
  [1.470,0],#24
  [2.885,0],#25
  [4.355,0],#26
  [5.770,0],#27
  [7.240,0],#28
  [8.640,0] #29
])

num_DOF = nodos.size   #nodos.size es 58

elementos = np.array([
  [1,2],#1
  [2,3],#2
  [1,9],#3
  [1,10],#4
  [2,10],#5
  [10,3],#6
  [3,11],#7
  [4,5],#8
  [5,6],#9
  [6,7],#10
  [7,8],#11
  [8,9],#12
  [9,10],#13
  [10,11],#14
  [11,12],#15
  [12,13],#16
  [13,14],#17
  [14,15],#18
  [15,16],#19
  [4,17],#20
  [4,18],#21
  [5,18],#22
  [5,19],#23
  [6,19],#24
  [6,20],#25
  [7,20],#26
  [7,21],#27
  [8,21],#28
  [8,22],#29
  [9,22],#30
  [9,23],#31
  [10,23],#32
  [23,11],#33
  [11,24],#34
  [24,12],#35
  [12,25],#36
  [25,13],#37
  [13,26],#38
  [26,14],#39
  [14,27],#40
  [27,15],#41
  [15,28],#42
  [28,16],#43
  [16,29],#44
  [17,18],#45
  [18,19],#46
  [19,20],#47
  [20,21],#48
  [21,22],#49
  [22,23],#50
  [23,24],#51
  [24,25],#52
  [25,26],#53
  [26,27],#54
  [27,28],#55
  [28,29] #56
  ])

L = []

for e in elementos:

  #Se extrae el índice de los nodos y se corrige para poder acceder en la lista de coordenadas
  ni = e[0] - 1
  nj = e[1] - 1

  xi, yi = nodos[ni][0], nodos[ni][1]
  xj, yj = nodos[nj][0], nodos[nj][1]
  L.append(math.sqrt((xj-xi)**2 + (yj-yi)**2)) 

#Vector de fuerzas EXTERNAS de toda la estructura
Q = np.array([
  [0,-W*(0.750+L[0]/2)],#1
  [0,-W*(L[0]/2+L[1]/2)],#2
  [0,-W*(0.750+L[2]/2)],#3
  ['R','R'],#4
  [0,-W*(L[7]/2+L[8]/2)],#5
  [0,-W*(L[8]/2+L[9]/2)],#6
  [0,-W*(L[9]/2+L[10]/2)],#7
  [0,-W*(L[10]/2+L[11]/2)],#8
  [0,-W*(L[11]/2)],#9
  [0,0],#10
  [0,-W*(L[14]/2)],#11
  [0,-W*(L[14]/2+L[15]/2)],#12
  [0,-W*(L[15]/2+L[16]/2)],#13
  [0,-W*(L[16]/2+L[17]/2)],#14
  [0,-W*(L[17]/2+L[18]/2)],#15
  ['R','R'],#16
  ['R','R'],#17
  [0,0],#18
  [0,0],#19
  [0,0],#20
  [0,0],#21
  [0,0],#22
  [0,0],#23
  [0,0],#24
  [0,0],#25
  [0,0],#26
  [0,0],#27
  [0,0],#28
  ['R','R'] #29
])

#Al convertir en un np.array, cuando hay strings y enteros, pasa todos los datos los pasa a string
Q=np.array(Q).flatten()

#1 es restringido, 0 es libre. Dónse hayan 1, habrán incógnitas
restricciones = np.array([
  [0,0],#1
  [0,0],#2
  [0,0],#3
  [1,1],#4
  [0,0],#5
  [0,0],#6
  [0,0],#7
  [0,0],#8
  [0,0],#9
  [0,0],#10
  [0,0],#11
  [0,0],#12
  [0,0],#13
  [0,0],#14
  [0,0],#15
  [1,1],#16
  [1,1],#17
  [0,0],#18
  [0,0],#19
  [0,0],#20
  [0,0],#21
  [0,0],#22
  [0,0],#23
  [0,0],#24
  [0,0],#25
  [0,0],#26
  [0,0],#27
  [0,0],#28
  [1,1] #29
])

#Buscar índices en el arreglo completo
rest_DoF = restricciones.flatten() #Para leer la matriz como un vector que se lee en orden de las filas
rest_index = np.where(rest_DoF != 0)[0] #[0,6,7,13], se escoje la posición [0] porque el np.where incluye una lista de tuplas con el primer elemento de cada tupla igual a la posición encontrada
free_index = np.where(rest_DoF == 0)[0]

theta = []
K_elem = []

def K_cercha(e,a,l,t):
    c=math.cos(t)
    s=math.sin(t)

    k=np.array([
        [c**2,c*s,-c**2,-c*s],
        [c*s,s**2,-c*s,-s**2],
        [-c**2,-c*s,c**2,c*s],
        [-c*s,-s**2,c*s,s**2],
    ])
    return (e*a/l)*k

fig, estructura = plt.subplots()

L=[]
cont=0
for e in elementos:
    ni=e[0]-1   #Le restamos 1 para volver a la convención inicial de nodos
    nj=e[1]-1

    xi,yi=nodos[ni][0], nodos[ni][1]
    xj,yj=nodos[nj][0], nodos[nj][1]

    estructura.plot([xi,xj],[yi,yj])

    L.append(math.sqrt((xj-xi)**2 + (yj-yi)**2))
    theta.append(math.atan2(yj-yi,xj-xi))
    K_elem.append(K_cercha(E, A[cont], L[-1], theta[-1]))
    cont+=1


#Matriz global de la estructura
K=np.zeros([num_DOF,num_DOF])

def pos_ini(i):  #Válido solo para cercha
    return 2*i   #Sería dof*i
def pos_fin(i):  #Solo para cercha
    return 2*i+1

for e,elem in enumerate(elementos): #Para cada elemento realizar ensamble
    ni = elem[0]-1
    nj = elem[1]-1

    pos_ini_i, pos_ini_j = pos_ini(ni),pos_ini(nj) #Pos.incial del dof del nodo i
    pos_fin_i, pos_fin_j = pos_fin(ni),pos_fin(nj)

    # i con i (slicing)
    K[pos_ini(ni): pos_fin(ni)+1,pos_ini(ni): pos_fin(ni)+1] += K_elem[e][0:2,0:2] #Sería de 0 a dof, solo para cercha
    K[pos_ini(ni): pos_fin(ni)+1,pos_ini(nj): pos_fin(nj)+1] += K_elem[e][0:2,2:4] #Solo en cercha
    K[pos_ini(nj): pos_fin(nj)+1,pos_ini(ni): pos_fin(ni)+1] += K_elem[e][2:4,0:2] #Solo en cercha
    K[pos_ini(nj): pos_fin(nj)+1,pos_ini(nj): pos_fin(nj)+1] += K_elem[e][2:4,2:4] #Solo en cercha

print(K[32][30])
#List comprehension
Kaa=[[K[i][j] for j in free_index] for i in free_index]
Kbb=[[K[i][j] for j in rest_index] for i in rest_index]
Kab=[[K[i][j] for j in rest_index] for i in free_index]
Kba=[[K[i][j] for j in free_index] for i in rest_index]

Qa = [float(Q[i]) for i in free_index]

Kaa=np.array(Kaa)
Kbb=np.array(Kbb)
Kba=np.array(Kba)
Kab=np.array(Kab)

## Resolver un sistema matricial del tipo Ax=b, con (A,b)
qa = np.linalg.solve(Kaa,Qa)
Qb = np.matmul(Kba,qa)

print(qa)  #Desplazamientos desconocidos
print(Qb)  #Reacciones

fig, ax = plt.subplots()
ax.spy(K, markersize=3)


estructura.set_aspect('equal','box')
plt.show()
