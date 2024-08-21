import math
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
fig,ax = plt.subplots()
fig,estructura = plt.subplots()

E = 210e6 #[KPa]
b = 0.03 #[m]
h = 0.15 #[m]

A = b*h #[m]^2

dof = 2 #cada nodo tiene movimiento en x y movimiento en y

nodos = np.array([
  [6 ,7],
  [13,7],
  [7 ,3],
  [3 ,3],
  [6 ,5],
  [10,5],
  [9 ,9],
])

num_Dof = int(dof * nodos.size / 2) #Numero de grados de libertad total de la estructura

elementos = np.array([
  [4,1],
  [4,3],
  [4,5],
  [1,7],
  [7,2],
  [7,6],
  [5,7],
  [1,6],
  [5,2],
  [5,1],
  [6,2],
  [3,6],
  [5,3],
])

#Vector de fuerzas EXTERNAS de toda la estructura
Q = np.array([
  ['R',0],
  [0,-25],
  [0,0],
  ['R','R'],
  [0,0],
  [30*math.sin(55*math.pi/180), -30*math.cos(55*math.pi/180)],
  [0,'R']
])

Q = Q.flatten()

q = np.array([
  [0,'I'],
  ['I','I'],
  ['I',',I'],
  [0,0],
  ['I','I'],
  ['I','I'],
  ['I',0]
])
#Este no se usará, ya que se verá mejor reflejado en la matriz de restricciones

restricciones = np.array([
  [1,0],
  [0,0],
  [0,0],
  [1,1],
  [0,0],
  [0,0],
  [0,1]
])

rest_Dof = restricciones.flatten()

#Buscar índices en el arreglo completo de la estructura, correspondientes a los grados de libertad libres y restringidos
rest_index = np.where(rest_Dof != 0)[0]
free_index = np.where(rest_Dof == 0)[0]

#Matrices globales elementos
L = []
theta = []
K_elem = []

def K_cercha(e,a,l,t):
  c = math.cos(t)
  s = math.sin(t)

  k = np.array([
  [c**2,c*s,-c**2,-c*s],
  [c*s,s**2,-c*s,-s**2],
  [-c**2,-c*s,c**2,c*s],
  [-c*s,-s**2,c*s,s**2],
  ])

  return (e*a/l)*k

for e in elementos:

  #Se extrae el índice de los nodos y se corrige para poder acceder en la lista de coordenadas
  ni = e[0] - 1
  nj = e[1] - 1

  xi, yi = nodos[ni][0], nodos[ni][1]
  xj, yj = nodos[nj][0], nodos[nj][1]

  estructura.plot([xi,xj],[yi,yj],'ko-') #Grafica elemento por elemento

  L.append(math.sqrt((xj-xi)**2 + (yj-yi)**2))
  theta.append(math.atan2(yj-yi, xj-xi))

  K_elem.append(K_cercha(E,A,math.sqrt((xj-xi)**2 + (yj-yi)**2),math.atan2(yj-yi, xj-xi)))

#Matriz global estructura
K = np.zeros((num_Dof, num_Dof))

def pos_ini(i): #Valido solo para cercha
  return 2*i

def pos_fin(i,dof):
  return dof*i + (dof-1)

for e,elem in enumerate(elementos): #Para cada elemento "elem" realizar ensamble

  ni = elem[0] - 1
  nj = elem[1] - 1

  pos_ini_i , pos_ini_j = pos_ini(ni) , pos_ini(nj) #Posicion inicial del dof de los nodos i y j
  pos_fin_i , pos_fin_j = pos_fin(ni,dof) , pos_fin(nj,dof) #Posicion final del dof de los nodos i y j

  #i con i (slicing)
  K[pos_ini_i:pos_fin_i+1 , pos_ini_i:pos_fin_i+1] += K_elem[e][0:2, 0:2] #Valido solo para cercha
  #i con j
  K[pos_ini_i:pos_fin_i+1 , pos_ini_j:pos_fin_j+1] += K_elem[e][0:2, 2:4] #Valido solo para cercha
  #j con i
  K[pos_ini_j:pos_fin_j+1 , pos_ini_i:pos_fin_i+1] += K_elem[e][2:4, 0:2] #Valido solo para cercha
  #j con j
  K[pos_ini_j:pos_fin_j+1 , pos_ini_j:pos_fin_j+1] += K_elem[e][2:4, 2:4] #Valido solo para cercha

#List comprehension
Kaa = [[K[i][j] for j in free_index] for i in free_index]
Kbb = [[K[i][j] for j in rest_index] for i in rest_index]
Kab = [[K[i][j] for j in rest_index] for i in free_index]
Kba = [[K[i][j] for j in free_index] for i in rest_index]

Qa = [float(Q[i]) for i in free_index]
qb = []

Kaa = np.array(Kaa)
Kbb = np.array(Kbb)
Kab = np.array(Kab)
Kba = np.array(Kba)

qa = np.linalg.solve(Kaa,Qa) #Desplazamientos
Qb = np.matmul(Kba,qa) #Reacciones

#print(Kaa)

ax.spy(K, markersize=3)

estructura.set_aspect('equal','box')
plt.show()

print(qa)
print(Qb)
