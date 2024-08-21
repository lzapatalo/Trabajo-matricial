#include <iostream>
#include <string.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <stdlib.h>

using namespace std;
using namespace Eigen;

const int dof = 2; //cada nodo tiene movimiento en x y movimiento en y

//Prototipos de las funciones (lo que hace cada una está en la definición de cada una al final)
int pos_ini(int i);
int pos_fin(int i, int dof);
Matrix<double, 2*dof, 2*dof> K_cercha(double e, double a, double l, double t);
void imprimirK (MatrixXd matriz);

int main(){
    //Propiedades
    double E = 210*pow(10,9); //[Pa]
    double W = 450*7.9*0.5*9.8; //[N/m] carga aferente
    //Áreas [m^2]
    double A7 = 9.44*pow(10,-4);
    double A8 = 5.3*pow(10,-4);
    double A9 = 2*9.27*pow(10,-4);
    double A10 = 2*7.66*pow(10,-4);
    double A11 = 2*3.4*pow(10,-4);

    const int Nnodos = 29;
    const int Nelementos = 56;

    Vector<double, Nelementos> A;
    A <<
        A7,//1
        A7,//2
        A7,//3
        A8,//4
        A7,//5
        A8,//6
        A7,//7
        A9,//8
        A9,//9
        A9,//10
        A9,//11
        A9,//12
        A9,//13
        A9,//14
        A9,//15
        A9,//16
        A9,//17
        A9,//18
        A9,//19
        A9,//20
        A11,//21
        A11,//22
        A11,//23
        A11,//24
        A11,//25
        A11,//26
        A11,//27
        A11,//28
        A11,//29
        A11,//30
        A11,//31
        A11,//32
        A11,//33
        A11,//34
        A11,//35
        A11,//36
        A11,//37
        A11,//38
        A11,//39
        A11,//40
        A11,//41
        A11,//42
        A11,//43
        A9,//44
        A10,//45
        A10,//46
        A10,//47
        A10,//48
        A10,//49
        A10,//50
        A10,//51
        A10,//52
        A10,//53
        A10,//54
        A10,//55
        A10; //56

    //Se define el vector de nodos como variable global
    MatrixXd nodos(Nnodos, dof);
    nodos <<
        -1.470,2.2653, //1
        0.0,2.4347, //2
        1.470,2.2653, //3
        -8.640,0.500, //4
        -7.240,0.6613, //5
        -5.770,0.8306, //6
        -4.355,0.9936, //7
        -2.885,1.163, //8
        -1.470,1.326, //9
        0.0,1.4954, //10
        1.470,1.326, //11
        2.885,1.163, //12
        4.355,0.9936, //13
        5.770,0.8306, //14
        7.240,0.6613, //15
        8.640,0.500, //16
        -8.640,0.0, //17
        -7.240,0.0, //18
        -5.770,0.0, //19
        -4.355,0.0, //20
        -2.885,0.0, //21
        -1.470,0.0, //22
        0.0,0.0, //23
        1.470,0.0, //24
        2.885,0.0, //25
        4.355,0.0, //26
        5.770,0.0, //27
        7.240,0.0, //28
        8.640,0.0;  //29

    //Numero de grados de libertad total de la estructura = #elementos de la matriz / bytes por elemento
    const int num_Dof = Nnodos *dof; 

    MatrixXi elementos(Nelementos, 2);
    elementos <<
        1,2, //1
        2,3, //2
        1,9, //3
        1,10, //4
        2,10, //5
        10,3, //6
        3,11, //7
        4,5, //8
        5,6, //9
        6,7, //10
        7,8, //11
        8,9, //12
        9,10, //13
        10,11, //14
        11,12, //15
        12,13, //16
        13,14, //17
        14,15, //18
        15,16, //19
        4,17, //20
        4,18, //21
        5,18, //22
        5,19, //23
        6,19, //24
        6,20, //25
        7,20, //26
        7,21, //27
        8,21, //28
        8,22, //29
        9,22, //30
        9,23, //31
        10,23, //32
        23,11, //33
        11,24, //34
        24,12, //35
        12,25, //36
        25,13, //37
        13,26, //38
        26,14, //39
        14,27, //40
        27,15, //41
        15,28, //42
        28,16, //43
        16,29, //44
        17,18, //45
        18,19, //46
        19,20, //47
        20,21, //48
        21,22, //49
        22,23, //50
        23,24, //51
        24,25, //52
        25,26, //53
        26,27, //54
        27,28, //55
        28,29;  //56

    Matrix<int, Nnodos, dof> restricciones;
    restricciones <<
        0,0, //1
        0,0, //2
        0,0, //3
        1,1, //4
        0,0, //5
        0,0, //6
        0,0, //7
        0,0, //8
        0,0, //9
        0,0, //10
        0,0, //11
        0,0, //12
        0,0, //13
        0,0, //14
        0,0, //15
        1,1, //16
        1,1, //17
        0,0, //18
        0,0, //19
        0,0, //20
        0,0, //21
        0,0, //22
        0,0, //23
        0,0, //24
        0,0, //25
        0,0, //26
        0,0, //27
        0,0, //28
        1,1; //29

    Vector<int, Nnodos*dof> rest_Dof;
    for(int i=0 ; i<Nnodos ; i++){
        for(int j=0; j<dof; j++){
            rest_Dof(i*dof + j) = restricciones(i,j);
        }
    }

    cout<<"Datos recibidos"<<endl;
    
    //Buscar índices en el arreglo completo de la estructura, correspondientes a los grados de libertad libres y restringidos
    vector<int> free_index;
    vector<int> rest_index;
    for (int i = 0; i < Nnodos*dof; i++) {
        if (rest_Dof(i) == 0) {
            free_index.push_back(i);
        }
        else {
            rest_index.push_back(i);
        }
    }

    //Matrices globales elementos
    Vector<double, Nelementos> L;
    Vector<double, Nelementos> theta;
    vector<Matrix<double, 2*dof, 2*dof>> K_elem;

    cout<<"Matrices creadas"<<endl;

    for(int n=0 ; n<Nelementos ; n++){
        //Se extrae el índice de los nodos y se corrige para poder acceder en la lista de coordenadas
        int ni = elementos(n,0) - 1;
        int nj = elementos(n,1) - 1;

        double xi = nodos(ni,0);
        double yi = nodos(ni,1);
        double xj = nodos(nj,0);
        double yj = nodos(nj,1);

        L(n) = pow((pow(xj-xi, 2) + pow(yj-yi, 2)), 0.5);
        theta(n) = atan2(yj-yi, xj-xi);

        K_elem.push_back(K_cercha(E , A(n) , L(n) , theta(n)));
    }

    cout<<"Listas matrices de cada elemento"<<endl;

    //Matriz global estructura
    Matrix<double, num_Dof, num_Dof> K;
    K.setZero(num_Dof, num_Dof);

    //Se definen las "submatrices" de la matriz de rigidez
    const int f = free_index.size();
    const int r = rest_index.size();
    MatrixXd Kaa(f,f); Kaa.setZero();
    MatrixXd Kbb(r,r); Kbb.setZero();
    MatrixXd Kab(f,r); Kab.setZero();
    MatrixXd Kba(r,f); Kba.setZero();

    cout<<"Creadas las submatrices"<<endl;

    for (int n=0 ; n<Nelementos ; n++){ //Para cada elemento "elem" realizar ensamble
        int ni = elementos(n,0) - 1;
        int nj = elementos(n,1) - 1;

        int pos_ini_i = pos_ini(ni); //Posicion inicial del dof del nodo i
        int pos_fin_i = pos_fin(ni,dof); //Posicion final del dof del nodo i
        int pos_ini_j = pos_ini(nj); //Posicion inicial del dof del nodo j
        int pos_fin_j = pos_fin(nj,dof); //Posicion final del dof del nodo j

        Matrix<double, 2*dof, 2*dof> matriz = K_elem[n];

        //i con i (slicing)
        K.block(pos_ini_i, pos_ini_i, dof, dof) += matriz.block(0, 0, dof, dof);
        //i con j
        K.block(pos_ini_i, pos_ini_j, dof, dof) += matriz.block(0, dof, dof, dof);
        //j con i
        K.block(pos_ini_j, pos_ini_i, dof, dof) += matriz.block(dof, 0, dof, dof);
        //j con j
        K.block(pos_ini_j, pos_ini_j, dof, dof) += matriz.block(dof, dof, dof, dof);
    }
    cout<<"Ensamble realizado"<<endl;   

    //List comprehension
    int row_free = 0;
    for (const auto& i : free_index) {
        int col_free = 0;
        for (const auto& j : free_index) {
            Kaa(row_free, col_free) = K(i, j);
            col_free++;
        }
        row_free++;
    }
    
    // Asignar índices de fila y columna para los grados de libertad restringidos
    int row_rest = 0;
    for (const auto& i : rest_index) {
        int col_rest = 0;
        for (const auto& j : rest_index) {
            Kbb(row_rest, col_rest) = K(i, j);
            col_rest++;
        }
        row_rest++;
    }
    
    // Asignar índices de fila y columna para los grados de libertad libres y restringidos
    row_free = 0;
    for (const auto& i : free_index) {
        int col_rest = 0;
        for (const auto& j : rest_index) {
            Kab(row_free, col_rest) = K(i, j);
            col_rest++;
        }
        row_free++;
    }
    
    // Asignar índices de fila y columna para los grados de libertad restringidos y libres
    row_rest = 0;
    for (const auto& i : rest_index) {
        int col_free = 0;
        for (const auto& j : free_index) {
            Kba(row_rest, col_free) = K(i, j);
            col_free++;
        }
        row_rest++;
    }

    cout<<"Submatrices definidas"<<endl;

    //Vector de fuerzas EXTERNAS de toda la estructura (sin tener en cuenta reacciones)
    Matrix<double, Nnodos, dof> Q;
    Q <<
        0,-W*(0.750+L(0)/2), //1
        0,-W*(L(0)/2+L(1)/2), //2
        0,-W*(0.750+L(2)/2), //3
        0,0, //4
        0,-W*(L(7)/2+L(8)/2), //5
        0,-W*(L(8)/2+L(9)/2), //6
        0,-W*(L(9)/2+L(10)/2), //7
        0,-W*(L(10)/2+L(11)/2), //8
        0,-W*(L(11)/2), //9
        0,0, //10
        0,-W*(L(14)/2), //11
        0,-W*(L(14)/2+L(15)/2), //12
        0,-W*(L(15)/2+L(16)/2), //13
        0,-W*(L(16)/2+L(17)/2), //14
        0,-W*(L(17)/2+L(18)/2), //15
        0,0, //16
        0,0, //17
        0,0, //18
        0,0, //19
        0,0, //20
        0,0, //21
        0,0, //22
        0,0, //23
        0,0, //24
        0,0, //25
        0,0, //26
        0,0, //27
        0,0, //28
        0,0; //29
    
    Vector<double, Nnodos*dof> Q_vector;
    for(int i=0 ; i<Nnodos ; i++){
        for(int j=0; j<dof; j++){
            Q_vector(i*dof + j) = Q(i,j);
        }
    }

    VectorXd Qa(f); //Fuerzas externas en los nodos
    int pos=0;
    for(const auto& q : free_index){
        Qa(pos) = Q_vector(q);
        pos++;
    }

    VectorXd qb(r); //Desplazamientos conocidos

    VectorXd qa(f);
    qa = Kaa.colPivHouseholderQr().solve(Qa); //Hallamos desplazamientos desconocidos
    cout<<"Desplazamientos hallados"<<endl;

    VectorXd Qb(r);
    Qb = Kba * qa; //Hallamos reacciones
    cout<<"Reacciones halladas"<<endl;

    cout<<"Matriz global:"<<endl;
    cout<<K<<endl;
    imprimirK(K);

    cout<<"Reacciones:"<<endl;
    cout<<Qb<<endl;

    cout<<"Desplazamientos:"<<endl;
    cout<<qa<<endl;

    system("pause");
    return 0;
}

Matrix<double, 2*dof, 2*dof> K_cercha(double e, double a, double l, double t){
    double c = cos(t);
    double pc = pow(c, 2);
    double s = sin(t);
    double ps = pow(s, 2);
    Matrix<double, 2*dof, 2*dof> k; 
    k <<
        pc, c*s, -pc, -c*s,
        c*s, ps, -c*s, -ps,
        -pc, -c*s,pc, c*s,
        -c*s, -ps,c*s, ps;

    Matrix<double, 2*dof, 2*dof> kcercha = k*(e*a/l);
    
    return kcercha;
}

//Halla la posición inicial del nodo dentro de la matriz K completa
int pos_ini(int i){ 
    return dof*i;
}

//Halla la posición inicial del nodo dentro de la matriz K completa
int pos_fin(int i, int dof){
    return dof*i + (dof-1);
}

void imprimirK (MatrixXd matriz){
    cout<<"|";
    for(int n=0;n<matriz.cols();n++){
        cout<<"--";
    }
    cout<<"|"<<endl;
    for(int i=0;i<matriz.rows();i++){
        cout<<"|";
        for(int j=0;j<matriz.cols();j++){
            if(matriz(i,j) == 0){
                cout<<"  ";
            }
            else{
                cout<<"# ";
            }
        }
        cout<<"|"<<endl;
    }
    cout<<"|";
    for(int m=0;m<matriz.cols();m++){
        cout<<"--";
    }
    cout<<"|"<<endl;
}