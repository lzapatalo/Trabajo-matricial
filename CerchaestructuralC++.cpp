//Este código funciona si todos los desplazamientos conocidos son cero. Añadir caso en el que hay asentamientos.
//Crear código para elementos viga sometidos a axial, cortante y flexión (es decir, crear los dos sistemas de superposición).
//Sobre esto último, añadir momentos de empotramiento para el sistema empotrado.
//Incluir códigos archivos externos que contengan los métodos que permitan crear las matrices de rigidez globales pora cada caso. Usar clases.
//Incluir gráficas de la estructura inicial y la estructura deformada en los nodos. El k.spy y las reacciones se pueden incluir en un archivo .txt. El archivo puede contener conclusiones sobre los resultados.  

#include<iostream>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

int dof = 2;
double E = 210*pow(10,9);
double W = 450*7.9*0.5*9.8;

double A7 = 9.44*pow(10,-4);
double A8 = A8 = 5.3*pow(10,-4);
double A9 = A9 = 2*9.27*pow(10,-4);
double A10 = A10 = 2*7.66*pow(10,-4);
double A11 = A11 = 2*3.4*pow(10,-4);

//A cada elemento se le atribuye un área de acuerdo a su sección

double A[] = {A7,A7,A7,A8,A7,A8,A7,A9,A9,A9,A9,A9,A9,A9,A9,A9,A9,A9,A9,A9,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11,A9,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10,A10 
};

vector<vector<double>> nodos = {
    {-1.470, 2.2653},
    {0, 2.4347},
    {1.470, 2.2653},
    {-8.640, 0.500},
    {-7.240, 0.6613},
    {-5.770, 0.8306},
    {-4.355, 0.9936},
    {-2.885, 1.163},
    {-1.470, 1.326},
    {0, 1.4954},
    {1.470, 1.326},
    {2.885, 1.163},
    {4.355, 0.9936},
    {5.770, 0.8306},
    {7.240, 0.6613},
    {8.640, 0.500},
    {-8.640, 0},
    {-7.240, 0},
    {-5.770, 0},
    {-4.355, 0},
    {-2.885, 0},
    {-1.470, 0},
    {0, 0},
    {1.470, 0},
    {2.885, 0},
    {4.355, 0},
    {5.770, 0},
    {7.240, 0},
    {8.640, 0}
};

int num_DOF = nodos.size()*nodos[0].size();                   

vector<vector<int>> elementos = {{1,2},
                      {2,3},
                      {1,9},
                      {1,10},
                      {2,10},
                      {10,3},
                      {3,11},
                      {4,5},
                      {5,6},
                      {6,7},
                      {7,8},
                      {8,9},
                      {9,10},
                      {10,11},
                      {11,12},
                      {12,13},
                      {13,14},
                      {14,15},
                      {15,16},
                      {4,17},
                      {4,18},
                      {5,18},
                      {5,19},
                      {6,19},
                      {6,20},
                      {7,20},
                      {7,21},
                      {8,21},
                      {8,22},
                      {9,22},
                      {9,23},
                      {10,23},
                      {23,11},
                      {11,24},
                      {24,12},
                      {12,25},
                      {25,13},
                      {13,26},
                      {26,14},
                      {14,27},
                      {27,15},
                      {15,28},
                      {28,16},
                      {16,29},
                      {17,18},
                      {18,19},
                      {19,20},
                      {20,21},
                      {21,22},
                      {22,23},
                      {23,24},
                      {24,25},
                      {25,26},
                      {26,27},
                      {27,28},
                      {28,29}
};
int num_ele = elementos.size();

vector<int> where( vector<vector<int>>matriz, bool what) {
    vector<int> places;
    for (int i = 0; i < num_DOF/dof; i++) {
        for (int j = 0; j < dof; j++) {
            if ((matriz[i][j] && what) || (!matriz[i][j] && !what)) {
                places.push_back(i * dof + j);
            }
        }
    }
    return places;
}
     
vector<vector<double>> K_cercha (double e, double a, double l, double t) {
    double c = cos(t);
    double s = sin(t);

    int n = 2 * dof;
    int m = 2 * dof;

    vector<vector<double>> K_el_glob(n,vector<double>(m)); //Si se crea sin definir, entonces hay que ponerle diensiones.

    K_el_glob[0][0] = pow(c, 2);
    K_el_glob[0][1] = c * s;
    K_el_glob[0][2] = -pow(c, 2);
    K_el_glob[0][3] = -c * s;
    K_el_glob[1][0] = c * s;
    K_el_glob[1][1] = pow(s, 2);
    K_el_glob[1][2] = -(c * s);
    K_el_glob[1][3] = -pow(s,2);
    K_el_glob[2][0] = -pow(c, 2);
    K_el_glob[2][1] = -(c * s);
    K_el_glob[2][2] = pow(c, 2);
    K_el_glob[2][3] = c * s;
    K_el_glob[3][0] = -(c * s);
    K_el_glob[3][1] = -pow(s, 2);
    K_el_glob[3][2] = c * s;
    K_el_glob[3][3] = pow(s,2);

    for (int i = 0; i < 2 * dof; i++) {
        for (int j = 0; j < 2 * dof; j++) {
            K_el_glob[i][j] = (double)((e*a)/l) * K_el_glob[i][j];
        }
    }

    return K_el_glob;
}    

int pos_ini(int i){
    return dof*i;
}
int pos_fin(int j){
    return (dof*j)+1;
}

vector<vector<int>> restricciones = { {0, 0},
                           { 0, 0 },
                           { 0, 0 },
                           { 1, 1 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 1, 1 },
                           { 1, 1 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 0, 0 },
                           { 1, 1 }};

int main(){
    double L[num_ele];
    int ni, nj; 
    double xi, yi, xj, yj;
    for(int i=0; i<num_ele; i++){
        ni = elementos[i][0] -1;
        nj = elementos[i][1] -1;
        xi = nodos[ni][0];
        yi = nodos[ni][1];
        xj = nodos[nj][0];
        yj = nodos[nj][1];
        L[i] = pow((pow(xj-xi,2)+ pow(yj-yi,2)),0.5);
    }

    double Q[] = {
        0,-W*(0.750+L[0]/2),
        0,-W*(L[0]/2+L[1]/2),
        0,-W*(0.750+L[2]/2),
        1,1,
        0,-W*(L[7]/2+L[8]/2),
        0,-W*(L[8]/2+L[9]/2),
        0,-W*(L[9]/2+L[10]/2),
        0,-W*(L[10]/2+L[11]/2),
        0,-W*(L[11]/2),
        0,0,
        0,-W*(L[14]/2),
        0,-W*(L[14]/2+L[15]/2),
        0,-W*(L[15]/2+L[16]/2),
        0,-W*(L[16]/2+L[17]/2),
        0,-W*(L[17]/2+L[18]/2),
        1,1,
        1,1,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        0,0,
        1,1
    };


    vector<int> rest_index = where(restricciones, true);
    vector<int> free_index = where(restricciones, false);

    double theta[num_ele];
    vector<vector<vector<double>>> K_elem;

    for (int i = 0; i < num_ele; i++) {
        ni = elementos[i][0] - 1;
        nj = elementos[i][1] - 1;
        xi = nodos[ni][0];
        yi = nodos[ni][1];
        xj = nodos[nj][0];
        yj = nodos[nj][1];
        theta[i] = atan2(yj - yi, xj - xi);
        K_elem.push_back(K_cercha(E, A[i], L[i], theta[i])); //es como .append
    }

    //Ensamble de la matriz global:

    vector<vector<double>> K(num_DOF, vector<double>(num_DOF)); 
    cout<<"Num_DOF= "<<num_DOF<<endl;

    //Se crea la matriz de ceros:
    for (int i=0; i<num_DOF; i++){
        for (int j=0; j<num_DOF; j++){
            K[i][j] = 0;
        }
    }

    for(int i=0; i<num_ele; i++){
        ni = elementos[i][0] -1;
        nj = elementos[i][1] -1;
        int pos_ini_i = pos_ini(ni);
        int pos_ini_j = pos_ini(nj);
        int pos_fin_i = pos_fin(ni);
        int pos_fin_j = pos_fin(nj);


        //Slicing:

        int cont_a = 0;
        int cont_b = 0;
        for (int a = pos_ini_i; a < pos_fin_i+1; a++) {
        cont_b = 0; // Reiniciar cont_b
            for (int b = pos_ini_i; b < pos_fin_i+1; b++) {
                K[a][b] += K_elem[i][cont_a][cont_b];
                cont_b++;
            }
        cont_a++;    
        }
        
        cont_a = 0;
        for (int a = pos_ini_i; a < pos_fin_i+1; a++) {
            cont_b = dof; // Reiniciar cont_b
            for (int b = pos_ini_j; b < pos_fin_j+1; b++) {
                K[a][b] += K_elem[i][cont_a][cont_b];
                cont_b++;
            }
            cont_a++;    
        }
            
        cont_a = dof;
        for (int a = pos_ini_j; a < pos_fin_j+1; a++) {
            cont_b = 0; // Reiniciar cont_b
            for (int b = pos_ini_i; b < pos_fin_i+1; b++) {
                K[a][b] += K_elem[i][cont_a][cont_b];
                cont_b++;
            }
            cont_a++;    
        }
            
        cont_a = dof;
        for (int a = pos_ini_j; a < pos_fin_j+1; a++) {
            cont_b = dof; // Reiniciar cont_b
            for (int b = pos_ini_j; b < pos_fin_j+1; b++) {
                K[a][b] += K_elem[i][cont_a][cont_b];
                cont_b++;
            }
            cont_a++;    
        }
    }   
         //Se organiza la matriz de rigidez:

    vector<vector<double>> Kaa;

    // Recorrer las filas seleccionadas
    for (int i : free_index) {
        // Crear una nueva fila para Kaa
        vector<double> row;
        // Recorrer las columnas seleccionadas
        for (int j : free_index) {
            // Agregar el elemento correspondiente de K a Kaa
            row.push_back(K[i][j]);
        }
        // Agregar la fila a Kaa
        Kaa.push_back(row);
    }

    vector<vector<double>> Kbb;

    for (vector<int>::iterator i = rest_index.begin(); i!=rest_index.end(); i++) { //i and j are iterators, they are similar to pointers.
        vector<double> row;
        for (vector<int>::iterator j = rest_index.begin(); j!=rest_index.end(); j++) {
            row.push_back(K[*i][*j]); //*i accesses the information i is pointing at. Remember: & means where
        }
        Kbb.push_back(row);
    }

    vector<vector<double>> Kab;

    for (int i : free_index) { //Another way to iterate through free_index and rest_index
        vector<double> row;
        for (int j : rest_index) {
            row.push_back(K[i][j]);
        }
        Kab.push_back(row);
    }
    
    vector<vector<double>> Kba;

    for (int i : rest_index) {
        vector<double> row;
        for (int j : free_index) {
            row.push_back(K[i][j]);
        }
        Kba.push_back(row);
    }

    vector<double> Qa;
    for (int i : free_index){
        Qa.push_back(Q[i]);
    }
    
    //Visualizamos la distribución de entradas de la matriz:
    for(int i = 0; i < num_DOF; i++) {
        for (int j = 0; j < num_DOF; j++) {
            if (K[i][j] != 0) {
                cout << "X ";
            } else {
                cout << "_ "; // Espacio en blanco si es nulo
            }
        }
        cout << endl; // Salto de línea después de cada fila
    }
    
    //Volvemos las matrices y vectores objetos operables por Eigen, similar a np.array:

    int sizerest = rest_index.size();
    int sizefree = free_index.size();
   
    // Definir y llenar Kaa como una matriz
    MatrixXd Kaa1(sizefree,sizefree);
    for (int a = 0; a < sizefree; ++a) {
        for (int b = 0; b < sizefree; ++b) {
            // Asignar valores de K[a][b] a Kaa(a, b)
            Kaa1(a, b) = Kaa[a][b];
        }
    }

    MatrixXd Kbb1(sizerest,sizerest);
    for (int a = 0; a < sizerest; ++a) {
        for (int b = 0; b < sizerest; ++b) {
            Kbb1(a, b) = Kbb[a][b];
        }
    }

    MatrixXd Kab1(sizefree,sizerest);
    for (int a = 0; a < sizefree; ++a) {
        for (int b = 0; b < sizerest; ++b) {
            Kab1(a, b) = Kab[a][b];
        }

    }

    MatrixXd Kba1(sizerest,sizefree);
    for (int a = 0; a < sizerest; ++a) {
        for (int b = 0; b < sizefree; ++b) {
            Kba1(a, b) = Kba[a][b];
        }
    }

    VectorXd Qa1(sizefree);
    for (int c = 0; c < sizefree; ++c) {
        // Asignar valores de Qa[c] a Qa(c)
        Qa1(c) = Qa[c];
    }

    VectorXd qa1 = Kaa1.fullPivLu().solve(Qa1);
    VectorXd Qb1 = Kba1 * qa1;

    cout<<"Desplazamientos [mm]: "<<endl << qa1*1000 <<endl;
    cout <<"Reacciones[N]: "<<endl<< Qb1 << endl;
}