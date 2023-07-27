function x = FEMJeffcott
% Encuentra la solucion de un modelo de Jeffcott modelado con elemento finito
clc;
clear;

%Definiendo la cantidad de elementos que va tener la malla
n=49;

%La lonjitud del eje
L=0.857;

%El modulo de Young
E=210e9;

%El diametro del eje
de=0.0127;

%La densidad del eje
densidad=7850;

%El area transversal del eje
Ae=(pi*de^2)/4;

%Los valores de alfa y beta se otienen de manera experimental
alfa=0;
beta=0.0000135;

%La velocidad angular del eje
omega=0;%1200*(2*pi)*(1/60);

%Defimi,os la inercia del eje
I=(pi/4)*(de/2)^4;

%Definimos la matriz de rigidez elastica
Ke=(E*I/(L^3))*[12, 0,    0,      6*L,    -12,  0,    0,     6*L;
                0,  12,   -6*L,   0,      0,    -12,  -6*L,  0;
                0,  -6*L, 4*L^2,  0,      0,    6*L,  2*L^2, 0;
                6*L, 0,   0,      4*L^2,  -6*L, 0,    0,     2*L^2,
                -12, 0,   0,      -6*L,   12,   0,    0,     -6*L;
                0,  -12, 6*L,     0,      0,    12,   6*L,   0;
                0,  -6*L, 2*L^2,  0,      0,    6*L,  4*L^2, 0;
                6*L, 0,   0,      2*L^2, -6*L,  0,    0,     4*L^2];
   

%Definimos la matriz de rigidez de amrtiguamiento interno.
Kai=(E*I/(L^3))*[0, 12,  -6*L,   0,     0,    -12,  -6*L,   0;
               -12, 0,   0,      -6*L,  12,   0,    0,      -6*L;
               6*L, 0 ,  0,      4*L^2, -6*L, 0,    0,      2*L^2;
               0,   6*L, -4*L^2, 0,     0,    -6*L, -2*L^2, 0;
               0,   -12, 6*L,    0,     0,    12,   6*L,    0;
               12,  0,   0,      6*L,   -12,  0,    0,      6*L;
               6*L, 0,   0,      2*L^2, -6*L, 0,    0,      4*L^2;
               0,   6*L, -2*L^2, 0,     0,    -6*L, -4*L^2, 0];

%Definimos la matriz de masa              
Me=(densidad*Ae*L/840)*[312,   0,     0,      44*L,   108,    0,      0,      -26*L;
                        0,     312,   -44*L,  0,      0,      108,    26*L,   0;
                        0,     -44*L, 8*L^2,  0,      0,      -26*L,  -6*L^2, 0;
                        44*L,  0,     0,      8*L^2,  26*L,   0,      0,      -6*L^2;
                        108,   0,     0,      26*L,   312,    0,      0,      -44*L;
                        0,     108,   -26*L,  0,      0,      312,    44*L,   0;
                        0,     26*L,  -6*L^2, 0,      0,      44*L,   8*L^2,  0;
                        -26*L, 0,     0,      -6*L^2, -44*L,  0,      0,      8*L^2];

%Definimos la matriz de inercia rotacional
Mir=((densidad*I)/(30*L))*[36,  0 ,   0,     3*L,   -36,  0,    0,    3*L;
                           0,   36,   -3*L,  0,     0,    -36,  -3*L, 0;
                           0,   -3*L, 4*L^2, 0,     0,    3*L,  -L^2, 0;
                           3*L, 0,    0,     4*L^2, -3*L, 0,    0,    -L^2;
                           -36, 0,    0,     -3*L,  36,   0,    0,    -3*L^2;
                           0,   -36,  3*L,   0,     0,    36,   3*L,  0;
                           0,   -3*L, -L^2,  0,     0,    3*L,  4*L^2,0;
                           3*L, 0,    0,     -L^2,  -3*L, 0,    0,    4*L^2];
 
%Definimos la matriz de efecto giroscopio 
Cg=(densidad*I/(15*L))*[0,  -36,  3*L,   0,     0,    36,  3*L,   0;
                      36,   0,    0,     3*L,   -36,  0,   0,     3*L;
                      -3*L, 0,    0,    -4*L^2, 3*L,  0,   0,     L^2;
                      0,    -3*L, 4*L^2, 0,     0,    3*L, -L^2,  0;
                      0,    36,   -3*L,  0,     0,    -36, -3*L,  0;
                      -36,  0,    0,     -3*L,  36,   0,   0      -3*L;
                      -3*L, 0,    0,     L^2,   3*L,  0,   0,     -4*L^2;
                      0,   -3*L,  -L^2,   0,    0,    3*L, 4*L^2, 0];


%definimos el epsilon como la relacion de extentricidad del cojinente.
epsilon=0.9;%L/(n-1);

%Es la fuerza ejercida por el peso del rotor sobre los cojinetes
f_cojinete=15;

%Definimos la holgura o juego raidal entre los mu;onesy los cojinetes
c_cojinete=0.15;

%definimos la velocidad angular del cojinete.
omega_cojinete=omega;

h0_cojinete=1/((pi^2)*(1-epsilon^2)+16*epsilon^2)^(3/2);

%Definiendo las componentes en x e y del componente de rigidez
kxx_cojinete=4*h0_cojinete*((pi^2)*(2-epsilon^2)+16*(epsilon^2));
kxy_cojinete=(h0_cojinete*pi*((pi^2)*(1-epsilon^2)^2+16*(epsilon^4)))/(epsilon*sqrt(1-epsilon^2));
kyx_cojinete=-1*(h0_cojinete*pi*(pi^2*(1-epsilon^2)*(1+2*epsilon^2)+32*(epsilon^2)*(1+epsilon^2)))/(epsilon*sqrt(1-epsilon^2));
kyy_cojinete=4*h0_cojinete*((pi^2)*(1+2*epsilon^2)+(32*epsilon^2*(1+epsilon^2))/(1-epsilon^2));

%definiendo las componentes x e y del amortiguamiento
cxx=h0_cojinete*(2*pi*(1-epsilon^2)^(1/2)*(pi^2)*(1+2*(epsilon^2))-16*(epsilon^2))/epsilon;
cxy=-8*h0_cojinete*(pi^2)*((pi^2)*(1+2*epsilon^2)+16*(epsilon^2));
cyx=cxy;
cyy=h0_cojinete*2*pi*(pi^2*(1-epsilon^2)^2+48*epsilon^2)/(epsilon*sqrt(1-epsilon^2));

%Definiendo las matrices de rigides del cojinete y de rigidez del soporte 

k_cojinetes=zeros(8, 8);
k_soportes=k_cojinetes;
c_cojinetes=k_cojinetes;
g_discos=k_cojinetes;
m_discos=k_cojinetes;

%Inicializo las matrices del rigidez y soporte de los  cojinetes y el amortiguamiento de los cojinetes.
k_cojinetes(1,1)=(f_cojinete/c_cojinete)*(kxx_cojinete);
k_cojinetes(1,2)=(f_cojinete/c_cojinete)*(kxy_cojinete);
k_cojinetes(2,1)=(f_cojinete/c_cojinete)*(kyx_cojinete);
k_cojinetes(2,2)=(f_cojinete/c_cojinete)*(kyy_cojinete); 

k_cojinetes(7,7)=(f_cojinete/c_cojinete)*(kxx_cojinete);
k_cojinetes(7,8)=(f_cojinete/c_cojinete)*(kxy_cojinete);
k_cojinetes(8,7)=(f_cojinete/c_cojinete)*(kyx_cojinete);
k_cojinetes(8,8)=(f_cojinete/c_cojinete)*(kyy_cojinete);         

k_soportes(1,1)=(f_cojinete/c_cojinete)*(kxx_cojinete);
k_soportes(2,2)=(f_cojinete/c_cojinete)*(kyy_cojinete);

k_soportes(7,7)=(f_cojinete/c_cojinete)*(kxx_cojinete);
k_soportes(8,8)=(f_cojinete/c_cojinete)*(kyy_cojinete);
      
%Definiendo la matriz de amortiguamiento del cojinete.     
      
c_cojinetes(1,1)=f_cojinete/(c_cojinete*omega_cojinete)*(cxx);
c_cojinetes(1,2)=f_cojinete/(c_cojinete*omega_cojinete)*(cxy);
c_cojinetes(2,1)=f_cojinete/(c_cojinete*omega_cojinete)*(cyx);
c_cojinetes(2,2)=f_cojinete/(c_cojinete*omega_cojinete)*cyy;  

c_cojinetes(7,7)=f_cojinete/(c_cojinete*omega_cojinete)*(cxx);
c_cojinetes(7,8)=f_cojinete/(c_cojinete*omega_cojinete)*(cxy);
c_cojinetes(8,7)=f_cojinete/(c_cojinete*omega_cojinete)*(cyx);
c_cojinetes(8,8)=f_cojinete/(c_cojinete*omega_cojinete)*cyy;
               
%Definiendo la matriz masa de aporte del disco               
% md: masa del disco detemrinada por la geometria y el material
% jdd: Inercia de la masa del disco

jdd=1;

%espesor del disco
e=0.00555;
%diametro del disco
Dd=0.1245;

md=pi*densidad*e*(Dd/4)^2;
%Diametro del eje


%Inercia polar del disco
ipd=(md/2)*(Dd^2-de^2)/4;

%caclulando la inercia de masa del disco 
jdd=0.5*ipd+(md*e^2)/12;

%agregando los valores de la matriz de masa del disco
m_discos(3,3)=md;
m_discos(4,4)=md;
m_discos(5,5)=jdd;
m_discos(6,6)=jdd;  

g_discos(5,6)=omega*ipd;
g_discos(6,5)=-omega*ipd;
  


%Definiendo la matriz de masa                      
M=Me+Mir;

%Definiendo la matriz de rigidez
K=Ke+beta*Kai*omega;

%Matriz de efecto giroscopio
Ce=alfa*M+beta*K;

%Matriz de amortiguamiento
C=beta*Ke+Cg*omega;


%Agregando los componentes del cojinete a la matriz de rigidez y de amortiguamiento
K=K+k_cojinetes+k_soportes;
C=C+c_cojinete;
M=M+m_discos+g_discos;

%Creacion de las matrices de masa global y la matriz de amortiguamiento y de rigidez
%Aqui definimos la matriz de masa y de amortiguamiento del sistema global.
M_Global=zeros(4+4*n,4+4*n);
C_Global=zeros(4+4*n,4+4*n);
K_Global=zeros(4+4*n,4+4*n);
%Las siguientes instrucciones son usadas para superponer la marices de masa y de rigidez en una matriz global del sistema
for w =0:4:4*n-4
  for(i=1:8)
    for(j=1:8)
        if((w+i)<=(n*4+4) && (w+j)<=(n*4+4))
          M_Global(w+i,w+j) = M_Global(w+i,w+j)+M(i,j);  % Superponer la matriz actual a la matriz superpuesta
          C_Global(w+i,w+j) = C_Global(w+i,w+j)+C(i,j);  % Superponer la matriz actual a la matriz superpuesta
          K_Global(w+i,w+j) = K_Global(w+i,w+j)+K(i,j);  % Superponer la matriz actual a la matriz superpuesta
        end
    end
  end
end  


%Para encontrar los modos normales de vibracion definimos el problema de autovalores y autovectores
Zeros=zeros(4+4*n,4+4*n);
Identidad=eye(4+4*n,4+4*n);


M_prima=-inv(M_Global)*K_Global;
C_prima=-inv(M_Global)*C_Global;

%La matriz A del cual queremos encontrar la solucion es la siguiente
%Contiene 4 elementos;
%Una matriz de Zeros
%Una matriz identidad
%Una matriz de rigidez K expresacda como el negativo de la  inversa de la matriz M_Global multiplicada por 
%Una matriz C expresada como el negativo de la inversa de la matriz global multiplicada por la C_Global
A=[Zeros, Identidad;
  M_prima, C_prima];   

%se utiliza la rutina de Octave para calcular los autovalores y autovectores.
[P,D]=eig(A);

%mostramos los modos normales de vibracion 
for k=1:(n*8)
 D(k,k)
end 
