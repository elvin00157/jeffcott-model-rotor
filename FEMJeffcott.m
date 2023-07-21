function x = FEMPJeffcott
% Encuentra la solucion de un modelo de Jeffcott modelado con elemento finito
clc;
clear;

%Definiendo la cantidad de elementos que va tener la malla
n=2;

%La lonjitud del eje
L=0.857;

%El modulo de Young
E=210000000000;

%El diametro del eje
de=0.0127;

%La densidad del eje
densidad=7850;

%El area transversal del eje
Ae=pi*de^2/4;

%Los valores de alfa y beta se otienen de manera experimental
alfa=0;
beta=0.0000135;

%La velocidad angular del eje
omega=1200*(2*pi)*(1/60);

%Defimi,os la inercia del eje
I=(pi/4)*(de/2)^4;

%Definimos la matriz de rigidez elastica
Ke=(E*I/L^3)*[12, 0, 0, 6*L, -12, 0, 0, 6*L;
              0, 12,-6*L, 0, 0, -12, -6*L, 0;
              0, -6*L, 4*L^2, 0, 0, 6*L, 2*L^2, 0;
              6*L, 0, 0, 4*L^2, -6*L, 0, 0, 2*L^2,
              -12, 0, 0, -6*L, 12, 0, 0, -6*L;
              0, -12, 6*L, 0, 0, 12, 6*L, 0;
              0, -6*L, 2*L^2, 0, 0, 6*L, 4*L^2, 0;
              6*L, 0, 0, 2*L^2, -6*L, 0,0,4*L^2];
   

%Definimos la matriz de rigidez de amrtiguamiento interno.
Kai=(E*I/L^3)*[0, 12, -6*L, 0, 0, -12, -6*L, 0;
              -12, 0, 0, -6*L, 12, 0, 0, -6*L;
              6*L, 0 , 0, 4*L^2, -6*L, 0,0,2*L^2;
              0, 6*L, -4*L^2, 0, 0, -6*L, -2*L^2, 0;
              0, -12, 6*L, 0, 0, 12, 6*L, 0;
              12, 0, 0, 6*L, -12, 0, 0, 6*L;
              6*L, 0, 0, 2*L^2,-6*L, 0, 0, 4*L^2;
              0, 6*L, -2*L^2, 0, 0, -6*L, -4*L^2,0];

%Definimos la matriz de masa              
Me=(densidad*Ae*L/840)*[312, 0, 0, 44*L, 108, 0, 0, -26*L;
                        0, 312, -44*L, 0, 0, 108, 26*L, 0;
                        0, -44*L, 8*L^2, 0, 0, -26*L, -6*L^2, 0;
                        44*L, 0, 0, 8*L^2, 26*L, 0, 0, -6*L^2;
                        108, 0, 0, 26*L, 312, 0, 0,-44*L;
                        0, 108, -26*L, 0, 0, 312, 44*L, 0;
                        0, 26*L, -6*L^2, 0, 0, 44*L, 8*L^2,0;
                        -26*L, 0, 0, -6*L^2, -44*L, 0, 0, 8*L^2];

%Definimos la matriz de inercia rotacional
Mir=(densidad*I/(30*L))*[36, 0 ,0, 3*L, -36, 0, 0, 3*L;
                        0, 36, -3*L, 0, 0, -36, -3*L, 0;
                        0, -3*L, 4*L^2, 0, 0, 3*L, -L^2, 0;
                        3*L, 0, 0, 4*L^2, -3*L, 0, 0, -L^2;
                        -36, 0, 0, -3*L,36, 0, 0, -3*L^2;
                        0, -36, 3*L, 0, 0, 36, 3*L, 0;
                        0, -3*L, -L^2, 0, 0, 3*L, 4*L^2, 0;
                        3*L, 0, 0, -L^2, -3*L, 0, 0, 4*L^2];
 
%Definimos la matriz de efecto giroscopio 
Cg=(densidad*I/15*L)*[0, -36, 3*L, 0, 0, 36, 3*L, 0;
                      36, 0, 0, 3*L, -36, 0, 0, 3*L;
                      -3*L, 0, 0, -4*L^2, 3*L, 0, 0, L^2;
                      0, -3*L, 4*L^2, 0, 0, 3*L, -L^2, 0;
                      0, 36, -3*L, 0, 0, -36, -32, 0;
                      -36, 0, 0, -3*L, 36, 0, 0 -3*L;
                      -3*L, 0, 0, L^2, 3*L, 0, 0, -4*L^2;
                      0, -3*L, -L^2, 0, 0, 3*L, -4*L^2,0];


%definimos el epsioln que es el tamanio del elemento de la malla.
epsilon=L/(n-1);

%Definimos la fuerza sobre el cojinete.
f_cojinete=30;

%Definimos la constante de amortiguamiento del cojinete
c_cojinete=0.001;

%definimos la velocidad angular del cojinete.
omega_cojinete=omega;

h0_cojinete=1/(pi^2*(1-epsilon^2)+16*epsilon^2)^(3/2);

%Definiendo las componentes en x e y del componente de rigidez
kxx_cojinete=4*h0_cojinete*(pi^2*(2-epsilon^2)+16*epsilon^2);
kxy_cojinete=(h0_cojinete*pi*pi^2*(1-epsilon^2)^2+16*epsilon^4)/(epsilon*sqrt(1-epsilon^2));
kyx_cojinete=-1*(h0_cojinete*pi*pi^2*(1-epsilon^2)*(1+2*epsilon^2)+32*epsilon^2*(1+epsilon^2))/(epsilon*sqrt(1-epsilon^2));
kyy_cojinete=4*h0_cojinete*(pi^2*(1+2*epsilon^2)+(32*epsilon^2*(1+epsilon)^2)/(1-epsilon^2));

%definiendo las componentes x e y del amortiguamiento
cxx=h0_cojinete*(2*pi*(1-epsilon^2)^(1/2)*(pi^2)*(1+2*epsilon^2)-16*epsilon^2)/epsilon;
cxy=-8*h0_cojinete*pi^2*(pi^2*(1+2*epsilon^2)+16*epsilon^2);
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

k_soportes(1,1)=(f_cojinete/c_cojinete)*(kxx_cojinete);
k_soportes(2,2)=(f_cojinete/c_cojinete)*(kyy_cojinete);
      
%Definiendo la matriz de amortiguamiento del cojinete.     
      
c_cojinetes(1,1)=f_cojinete/(c_cojinete*omega_cojinete)*(cxx);
c_cojinetes(1,2)=f_cojinete/(c_cojinete*omega_cojinete)*(cxy);
c_cojinetes(2,1)=f_cojinete/(c_cojinete*omega_cojinete)*(cyx);
c_cojinetes(2,2)=f_cojinete/(c_cojinete*omega_cojinete)*cyy;  
               
%Definiendo la matriz masa de aporte del disco               
% md: masa del disco detemrinada por la geometria y el material
% jdd: Inercia de la masa del disco
md=1;
jdd=1;

%espesor del disco
e=0.5;
%diametro del disco
Dd=1;

%Diametro del eje
De=0.5;

%Inercia polar del disco
ipd=(md/2)*(Dd^2-De^2)/4;

%caclulando la inercia de masa del disco 
jdd=0.5*ipd+(md*e^2)/12;

%agregando los valores de la matriz de masa del disco
m_discos(1,1)=md;
m_discos(2,2)=md;
m_discos(3,3)=jdd;
m_discos(4,4)=jdd;  

g_discos(3,4)=omega*ipd;
g_discos(4,3)=-omega*ipd;
  


%Definiendo la matriz de masa                      
M=Me+Mir;

%Definiendo la matriz de rigidez
K=Ke+beta*Kai*omega;

%Matriz de efecto giroscopio
Cg=alfa*M+beta*K;

%Matriz de amortiguamiento
C=beta*Ke+Cg*omega+c_cojinete;

%Agregando los componentes del cojinete a la matriz de rigidez y de amortiguamiento
K=K+k_cojinetes+k_soportes;
C=C+c_cojinete;
M=m_discos+g_discos;

%Creacion de las matrices de masa global y la matriz de amortiguamiento y de rigidez
%Aqui definimos la matriz de masa y de amortiguamiento del sistema global.
M_Global=zeros(8*n,8*n);
C_Global=zeros(8*n,8*n);
K_Global=zeros(8*n,8*n);

%Las siguientes instrucciones son usadas para superponer la marices de masa y de rigidez en una matriz global del sistema
for w =0:8*(n)
  for(i=1:8)
    for(j=1:8)
        if((w+i)<=(n*8) && (w+j)<=(n*8))
          M_Global(w+i,w+j) = M_Global(w+i,w+j)+M(i,j);  % Superponer la matriz actual a la matriz superpuesta
          C_Global(w+i,w+j) = C_Global(w+i,w+j)+C(i,j);  % Superponer la matriz actual a la matriz superpuesta
          K_Global(w+i,w+j) = K_Global(w+i,w+j)+K(i,j);  % Superponer la matriz actual a la matriz superpuesta
        end
    end
  end
end  

%Para encontrar los modos normales de vibracion definimos el problema de autovalores y autovectores
Zeros=zeros(8*n,8*n);
Identidad=eye(8*n,8*n);

%Agregando los componentes del cojinete a la matriz de rigidez y de amortiguamiento
%K_Global=K_Global+k_cojinetes+k_soportes;
%C_Global=C_Global+c_cojinete;

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
real(D)