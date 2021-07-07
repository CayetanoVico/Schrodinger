#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "gsl_rng.h" //Libreria para generación de números aleatorios

#define N 500          //Número de puntos
#define nciclos 20      //A mayor número de ciclos más lento va el gif
#define lambdao 0.1
#define PI 3.141592
#define simulaciones 10000

gsl_rng *tau;

int main()
{
    int i,j,k,l,resto,iteracion;             //Contadores
    double s,ko,Vj[N],aux;
    double intervalo1,intervalo2;
    double mod[N],cnormal,sumpar,sumimpar,integral;
    double P_D,P_I,x;
    double lambda;
    int medicion,mt;
    int nd;
    int interPD1,interPD2,interPI1,interPI2;
    int semilla=679578;       //Semilla de generación de números
    extern gsl_rng *tau;    //Puntero al estado del número aleatorio
    FILE *f1;  
    fcomplex phi[N],alpha[N],beta[N],chi[N];
    fcomplex aux2,aux3,aux4,aux5,aux6,aux7,aux8,aux9;        //Variables auxiliares para facilitar las cuentas    
    tau=gsl_rng_alloc(gsl_rng_taus);    //Inicializamos el puntero
    gsl_rng_set(tau,semilla);           //Inicializamos la semilla
    //Vamos a hacer un for más grande para analizar varias nd en un único programa
    f1=fopen("datos.txt","w");
for(nd=10;nd<151;)
{
    lambda=lambdao;
    fprintf(f1, "%i\tnd:%i\n", semilla,nd);
for(l=0;l<6; l++)
{
    if(lambda>5.1)
    {
        lambda=10;
    }
    else
    {
        if(lambda>1.1)
        {
            lambda=5;
        }
        else
        {
            if(lambda>0.6)
            {
                lambda=1;
            }
        }
    }
    interPD1=0.8*N;
    interPD2=N;
    interPI1=0;
    interPI2=0.2*N;
    mt=0;
    ko=2*PI*nciclos/N;
    s=1/(4.0*ko*ko);
    //Imponemos las condiciones de contorno, phi en los extremos es 0
    phi[0].r=0.0;
    phi[0].i=0.0;
    //Podemos hacerlo de otra forma
    phi[N-1]=Complex(0.0,0.0);
    mod[0]=0.0;
    mod[N-1]=0.0;
    //Sabiendo ko y N podemos calcular el resto de las phi
    for(i=1;i<N-1; i++)
    {
        aux=exp((-8.0*(4.0*i-N)*(4.0*i-N))/(N*N));
        phi[i].r=cos(ko*i)*aux;
        phi[i].i=sin(ko*i)*aux;
        mod[i]=phi[i].r*phi[i].r+phi[i].i*phi[i].i;
    }
    //Calculamos la integral de esta primera phi y la normalizamos.
    sumpar=0.0;
    sumimpar=0.0;
    for(i=1;i<N-1; i++)
    {
        resto=i%2;
        if(resto==0)
        {
            sumpar=sumpar+mod[i];
        }
        else
        {
            sumimpar=sumimpar+mod[i];
        }
    }
    integral=mod[0]+mod[N-1]+4*sumpar+2*sumimpar;
    integral=integral/3.0;
    cnormal=1.0/integral;       //Esta es la cte de normalización de la phi inicial, la cual usaremos para verificar que la norma se conserva
    //Calculamos ahora el potencial
    intervalo1=(2.0*N)/5;   //El extremo inferior
    intervalo2=(3.0*N)/5;   //El extremo superior
    for(i=0;i<N; i++)
    {
        if(i<intervalo1 || i>intervalo2)
        {
            Vj[i]=0.0;
        }
        else
        {
            Vj[i]=lambda*ko*ko;
        }
    }
    //Una vez calculado s y el potencial podemos calcular los alpha, que serán ctes en todo el proceso
    //Esta es una recurrencia inversa, partimos de alpha_N=0 y vamos calculando hacia atrás
    alpha[N-1]=Complex(0.0,0.0);
    aux2=Complex(0.0,0.0);
    aux3=Complex(0.0,0.0);
    aux4=Complex(0.0,0.0);
    aux5=Complex(0.0,0.0);
    aux6=Complex(-1.0,0.0);
    for(i=N-2;i>=0; i--)
    {
        aux2.r=-2-Vj[i+1];
        aux3.i=2/s;
        aux4=Cadd(aux3,alpha[i+1]);
        aux5=Cadd(aux2,aux4);
        alpha[i]=Cdiv(aux6,aux5);
    }
    //Una vez calculadas las alpha ya podemos comenzar a iterar y calcular las
    //Beta, las chi y la nueva phi
    //Vuelvo a inicializar mis variables auxiliares a 0 para volver a utilizarlas y no definir más variables
    aux2=Complex(0.0,0.0);
    aux3=Complex(0.0,0.0);
    aux4=Complex(0.0,0.0);
    aux5=Complex(0.0,0.0);
    aux6=Complex(0.0,0.0);
    aux7=Complex(0.0,0.0);
    aux8=Complex(0.0,0.0);
    beta[N-1]=Complex(0.0,0.0); //Inicializo esta beta y chi fuera del bucle porque siempre van a ser 0
    chi[0]=Complex(0.0,0.0);
    aux3.i=2/s;     //Definimos estas dos variables fuera del bucle porque de la forma en la que estamos realizando las cuentas
    aux6.i=4/s;     //Estas variables siempre van a ser constantes
    i=0;
    iteracion=1;
    for(i=0;i<simulaciones;)
    {
        aux2=Complex(0.0,0.0);
        aux4=Complex(0.0,0.0);
        aux5=Complex(0.0,0.0);
        aux7=Complex(0.0,0.0);
        aux8=Complex(0.0,0.0);
        //Primero vamos a calcular las beta
        for(j=N-2;j>0; j--)
        {
            aux2.r=-2-Vj[j+1];
            aux4=Cadd(aux3,alpha[j+1]);
            aux5=Cadd(aux2,aux4);   //Este es el denominador de la beta
            aux7=Cmul(aux6,phi[j+1]);
            aux8=Csub(aux7,beta[j+1]);  //Este es el numerador
            beta[j]=Cdiv(aux8,aux5);
        }
        //A continuación calculamos las chi. Esta vez no se calcula a la inversa, sino que vamos hacia delante
        for (j=1;j<N; j++)
        {
            aux9=Cmul(alpha[j-1],chi[j-1]);
            chi[j]=Cadd(aux9,beta[j-1]);
        }
        //Finalmente calculamos las nuevas phi
        for(j=0;j<N; j++)
        {
            phi[j]=Csub(chi[j],phi[j]);
        }
        //Calculamos el módulo de phi
        for(j=1;j<N; j++)
        {
            mod[j]=phi[j].r*phi[j].r+phi[j].i*phi[j].i;
        }
        //Calculamos la norma para ver que se conserva
        sumpar=0.0;
        sumimpar=0.0;
        for(j=1;j<N-1; j++)
        {
            resto=j%2;
            if(resto==0)
            {
                sumpar=sumpar+mod[j];
            }
            else
            {
                sumimpar=sumimpar+mod[j];
            }
        }   
        integral=mod[0]+mod[N-1]+4*sumpar+2*sumimpar;
        integral=integral/3.0;
        //A diferencia del problema obligatorio, en este problema calcularemos siempre una nueva constante de 
        //normalizacion para asrgurarnos de que siempre es 1
        cnormal=(1.0)/integral;
        integral=integral*cnormal;
        iteracion=iteracion+1;
        //Hay que poner un filtro porque por como está hecho el programa, puede pasar que la función de onda 
        //entera desaparezca sin haver sido detectada y el programa seguiría corriendo sin parar nunca
        if(iteracion>25000)
        {
            //generamos una nueva función de onda
            for(k=0;k<N; k++)
            {
                aux=exp((-8.0*(4.0*k-N)*(4.0*k-N))/(N*N));
                phi[k].r=cos(ko*k)*aux;
                phi[k].i=sin(ko*k)*aux;
                mod[k]=phi[k].r*phi[k].r+phi[k].i*phi[k].i;
                iteracion=1;
            }
        }
        medicion=iteracion%nd;
        if(medicion==0)
        {
            P_D=0.0;
            //En el pdf pone que hagamos sumatorias como un método aproximado, pero nosotros vamos a hacer
            //la integral para hacerlo más exacto
            sumpar=0.0;
            sumimpar=0.0;
            for(j=interPD1+1;j<interPD2-1; j++)
            {
                resto=j%2;
                if(resto==0)
                {
                    sumpar=sumpar+mod[j];
                }
                else
                {
                    sumimpar=sumimpar+mod[j];
                }
            }
            P_D=mod[interPD1]+mod[interPD2-1]+4*sumpar+2*sumimpar;
            P_D=P_D/3.0;
            P_D=P_D*cnormal;
            x=gsl_rng_uniform(tau);  //número aleatorio real [0,1]
            if(x<P_D)
            {
                mt=mt+1;
                //Realizamos otra simulación
                i=i+1;
                iteracion=1;
                //Y generamos una nueva función de onda
                for(k=0;k<N; k++)
                {
                    aux=exp((-8.0*(4.0*k-N)*(4.0*k-N))/(N*N));
                    phi[k].r=cos(ko*k)*aux;
                    phi[k].i=sin(ko*k)*aux;
                    mod[k]=phi[k].r*phi[k].r+phi[k].i*phi[k].i;
                }
            }
            else
            {
                for(j=interPD1;j<interPD2; j++)
                {
                    phi[j]=Complex(0.0,0.0);
                }
                //Volvemos a calcular la cte de normalizacion para calcular P_I
                for(j=1;j<N; j++)
                {
                    mod[j]=phi[j].r*phi[j].r+phi[j].i*phi[j].i;
                }
                sumpar=0.0;
                sumimpar=0.0;
                for(j=1;j<N-1; j++)
                {
                    resto=j%2;
                    if(resto==0)
                    {
                        sumpar=sumpar+mod[j];
                    }
                    else
                    {
                        sumimpar=sumimpar+mod[j];
                    }
                }   
                integral=mod[0]+mod[N-1]+4*sumpar+2*sumimpar;
                integral=integral/3.0;
                cnormal=1.0/integral;
                //Calculamos P_I
                sumpar=0.0;
                sumimpar=0.0;
                for(j=interPI1+1;j<interPI2-1; j++)
                {
                    resto=j%2;
                    if(resto==0)
                    {
                        sumpar=sumpar+mod[j];
                    }
                    else
                    {
                        sumimpar=sumimpar+mod[j];
                    }
                }
                P_I=mod[interPI1]+mod[interPI2-1]+4*sumpar+2*sumimpar;
                P_I=P_I/3.0;
                P_I=P_I*cnormal;
                x=gsl_rng_uniform(tau);  //número aleatorio real [0,1]
                if(x<P_I)
                {
                    //Realizamos otra simulación
                    i=i+1;
                    iteracion=1;
                    //Y generamos una nueva función de onda
                    for(k=0;k<N; k++)
                    {
                        aux=exp((-8.0*(4.0*k-N)*(4.0*k-N))/(N*N));
                        phi[k].r=cos(ko*k)*aux;
                        phi[k].i=sin(ko*k)*aux;
                        mod[k]=phi[k].r*phi[k].r+phi[k].i*phi[k].i;
                    }
                }
                else
                {
                    for(j=interPI1;j<interPI2; j++)
                    {
                        phi[j]=Complex(0.0,0.0);
                    }
                }
            }
        }       
    }
    fprintf(f1, "%lf\t%i\n", lambda, mt);
    lambda=lambda+0.2;
}
    fprintf(f1,"\n");
    //Cambiamos de semilla
    semilla=gsl_rng_uniform_int(tau,500000);     //Asigno a semilla un valor entero entre [0,N-1]
    semilla=semilla+10000;                      //Sumamos esto para que la semilla no sea muy pequeña
    gsl_rng_set(tau,semilla);           //Inicializamos la semilla otra vez
    nd=nd+5;
}
    fclose(f1);
    return(0);
}
