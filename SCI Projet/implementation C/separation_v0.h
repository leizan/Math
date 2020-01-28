#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define S 10

/*********************************************************************************
  *Function: dot_product
  *Description: Calculate inner product of two vectors
  *Calls:  --
  *Called By: main
  *Input: Pointer of two arrays(vec_1, vec_2) and the length of the array(sizeVec)
  *Return:  Inner product(type double)
**********************************************************************************/
double dot_product(double * vec_1,double * vec_2, int sizeVec){
    double res=0;
    for(int i=0; i<sizeVec; ++i){
        res+= vec_1[i]*vec_2[i];
    }
    return res;
}

/*********************************************************************************
  *Function: norm
  *Description: Calculate norm of order 2
  *Calls:  --
  *Called By: main
  *Input: Pointer of one array(vec_1) and the length of the array(sizeVec)
  *Return:  Norm of order 2(type double)
**********************************************************************************/
double norm(double *vec_1, int sizeVec){
    double res=0;
    for(int i=0; i<sizeVec; ++i){
        res+= pow(vec_1[i],2);
    }
    return sqrt(res);
}

/*********************************************************************************
  *Function: off
  *Description: Calculate off of an array
  *Calls:  --
  *Called By: separation_v0.h
  *Input: One array(vec[S][S] S is set to 10) and the size of the array is sizeVec* sizeVec
  *Return:  Off of an array(type double)
**********************************************************************************/
double off(double vec[S][S],int sizeVec){
    double off=0;
    for(int m=0; m<sizeVec; ++m){
        for(int n=0; n<sizeVec; ++n){
            if(m==n){
                    continue;
            }else{
                     off+=vec[m][n];
            }
        }
    }
    off/=(sizeVec*(sizeVec-1));
    return off;
}

/*********************************************************************************
  *Function: tr
  *Description: Calculate trace of an array
  *Calls:  --
  *Called By: separation_v0.h
  *Input: One array(vec[S][S] S is set to 10) and the size of the array is sizeVec* sizeVec
  *Return:  Trace of an array(type double)
**********************************************************************************/
double tr(double vec[S][S] ,int sizeVec){
    double tr=0;
    for(int m=0; m<sizeVec; ++m){
        tr+= vec[m][m];
    }
    tr/=sizeVec;
    return tr;
}


/*********************************************************************************
  *Function: separation_v0
  *Description: Seperate two estimated original signals(output_1, output_2) from two mixed signals(input_1, input_2)
  *Calls:  --
  *Called By: main
  *Input: input_1,input_2 pointer of two mixed signal arrays, output_1, output_2 pointer of two seperated signal array where they are stored, nb_samples: the number of samples of the signal,
  time_delay: the size of the window
  *Return: NULL
**********************************************************************************/
void separation_v0(double *input_1, double *input_2, double *output_1, double *output_2, int nb_samples, int time_delay){
    int N=time_delay;
    int K = nb_samples/time_delay;
    double R11[N][N], R12[N][N], R21[N][N], R22[N][N];
    double T11=0, F11=0, T12=0, F12=0, T22=0, F22=0;
    int conNum = 0;

    //Initialization of correlation matrix
    for(int m=0; m<N; ++m){
            for(int n=0; n<N; ++n){
                R11[m][n]=0;
                R12[m][n]=0;
                R21[m][n]=0;
                R22[m][n]=0;
            }
        }


    //iteration  for calculating auto correlation by going through all time slot
    for(int k= 0; k<K; ++k){
        double x1k[N], x2k[N];
        for(int i=0; i<N; ++i){
            x1k[i] = *(input_1+k*N+i);
            conNum+=1;
            x2k[i] = *(input_2+k*N+i);
        }

        //accumulation correlation matrix of each time slot
        for(int m=0; m<N; ++m){
            for(int n=0; n<N; ++n){
                R11[m][n]+=x1k[m]*x1k[n];
                R12[m][n]+=x1k[m]*x2k[n];
                R21[m][n]+=x2k[m]*x1k[n];
                R22[m][n]+=x2k[m]*x2k[n];
            }
        }

    }
    //get the average of the auto correlation matrix
    for(int m=0; m<N; ++m){
            for(int n=0; n<N; ++n){
                R11[m][n]/=K;
                R12[m][n]/=K;
                R21[m][n]/=K;
                R22[m][n]/=K;
            }
        }
    //calculating Tr and off

    T11=tr(R11, N);
    F11=off(R11, N);
    T12=tr(R12, N);
    F12=off(R12, N);
    T22=tr(R22, N);
    F22=off(R22, N);

    // On the paper
    int variance=0;
    double alpha = 2*F12*T12-(F11*(T22-variance)+F22*(T11-variance));
    double beta=2*(pow(T12,2)-(T11-variance)*(T22-variance));
    double gama=pow((F11*(T22-variance)-F22*(T11-variance)),2)+4*(F12*(T22-variance)-T12*F22)*(F12*(T11-variance)-T12*F11);

    double d1=alpha- sqrt(gama);
    double d2=alpha+ sqrt(gama);


    double A[2][2] = {
        {beta*F11-(T11-variance)*d1, beta*F12-T12*d2} ,
        {beta*F12-T12*d1, beta*F22-(T22-variance)*d2}
    };


   //Don't need normalisation
   /*double col1 = A[0][0]+A[0][1];
   double col2 = A[1][0]+A[1][1];
   for(int u=0; u<2; ++u){
        for(int v=0; v<2; v++){
            if(u==0){
                A[u][v]/= col1;
            }else{
                A[u][v]/= col2;
            }
        }
   }
   */


   double detA =  A[0][0]*A[1][1]-A[1][0]*A[0][1];
   double invA[2][2] = {
        {A[1][1]/detA, -(A[0][1]/detA)},
        {-(A[1][0]/detA), A[0][0]/detA}
   };

   /*printf("%f", A[0][0]*invA[0][0]+A[0][1]*invA[1][0]);
   printf("%f", A[1][0]*invA[0][1]+A[1][1]*invA[1][1]);
   printf("%f", A[0][0]*invA[0][1]+A[0][1]*invA[1][1]);
   printf("%f", A[1][0]*invA[0][0]+A[1][1]*invA[1][0]);*/
   for(int w=0; w< nb_samples; ++w){
        output_1[w] = invA[0][0]*input_1[w] + invA[0][1]*input_2[w];
        output_2[w] = invA[1][0]*input_1[w] + invA[1][1]*input_2[w];
   }

}
