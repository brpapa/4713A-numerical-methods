#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 100
#define abs(n) ((n) > 0 ? (n) : -(n))
//erro absoluto e relativo
#define erro(vExato, vAprox) ((abs(vExato) > 1) ? (abs(vExato - vAprox) / abs(vExato)) : (abs(vExato - vAprox)))

typedef enum {
   false, //valor 0 por padrão
   true   //valor 1 por padrão
} bool;

void impVetor(int n, double b[]) {
   for (int i = 0; i < n; i++)
      printf("%10.4lf\n", b[i]);
   printf("\n");
}
void impMatriz(int n, double a[][MAX]) {
   for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++)
         printf("%10.4lf ", a[i][j]);
      printf("\n");
   }
   printf("\n");
}

void leVetor(int n, double v[]) {
   for (int i = 0; i < n; i++)
      scanf("%lf", &v[i]);
}
void leMatriz(int n, double a[][MAX]) {
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         scanf("%lf", &a[i][j]);
}