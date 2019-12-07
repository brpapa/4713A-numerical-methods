//integração numérica e resolução de sistemas não lineares
//FALTA SISTEMA NÃO LINEAR E O MENU

#include "main.h"
#define MAX_GRAU_FUNCAO 3

//n: qte de subintervalos
//a: limite inferior
//b: limite superior

typedef struct {
   double coef[MAX_GRAU_FUNCAO+1]; //coeficiente do termo de grau i, variando de 0 a MAX_GRAU_FUNCAO
} funcao;

void iniciaFuncao(funcao *f, int c0, int c1, int c2, int c3) {
   (*f).coef[0] = c0;
   (*f).coef[1] = c1;
   (*f).coef[2] = c2;
   (*f).coef[3] = c3;
}
double img(funcao f, double x) {
   double y = 0;
   for (int i = 0; i <= MAX_GRAU_FUNCAO; i++)
      y += f.coef[i] * pow(x, i);
   return y;
}

double trapezio_1ap(double y0, double y1, double h) {
   return (h * 1/2) * (y0 + y1);
}
double umTercoSimpson_1ap(double y0, double y1, double y2, double h) {
   return (h * 1/3) * (y0 + 4*y1 + y2);
}
double tresOitavosSimpson_1ap(double y0, double y1, double y2, double y3, double h) {
   return (h * 3/8) * (y0 + 3*y1 + 3*y2 + y3);
}

//calcula o vetor y = {f(x), f(x+h), f(x+2h), ... , f(x+n*h)} de tamanho n+1
void calculaY(int n, double h, double x, funcao f, double y[MAX]) {
   for (int i = 0; i <= n; i++)
      y[i] = img(f, x + i*h);
}

double trapezio(int n, double a, double b, funcao f) {
   double res = 0, y[MAX];
   double h = (b - a) / n;

   calculaY(n, h, a, f, y);

   //n aplicações
   for (int i = 0; i + 1 <= n; i++)
      res += trapezio_1ap(y[i], y[i+1], h);

   return res;
}

double umTercoSimpson(int n, double a, double b, funcao f) {
   double res = 0, y[MAX];
   double h = (b - a) / n;

   calculaY(n, h, a, f, y);

   // n/2 aplicações, logo n precisa ser múltiplo de 2
   for (int i = 0; i + 2 <= n; i += 2)
      res += umTercoSimpson_1ap(y[i], y[i+1], y[i+2], h);

   if (n%2 != 0) 
      res += trapezio_1ap(y[n-1], y[n], h); // faltou 1 intervalo

   return res;
}

double tresOitavosSimpson(int n, double a, double b, funcao f) {
   double res = 0, y[MAX];
   double h = (b - a) / n;

   calculaY(n, h, a, f, y);

   // n/3 aplicações, logo n precisa ser múltiplo de 3
   for (int i = 0; i + 3 <= n; i += 3)
      res += tresOitavosSimpson_1ap(y[i], y[i+1], y[i+2], y[i+3], h);

   if (n%3 != 0) {
      if (n%2 == 0) 
         res += umTercoSimpson_1ap(y[n-2], y[n-1], y[n], h); // faltaram 2 intervalos
      else res += trapezio_1ap(y[n-1], y[n], h); // faltou 1 intervalo
   }

   return res;
}

int main() {
   funcao f1, f2, f3;
   iniciaFuncao(&f1, 0, 0, 1, 0);
   iniciaFuncao(&f2, 1, 0, 0, 1);

   printf("%lf\n", tresOitavosSimpson(7, -1, 1, f2));

   return 0;
}