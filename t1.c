#include "header.h"

#define MAX_GRAU_FUNCAO 5 //se mudar, alterar tambem a quantidade de parâmetros de iniciaFuncao()
#define MAX_ERRO (int)0x3F3F3F3F
#define QTE_FUNCOES 4

typedef enum {
   polinomial,
   senoide
} tipoFuncao;

typedef struct {
   double coef[MAX_GRAU_FUNCAO + 1]; //coeficiente do termo de grau i, variando de 0 a MAX_GRAU_FUNCAO
   tipoFuncao tipo;
} funcao;

double img(funcao f, double x) {
   double y;
   y = f.coef[0];
   for (int i = 1; i <= MAX_GRAU_FUNCAO; i++) {
      if (f.tipo == polinomial)
         y += f.coef[i] * pow(x, i);
      else if (f.tipo == senoide)
         y += f.coef[i] * pow(sin(x), i);
   }
   return y;
}
double df(funcao f, double x, double e, int maxIte) {
   double erroAtual = MAX_ERRO, erroAnt, h = 1;
   double dfAtual, dfAnt; //derivada no ponto x

   dfAtual = (img(f, x + h) - img(f, x - h)) / (2 * h); //1a iteracao
   for (int i = 2; i <= maxIte; i++) {
      h /= 2;
      dfAnt = dfAtual;
      dfAtual = (img(f, x + h) - img(f, x - h)) / (2 * h);

      erroAnt = erroAtual;
      erroAtual = erro(dfAtual, dfAnt);

      if (erroAtual < e || erroAtual > erroAnt)
         return dfAtual;
   }
   return dfAtual;
}
double df2(funcao f, double x, double e, int maxIte) {
   double erroAtual = MAX_ERRO, erroAnt, h = 1;
   double dfAtual, dfAnt; //derivada no ponto x

   dfAtual = (img(f, x + h) - img(f, x - h)) / (2 * h); //1a iteracao
   for (int i = 2; i <= maxIte; i++) {
      h /= 2;
      dfAnt = dfAtual;
      dfAtual = (img(f, x + 2 * h) - 2 * img(f, x) + img(f, x - 2 * h)) / pow(2 * h, 2);

      erroAnt = erroAtual;
      erroAtual = erro(dfAtual, dfAnt);

      if (erroAtual < e || erroAtual > erroAnt)
         return dfAtual;
   }
   return dfAtual;
}
bool bissecao(funcao f, double a, double b, double e, int maxIte, double *x, int *ite) {
   if (img(f, a) * img(f, b) > 0)
      return false; //intervalo inicial invalido

   for (*ite = 1; *ite <= maxIte; (*ite)++) {
      //media aritmetica entre os limites do intervalo
      *x = (a + b) / 2;

      if (img(f, a) * img(f, *x) < 0)
         b = *x;
      else
         a = *x;

      if (abs(img(f, *x)) < e || abs(b - a) < e)
         return true;
   }
   (*ite)--; //estourou numero maximo de iterações
   return true;
}
bool posicaoFalsa(funcao f, double a, double b, double e, int maxIte, double *x, int *ite) {
   if (img(f, a) * img(f, b) > 0) //intervalo inicial invalido
      return false;

   for (*ite = 1; *ite <= maxIte; (*ite)++) {
      double moduloImgB = abs(img(f, b));
      double moduloImgA = abs(img(f, a));
      //media ponderada entre os limites do intervalo
      *x = (a * moduloImgB + b * moduloImgA) / (moduloImgB + moduloImgA);

      if (img(f, a) * img(f, *x) < 0)
         b = *x;
      else
         a = *x;

      if (abs(img(f, *x)) < e || abs(b - a) < e)
         return true;
   }
   (*ite)--; //estourou numero maximo de iterações
   return true;
}
bool newton(funcao f, double x_ant, double e, int maxIte, double *x, int *ite) {
   double df_1a; //1a derivada no ponto x_ant
   for (*ite = 1; *ite <= maxIte; (*ite)++) {
      if ((df_1a = df(f, x_ant, e, maxIte)) == 0) //nao converge
         return false;

      *x = x_ant - img(f, x_ant) / df_1a;

      if (abs(img(f, *x)) < e || erro(*x, x_ant) < e)
         return true;
      x_ant = *x;
   }
   (*ite)--; //estourou numero maximo de iterações
   return true;
}

void iniciaFuncao(funcao *f, int c0, int c1, int c2, int c3, int c4, int c5, tipoFuncao tipo) {
   (*f).tipo = tipo;
   (*f).coef[0] = c0;
   (*f).coef[1] = c1;
   (*f).coef[2] = c2;
   (*f).coef[3] = c3;
   (*f).coef[4] = c4;
   (*f).coef[5] = c5;
}
void mostraFuncao(funcao f, int op) {
   printf("%d › ", op);
   for (int j = MAX_GRAU_FUNCAO; j >= 1; j--)
      if (f.coef[j] != 0) {
         if (f.tipo == polinomial)
            printf("(%.2f).x^%d + ", f.coef[j], j);
         else if (f.tipo == senoide)
            printf("(%.2f).sen(x)^%d + ", f.coef[j], j);
      }
   printf("(%.2f)\n", f.coef[0]);
}
void leDouble(char *str, double *var) {
   printf("%s: ", str);
   fflush(stdin);
   scanf("%lf", var);
}
void leInt(char *str, int *var) {
   printf("%s: ", str);
   fflush(stdin);
   scanf("%d", var);
}

int main() {
   int opF, opM, ite, maxIte;
   double a, b, e, x, x_inicial;
   char r;
   funcao f[QTE_FUNCOES];
   iniciaFuncao(&f[0], 6, 1, -4, 1, 0, 0, polinomial);
   iniciaFuncao(&f[1], -10, -7, 4, 1, 0, 0, polinomial);
   iniciaFuncao(&f[2], -3, -1, 0, 0, 3, 0, polinomial);
   iniciaFuncao(&f[3], 0, 1, 0, 0, 0, 0, senoide);
   do {
      system("clear");
      printf("   Cálculo de zero de funções   \n");
      printf("    e diferenciação numérica    \n");
      printf("................................\n");
      printf("\nQual a função?\n");
      for (int i = 0; i < QTE_FUNCOES; i++)
         mostraFuncao(f[i], i);
      {
         leInt("Resp.", &opF);
      }
      while (opF < 0 || opF >= QTE_FUNCOES)
         ;

      printf("\nQual o método?\n");
      printf("0 › Bisseção\n");
      printf("1 › Posição Falsa\n");
      printf("2 › Newton\n");
      printf("3 › Derivada primeira no ponto\n");
      printf("4 › Derivada segunda no ponto\n");
      do {
         leInt("Resp.", &opM);
      } while (opM < 0 || opM > 4);

      printf("\n");
      leDouble("Defina um erro", &e);
      leInt("Defina um número máximo de iterações", &maxIte);
      if (opM == 0 || opM == 1) {
         leDouble("Defina um limite inferior", &a);
         leDouble("Defina um limite superior", &b);
         if (opM == 0 && bissecao(f[opF], a, b, e, maxIte, &x, &ite))
            printf("A solução é %.4lf em %d iterações\n", x, ite);
         else if (opM == 1 && posicaoFalsa(f[opF], a, b, e, maxIte, &x, &ite))
            printf("A solução é %.4lf em %d iterações\n", x, ite);
         else
            printf("\nNão foi possivel aplicar o método!\n");
      } else if (opM == 2) {
         leDouble("Defina um ponto inicial", &x_inicial);
         if (newton(f[opF], x_inicial, e, maxIte, &x, &ite))
            printf("A solução é %.4lf em %d iterações\n", x, ite);
         else
            printf("\nNão foi possivel aplicar o método!\n");
      } else if (opM == 3 || opM == 4) {
         leDouble("Defina um ponto", &x);
         if (opM == 3)
            printf("A solução é %.4lf\n", df(f[opF], x, e, maxIte));
         else if (opM == 4)
            printf("A solução é %.4lf\n", df2(f[opF], x, e, maxIte));
      }
      do {
         printf("\nRetornar ao menu [s/n]? ");
         fflush(stdin);
         scanf(" %c", &r);
      } while (r != 'n' && r != 's');
   } while (r == 's');
   return 0;
}