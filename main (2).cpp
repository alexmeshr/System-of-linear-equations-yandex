#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const double EPS = 0.000000001;
void mul (int m, int n, int q, double **a, double **b, double **c)//a[m][n], b[n][q], c[m][q]
{
    for(int i = 0; i < m; i++)
        for(int j = 0; j < q; j++){
            c[i][j] = 0;
            for(int k = 0; k < n; k++)
                    c[i][j] += a[i][k] * b[k][j];

    }

}

void mulvector (int m, int n, double **a, double *b, double *c)
{
for(int i = 0; i < m; i++){
        c[i] = 0;
        for(int k = 0; k < n; k++)
              c[i] += a[i][k] * b[k];
}
}


double norm(int n, double **summ, double *x, double *c)
{
    double sqr_s = 0;
    for (int i = 0; i < n; i++){
        double sum = 0;
        for (int j = 0; j < n; j++){
                sum += summ[i][j] * x[j];
        }
        sqr_s += (sum - c[i]) * (sum - c[i]);
    }
    printf("%lf\n", sqrt(sqr_s));
    return sqrt(sqr_s);
}

long long relaxation(int n, double **summ, double *c, double *x, double w)//summ[n][n], c[n], x[n],
{
    for (int i = 0; i < n; i++){
        x[i] = 0;
    }
    long long cnt = 0;
    double *x_prev = (double *)calloc(n, sizeof(*x_prev));
    do
    {
        cnt++;
        for (int i = 0; i < n; i++)
            x_prev[i] = x[i];
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (summ[i][j] * x[j]);
            for (int j = i; j < n; j++)
                sum += (summ[i][j] * x_prev[j]);

            if(summ[i][i] != 0)
                x[i] = w * (c[i] - sum) / summ[i][i] + x_prev[i];
        }
    } while (norm(n, summ, x, c) > EPS);//printf("| %.*8lf |\n",distance(n, x, x_prev)),
    free(x_prev);
    return cnt;
}

int main(int argc, const char *argv[]) {
    printf("Size of matr:\n");
    printf("Parameter w:\n");
    int n;
    double w;
    scanf("%d%lf", &n, &w);
    double **matr = (double **)calloc(n, sizeof(*matr));
    for (int i = 0; i < n; i++) {
        matr[i] = (double *)calloc(n, sizeof(**matr));
    }
    double **transp = (double **)calloc(n, sizeof(*transp));
    for (int i = 0; i < n; i++) {
        transp[i] = (double *)calloc(n , sizeof(**transp));
    }
    double *b = (double *)calloc(n, sizeof(*b));
    printf("Choose the option(1 or 2):\n");
    int opt = 0;
    scanf("%d", &opt);
    Begin:
    if (opt == 1) {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                printf("a[%d][%d] = ", i, j);
                scanf("%lf", &matr[i][j]);
                transp[j][i] = matr[i][j];
            }
            printf("y[%d] = ", i);
            scanf("%lf", &b[i]);
            printf("\n");
        }
    } else if (opt == 2) {
        const double q = 1.001 - 2 * 6 * 0.001;
        double x0;
        n = 100;
        free(matr);
        free(transp);
        free(b);
        matr = (double **)calloc(n, sizeof(*matr));
        for (int i = 0; i < n; i++) {
            matr[i] = (double *)calloc(n, sizeof(**matr));
        }
        transp = (double **)calloc(n, sizeof(*transp));
        for (int i = 0; i < n; i++) {
            transp[i] = (double *)calloc(n, sizeof(**transp));
        }
        b = (double *)calloc(n, sizeof(*b));
        printf("n  is auto-resized to 100\n");
        printf("x:\n");
        scanf("%lf", &x0);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j)
                    matr[i][j] = pow(q, i + 1 + j + 1) + 0.1 * (j - i);
                else
                    matr[i][j] = pow(q - 1, i + 1 + j + 1);
                transp[j][i] = matr[i][j];
                printf("%.*lf ", 8, matr[i][j]);
            }
            b[i] = x0 * exp(x0 / (i+1)) * cos(x0 / (i+1));
            printf("| %.*lf\n", 8, matr[i][n]);
        }
        printf("\n");

    } else {
        printf("It's not correct, try again:\n");
        scanf("%d", &opt);
        goto Begin;
    }

    double **summ = (double **)calloc(n, sizeof(*summ));
    for (int i = 0; i < n; i++) {
        summ[i] = (double *)calloc(n, sizeof(**summ));
    }
    double *x = (double *)calloc(n, sizeof(*x));
    double *c = (double *)calloc(n, sizeof(*c));
    mul(n, n, n, transp, matr, summ);
    mulvector(n, n, transp, b, c);
    long long cnt = relaxation(n, summ, c, x, w);

    for (int i = 0; i < n; i++) {
        printf("x%d = %.8lf ", i, x[i]);
        free(matr[i]);
        free(summ[i]);
        free(transp[i]);
    }

    printf("\niteration count = %lld\n", cnt);
    free(matr);
    free(x);
    free(summ);
    free(b);
    free(c);
    free(transp);
    return 0;
}




