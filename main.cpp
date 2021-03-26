#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double EPS = 0.000001;

void gauss(int n, double **sourse_matr, double *x, int mod)//sourse_matr[n][n + 1], x[n]
{
    double matr[n][n + 1];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n + 1; j++)
            matr[i][j] = sourse_matr[i][j];
    if(mod == 1){
        for (int k = 0; k < n; k++) {
            for (int i = k; i < n; i++) {
                double tmp = matr[i][k];
                if (fabs(tmp) < EPS)
                    continue;
                for (int j = 0; j < n + 1; j++)
                    matr[i][j] = matr[i][j] / tmp;
                if (i == k)
                    continue;
                for (int j = 0; j < n + 1; j++)
                    matr[i][j] = matr[i][j] - matr[k][j];
            }
        }
    } else {
        int max_index;
        double max_elem;
        for (int k = 0; k < n; k++) {
            max_elem = fabs(matr[k][k]);
            max_index = k;
            for (int i = k + 1; i < n; i++)
                if (max_elem < fabs(matr[i][k])) {
                    max_elem = fabs(matr[i][k]);
                    max_index = i;
                }
            if (max_elem < EPS)
                continue;
            for (int j = 0; j < n + 1; j++) {
                double tmp = matr[k][j];
                matr[k][j] = matr[max_index][j];
                matr[max_index][j] = tmp;
            }
            for (int i = k; i < n; i++) {
                double tmp = matr[i][k];
                if (fabs(tmp) < EPS)
                    continue;
                for (int j = 0; j < n + 1; j++)
                    matr[i][j] = matr[i][j] / tmp;
                if (i == k)
                    continue;
                for (int j = 0; j < n + 1; j++)
                    matr[i][j] = matr[i][j] - matr[k][j];
            }
        }
    }
    for (int k = n - 1; k >= 0; k--)
    {
        x[k] = matr[k][n];
        for (int i = 0; i < k; i++)
            matr[i][n] = matr[i][n] - matr[i][k] * x[k];
    }
}


double deter(int n, double **sourse_matr)
{
    double matr[n][n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matr[i][j] = sourse_matr[i][j];
    double det = 1;
    for (int k = 0; k < n; k++) {
        for (int i = k; i < n; i++) {
            double tmp = matr[i][k];
            if (fabs(tmp) < EPS)
                continue;
            for (int j = 0; j < n; j++)
                matr[i][j] = matr[i][j] / tmp;
            det *= tmp;
            if (i == k)
                continue;
            for (int j = 0; j < n; j++)
            matr[i][j] = matr[i][j] - matr[k][j];
        }
    }
    for (int i = 0; i < n; i++)
        det *= matr[i][i];
    return det;
}

void inverse(int n, double **sourse_matr, double **invrerse_matr)//invrerse_matr[n][n]
{
    double matr[n][n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matr[i][j] = sourse_matr[i][j];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i == j)
                invrerse_matr[i][j] = 1;
		 else
		     invrerse_matr[i][j] = 0;
    for (int k = 0; k < n; k++) {
        for (int i = k; i < n; i++) {
            double tmp = matr[i][k];
            if (fabs(tmp) < EPS)
                continue;
            for (int j = 0; j < n; j++) {
                matr[i][j] = matr[i][j] / tmp;
                invrerse_matr[i][j] = invrerse_matr[i][j] / tmp;
            }
            if (i == k)
            continue;
            for (int j = 0; j < n; j++) {
                matr[i][j] = matr[i][j] - matr[k][j];
                invrerse_matr[i][j] = invrerse_matr[i][j] - invrerse_matr[k][j];
            }
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        for (int i = n - 1; i >= 0; i--) {
            double tmp = matr[i][k];
            if (fabs(tmp) < EPS)
                continue;
            for (int j = 0; j < n; j++) {
                matr[i][j] = matr[i][j] / tmp;
                invrerse_matr[i][j] = invrerse_matr[i][j] / tmp;
            }
            if (i == k)
            continue;
            for (int j = 0; j < n; j++) {
                matr[i][j] = matr[i][j] - matr[k][j];
                invrerse_matr[i][j] = invrerse_matr[i][j] - invrerse_matr[k][j];
            }
        }
    }
}

double cond(int n, double **sourse_matr)
{
    double cond1 = 0;
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++)
            sum += fabs(sourse_matr[i][j]);
        if ((i == 0) || sum > cond1)
            cond1 = sum;
    }

    double **invrerse_matr = (double **)calloc(n, sizeof(*invrerse_matr));
    for (int i = 0; i < n; i++) {
        invrerse_matr[i] = (double *)calloc(n, sizeof(**invrerse_matr));
    }
    inverse(n, sourse_matr, invrerse_matr);
    double cond2 = 0;
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++)
            sum += fabs(invrerse_matr[i][j]);
        if ((i == 0) || sum > cond2)
            cond2 = sum;
    }
    free(invrerse_matr);
    return cond2 * cond1;
}


int main(int argc, const char *argv[]) {
    printf("Size of matr:\n");
    int n;
    scanf("%d", &n);
    double **matr = (double **)calloc(n, sizeof(*matr));
    for (int i = 0; i < n; i++) {
        matr[i] = (double *)calloc(n + 1, sizeof(**matr));
    }

    printf("Choose the option(1 or 2):\n");
    int opt = 0;
    scanf("%d", &opt);
    if (opt == 1) {
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n + 1; j++){
                printf("a[%d][%d] = ", i, j);
                scanf("%lf", &matr[i][j]);
            }
            printf("\n");
        }
    } else if (opt == 2) {
        const double q = 1.001 - 2 * 6 * 0.001;
        double x0;
        n = 90;
        free(matr);
        matr = (double **)calloc(n, sizeof(*matr));
        for (int i = 0; i < n; i++) {
            matr[i] = (double *)calloc(n + 1, sizeof(**matr));
        }
        printf("n  is auto-resized to 100\n");
        printf("x:\n");
        scanf("%lf", &x0);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j)
                    matr[i][j] = pow(q, i + 1 + j + 1) + 0.1 * (j - i);
                else
                    matr[i][j] = pow(q - 1, i + 1 + j + 1);
               // printf("%.*lf ", 8, matr[i][j]);
            }
            matr[i][n] = x0 * exp(x0 / (i+1)) * cos(x0 / (i+1));
            //printf("| %.*lf\n", 8, matr[i][n]);
        }
        printf("\n");
    } else {
        printf("It's not correct, try again:\n");
        scanf("%d", &opt);
    }
    double *x1 = (double *)calloc(n, sizeof(*x1));
    double *x2 = (double *)calloc(n, sizeof(*x2));
    double det = deter(n, matr);
    double **invrerse_matr = (double **)calloc(n, sizeof(*invrerse_matr));
    for (int i = 0; i < n; i++) {
        invrerse_matr[i] = (double *)calloc(n, sizeof(**invrerse_matr));
    }
    gauss(n, matr, x1, 1);//Gaussian
    gauss(n, matr, x2, 2);//with pivot selection
    if(det != 0){
        inverse(n, matr, invrerse_matr);
         printf("Inverse matr:\n");
        for (int i = 0; i < n; i++) {
            //for (int j = 0; j < n; j++)
            //    printf("%.*lf ", 5 , invrerse_matr[i][j]);
            //printf("\n");
        }
    }
    printf("Determinant:\ndet A = %f\n", det);
    printf("Gaussian method:\n");
    for (int i = 0; i < n; i++)
        printf("x%d = %.*f ", i + 1, 9, x1[i]);
    printf("\n");
    printf("Gaussian method with pivot selection:\n");
    for (int i = 0; i < n; i++)
        printf("x%d = %.*f ", i + 1, 9, x2[i]);
    printf("\n");
    double sr=0;
    for (int i = 0; i < n; i++)
        sr+=fabs(x2[i]-x1[i]);
    sr/=(double)n;
    printf("\n%f\n",sr);
    for (int i = 0; i < n; i++) {
        free(matr[i]);
        free(invrerse_matr[i]);
    }
    printf("Conditionality number of the matr:\nM = %lf\n", cond(n, matr));
    free(matr);
    free(x1);
    free(x2);
    free(invrerse_matr);
    return 0;
}




