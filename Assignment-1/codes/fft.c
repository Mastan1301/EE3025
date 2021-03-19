# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <complex.h>
# include <math.h>
typedef double complex cd;

double PI = acos(-1);

void fft(cd** x, int n){
    if (n == 1)
        return;

    cd* odd = malloc(n/2 * sizeof(cd)), *even = malloc(n/2 * sizeof(cd));
    for (int i = 0; 2 * i < n; i++) {
        even[i] = (*x)[2 * i];
        odd[i] = (*x)[2 * i + 1];
    }

    fft(&even, n/2);
    fft(&odd, n/2);

    double ang = 2 * PI / n;
    cd w = CMPLX(1, 0), wn = CMPLX(cos(ang), sin(ang));
    for (int i = 0; 2 * i < n; i++) {
        (*x)[i] = even[i] + w * odd[i];
        (*x)[i + n/2] = even[i] - w * odd[i];
        w *= wn;
    }
    free(even);
    free(odd);
}

void ifft(cd** x, int n){
    if (n == 1)
        return;

    cd* odd = malloc(n/2 * sizeof(cd)), *even = malloc(n/2 * sizeof(cd));
    for (int i = 0; 2 * i < n; i++) {
        even[i] = (*x)[2 * i];
        odd[i] = (*x)[2 * i + 1];
    }

    ifft(&even, n/2);
    ifft(&odd, n/2);

    double ang =  -2 * PI / n;
    cd w = CMPLX(1, 0), wn = CMPLX(cos(ang), sin(ang));
    for (int i = 0; 2 * i < n; i++) {
        (*x)[i] = (even[i] + w * odd[i])/2;
        (*x)[i + n/2] = (even[i] - w * odd[i])/2;
        w *= wn;
    }
    free(even);
    free(odd);
}

int main(){
    int n = (1 << 10);
    cd* x = malloc(n * sizeof(cd));
    FILE *fin, *fout;

    fin = fopen("./x.dat", "r");
    // Generating random data and storing in a dat file
    for(int i = 0; i < n; i++){
        double val = (double)rand()/RAND_MAX * 2.0 - 1.0; // float in range -1 to 1
        x[i] = CMPLX(val, 0);
        fprintf(fin, "%f\n", creal(val));
        // printf("%f\n", creal(x[i]));
    }

    cd* X = malloc(n * sizeof(cd));
    cd* x1 = malloc(n * sizeof(cd));
    for(int i = 0; i < n; i++){
        X[i] = x[i];
    }

    fft(&X, n);

    for(int i = 0; i < n; i++){
        x1[i] = X[i];
    }

    ifft(&x1, n);

    double error = 0, norm = 0;
    for(int i = 0; i < n; i++){
        error += creal((x[i] - x1[i]) * (x[i] - x1[i]));
        norm += creal(x[i] * x[i]);
    }
    error = sqrt(error);
    norm = sqrt(norm);
    error /= norm;
    printf("Relative error: %f\n", error);

    fout  = fopen("./X.dat", "w");
    for(int i = 0; i < n; i++){
        fprintf(fout, "%f + %fi\n", creal(X[i]), cimag(X[i]));
    }
    fclose(fout);
    fclose(fin);

    return 0;
}