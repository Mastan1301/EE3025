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
        (*x)[i] = (even[i] + w * odd[i]);
        (*x)[i + n/2] = (even[i] - w * odd[i]);
        w *= wn;
    }
    free(even);
    free(odd);
}

void convolve(cd** Y, cd* X, cd* H, int n){
    for(int i = 0; i < n; i++){
        (*Y)[i] = H[i] * X[i];
    }
}

void shift(cd** X, int n){
    cd* temp = malloc(n * sizeof(cd));
    for(int i = 0; i < n; i++){
        temp[i] = (*X)[i];
    }
    for(int i = 0; i < n; i++){
        (*X)[i] = temp[(i + (n >> 1)) % n];
    }
}

int main(){
    int n = (1 << 20);
    cd* x = malloc(n * sizeof(cd)), *H = malloc(n * sizeof(cd)), *Y = malloc(n * sizeof(cd)), *y = malloc(n * sizeof(cd));
    FILE *fin1, *fout1, *fout2, *fin2;

    fin1 = fopen("../data/x.dat", "r");
    if(!fin1)
    {
        perror("Error opening file");
        return -1;
    }

    int count = 0;
    while (!feof(fin1) && count < n) 
    {
        double val;
        fscanf(fin1, "%lf", &val);
        x[count++] = CMPLX(val, 0);
    }

    cd* X = malloc(n * sizeof(cd));
    for(int i = 0; i < n; i++){
        X[i] = x[i];
    }

    fft(&X, n);

    fin2 = fopen("../data/H.dat", "r");
    if(!fin2)
    {
        perror("Error opening file");
        return -1;
    }

    double r, i;
    count = 0;
    while (!feof(fin2) && count < n) 
    {
        fscanf(fin2, "%lf %lf", &r, &i);
        H[count++] = CMPLX(r, i);
    }

    convolve(&Y, X, H, n);
    
    for(int i = 0; i < n; i++){
        y[i] = Y[i];
    }

    ifft(&y, n);
    
    shift(&Y, n);

    fout1  = fopen("../data/Y.dat", "w");
    if(!fout1)
    {
        perror("Error opening file");
        return -1;
    }
    for(int i = 0; i < n; i++){
        fprintf(fout1, "%lf+%lfi\n", creal(Y[i]), cimag(Y[i]));
    }

    fout2  = fopen("../data/y.dat", "w");
    if(!fout2)
    {
        perror("Error opening file");
        return -1;
    }
    for(int i = 0; i < n; i++){
        y[i] /= n;
        fprintf(fout2, "%lf+%lfi\n", creal(y[i]), cimag(y[i]));
    }

    fclose(fout1);
    fclose(fout2);
    fclose(fin1);
    fclose(fin2);

    return 0;
}