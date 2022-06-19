#include "fft.cpp"
using namespace std;
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <numeric>
#include <algorithm>

void pmult(complex<double>* res, complex<double>* p1, complex<double>* p2, int n){
    for (int i = 0; i < n; i ++){
        res[i] = p1[i] * p2[i];
    }
}

void polynomial_multiplication(complex<double>* res, complex<double>* p1, complex<double>* p2, int n){
    complex<double> p1_fft[2 * n], p2_fft[2 * n], fft_res[2 * n];
    fft(p1_fft, p1, 2 * n);
    fft(p2_fft, p2, 2 * n);
    pmult(fft_res, p1_fft, p2_fft, 2 * n);
    ifft(res, fft_res, 2 * n);
};

int main(){
    const int n = 4;
    complex<double> a[2 * n] = {1,2,3,4}, b[2 * n] = {5,6,7,8}, res[2 * n];
    polynomial_multiplication(res, a, b, 4);
    for (int i = 0; i < 2 * n; i ++){
        cout << res[i].real() << ", ";
    }
    // we input (1 + 2x + 3x^2 + 4x^3) * (5 + 6x + 7x^2 + 8x^3)
    // as we can see this algorithm outputs: 5 + 16x + 34x^2 + 60x^3 + 61x^4 + 52x^5 + 32x^6
    // we can easily verify this to be the right answer
    cout << endl;

}