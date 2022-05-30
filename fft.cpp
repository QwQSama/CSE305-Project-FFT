#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <thread>
#define pi  3.14159265358979323846

void fft(std::complex<double>  FFT_p[], std::complex<double>  p[], int n){
    if (n==1){
        for (int i=0; i<n; i++){
            FFT_p[i] = p[i];
        }
        return;
    }
    std::complex<double>  U[n/2],V[n/2],FFT_U[n/2],FFT_V[n/2];
    for (int i=0; i<n/2; i++){
        U[i] = p[2*i];
        V[i] = p[2*i+1];
    }

    fft(FFT_U, U, n/2);
    fft(FFT_V, V, n/2);
    
    std::complex<double>  w,w_n;
    w_n = cos(2*pi/n) + 1j * sin(2*pi/n);
    w = 1;

    for (int i = 0; i<n/2; i++){
        FFT_p[i] = FFT_U[i] + w*FFT_V[i];
        FFT_p[i + n/2] = FFT_U[i] - w*FFT_V[i];
        w = w * w_n;
    }
}

int reverseBits(int n) {
    int res = 0;
    int digit = log2(n);
    for(int i = 0; i < digit; i++){
        res = (res << 1) + (n & 1);
        n = n >> 1;
    }
   return res;
}

void thread_fft(int begin, int end, std::complex<double>  FFT_p[], std::complex<double>  p[]){
    int n = end - begin;
    std::complex<double> newp[n], newfftp[n];
    for (int i = 0; i < n; i++){
        newp[i] = p[i+begin];
    }
    fft(newfftp, newp, n);

    for (int i = 0; i < n; i++){
        FFT_p[i+begin] = newfftp[i];
    }
}

void pfft(std::complex<double>  FFT_p[], std::complex<double>  p[], int n, int num_threads){
    int order[n];
    for (int i = 0; i < n; i++){
        order[i] = reverseBits(i);
    }

    std::complex<double>  ordered_p[n];
    for (int i = 0; i < n; i++){
        ordered_p[i] = p[order[i]];
    }

    int thread_n = n / num_threads;
    std::vector<std::thread> l_thread(num_threads - 1);
    int begin,end;
    begin = 0;
    end = thread_n;
    for (int i = 0; i < num_threads-1; i++){
        l_thread[i] = std::thread(&thread_fft, begin, end, FFT_p, ordered_p);
        begin += thread_n;
        end += thread_n;
    }

    for (int i = 0; i < num_threads-1; i++){
        l_thread[i].join();
    }
    thread_fft(begin, end, FFT_p, ordered_p);
    
}


int main(){
    int n=8;
    std::complex<double>  p[n]{std::complex<double>(0,0),std::complex<double>(1,1),std::complex<double>(3,3),std::complex<double>(4,4),
                               std::complex<double>(4,4),std::complex<double>(3,3),std::complex<double>(1,1),std::complex<double>(0,0),};
    std::complex<double>  fft_p[n];
    std::cout << "origin" << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << p[i] << std::endl;
    }
    fft(fft_p,p,n);
    std::cout << "fft" << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << fft_p[i] << std::endl;
    }

    thread_fft(0,8,fft_p,p);
    std::cout << "thread_fft" << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << fft_p[i] << std::endl;
    }
    return 0;

}