#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <thread>
#include <chrono>
#define pi  3.14159265358979323846

void fft(std::complex<double>  *FFT_p, std::complex<double>  *p, int n){
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

    std::complex<double>  w;

    for (int i = 0; i<n/2; i++){
        w = std::polar(1.0, -2 * pi * i / n);
        FFT_p[i] = FFT_U[i] + w*FFT_V[i];
        FFT_p[i + n/2] = FFT_U[i] - w*FFT_V[i];
    }
}

void order_fft(std::complex<double>  order_p[], std::complex<double>  p[], int n, int m){
    if (m == 1){
        for (int i=0; i<n; i++){
            order_p[i] = p[i];
        }
        return;
    }
    std::complex<double>  U[n/2],V[n/2],order_U[n/2],order_V[n/2];
    for (int i=0; i<n/2; i++){
        U[i] = p[2*i];
        V[i] = p[2*i+1];
    }
    order_fft(order_U,U,n/2,m/2);
    order_fft(order_V,V,n/2,m/2);

    for (int i=0; i<n/2; i++){
        order_p[i] = order_U[i];
        order_p[i+n/2] = order_V[i];
    }
}

int reverseBits(int n, int digit) {
    int res = 0;
    //std::cout << n << "digit" << digit << std::endl;
    for(int i = 0; i < digit; i++){
        res = (res << 1) | (n & 1);
        n = n >> 1;
    }
    //std::cout << res << std::endl;
    return res;
}

void thread_fft(int begin, int end, std::complex<double>*  FFT_p, std::complex<double>*  p){
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

void pfft(std::complex<double>*  FFT_p, std::complex<double>*  p, int n, int num_threads){
    /*
    if (num_threads ==1){
        fft(FFT_p,p,n);
        return;
    } */

    std::complex<double>  ordered_p[n];
    order_fft(ordered_p,p,n,num_threads);

    /**
    for (int i = 0; i < n; i++){
        std::cout << ordered_p[i];
    }
    std::cout << "order p" << std::endl;
    **/

    int thread_n = n / num_threads;
    std::vector<std::thread> l_thread(num_threads - 1);
    int begin,end;
    begin = 0;
    end = thread_n;
    for (int i = 0; i < num_threads-1; i++){
        l_thread[i] = std::thread(&thread_fft, begin, end, &FFT_p[0], &ordered_p[0]);
        begin += thread_n;
        end += thread_n;
    }

    for (int i = 0; i < num_threads-1; i++){
        l_thread[i].join();
    }
    thread_fft(begin, end, &FFT_p[0], &ordered_p[0]);

    int m = n/num_threads*2;
    for (int thread=0; thread < log2(num_threads); thread++){
        int begin,end;
        begin = 0;
        end = m;
        while (end <= n){
            std::complex<double>  FFT_U[m/2],FFT_V[m/2];
            for (int i=0; i<m/2; i++){
                FFT_U[i] = FFT_p[begin+i];
                FFT_V[i] = FFT_p[begin+i+m/2];
            }

            std::complex<double>  w;

            for (int i = 0; i<m/2; i++){
                w = std::polar(1.0, -2 * pi * i / n);
                FFT_p[begin+i] = FFT_U[i] + w*FFT_V[i];
                FFT_p[begin+i + m/2] = FFT_U[i] - w*FFT_V[i];
            }

            begin +=m;
            end +=m;
        }
        m *=2;
    }

}

void ifft(std::complex<double> *ifftp, std::complex<double> *fft_p, int n){
    std::complex<double> con_fft_p[n];
    for (int i =0; i < n; i++){
        con_fft_p[i] = std::conj(fft_p[i]);
    }

    fft(ifftp,con_fft_p,n);

    for (int j = 0; j <= n; j++){
        ifftp[j] = ifftp[j]/(n*1.0);
    }


    for (int i =0; i < n; i++){
        ifftp[i] = std::conj(ifftp[i]);
    }

}

void ifft_parallel(std::complex<double> *ifftp, std::complex<double> *fft_p, int n, int p ){
    std::complex<double> con_fft_p[n];
    for (int i =0; i < n; i++){
        con_fft_p[i] = std::conj(fft_p[i]);
    }

    pfft(ifftp,con_fft_p,n,p);

    for (int j = 0; j <= n; j++){
        ifftp[j] = ifftp[j]/(n*1.0);
    }


    for (int i =0; i < n; i++){
        ifftp[i] = std::conj(ifftp[i]);
    }

}


int main2(){
    int n=pow(2,3);
    int num_thread = 2;
    std::complex<double>  p0[8]{std::complex<double>(0,0),std::complex<double>(1,1),std::complex<double>(3,3),std::complex<double>(4,4),
                                std::complex<double>(4,4),std::complex<double>(3,3),std::complex<double>(1,1),std::complex<double>(0,0)};
    std::complex<double>  fft_p[n],p[n],pfft_p[n],offt_p[n],ifft_p[n],ifft_p_s[n];

    std::cout << "hi" << std::endl;
    int js=1;
    for (int i=0; i<n; i++){
        p[i] = std::complex<double>(js,0);
        //if (n/2-1 ==i){js = 0;}
        js++;
        p[i] = p0[i];
    }

    auto start = std::chrono::steady_clock::now();
    fft(fft_p,p,n);
    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for fft is " << elapsed << " microseconds" << std::endl;

    start = std::chrono::steady_clock::now();
    pfft(pfft_p,p,n,num_thread);
    finish = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for pfft is " << elapsed << " microseconds" << std::endl;

    start = std::chrono::steady_clock::now();
    ifft_parallel(ifft_p,fft_p,n,num_thread);
    finish = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for ifft_parallel is " << elapsed << " microseconds" << std::endl;

    start = std::chrono::steady_clock::now();
    ifft(ifft_p_s,fft_p,n);
    finish = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for ifft_simple is " << elapsed << " microseconds" << std::endl;

    std::cout << "test" << std::endl;
    bool f = true;
    for (int i = 0; i<n; i++){
        if (abs(ifft_p[i] - p[i])> 0.00001){
            f = false;
            std::cout <<  ifft_p_s[i] << p[i]<< std::endl;
        }
        std::cout <<  ifft_p_s[i] << p[i]<< std::endl;
    }
    if (!f){
        std::cout << "test not equal" << std::endl;
    }
    else{
        std::cout << "test pass" << std::endl;
    }
    return 0;

}