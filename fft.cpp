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

void order_fft(std::complex<double>  FFT_p[], std::complex<double>  p[], int n){
    if (n==1){
        for (int i=0; i<n; i++){
            FFT_p[i] = p[i];
        }
        return;
    }
    std::complex<double>  U[n/2],V[n/2],FFT_U[n/2],FFT_V[n/2];
    for (int i=0; i<n/2; i++){
        U[i] = p[i];
        V[i] = p[i+n/2];
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

    /**
    std::cout << "thread_fft newp" << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << newp[i];
    }
    std::cout << ' '<< std::endl;
    **/

    order_fft(newfftp, newp, n);

    /**
    std::cout << "thread_fft newfftp" << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << newfftp[i];
    }
    std::cout << ' '<< std::endl;
    **/

    for (int i = 0; i < n; i++){
        FFT_p[i+begin] = newfftp[i];
    }
}

void pfft(std::complex<double>*  FFT_p, std::complex<double>*  p, int n, int num_threads){
    
    if (num_threads ==1){
        fft(FFT_p,p,n);
        return;
    } 

    int order[n];
    for (int i = 0; i < n; i++){
        order[i] = reverseBits(i, log2(n));
    }

    /**
    for (int i = 0; i < n; i++){
        std::cout << order[i];
    }
    std::cout << "order" << std::endl;
    **/
    
    std::complex<double>  ordered_p[n];
    for (int i = 0; i < n; i++){
        ordered_p[i] = p[order[i]];
    }

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
    thread_fft(begin, end, FFT_p, ordered_p);
    
    /**
    for (int i = 0; i < n; i++){
        std::cout << FFT_p[i];
    }
    std::cout << "pfft FFT p" << std::endl;
    **/
    
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

            std::complex<double>  w,w_n;
            w_n = cos(2*pi/m) + 1j * sin(2*pi/m);
            w = 1;

            for (int i = 0; i<m/2; i++){
                FFT_p[begin+i] = FFT_U[i] + w*FFT_V[i];
                FFT_p[begin+i + m/2] = FFT_U[i] - w*FFT_V[i];
                w = w * w_n;
            }

            begin +=m;
            end +=m;
        }
        m *=2;
    }
    
}


int main(){
    int n=8;
    std::complex<double>  p[n]{std::complex<double>(0,0),std::complex<double>(1,1),std::complex<double>(3,3),std::complex<double>(4,4),
                               std::complex<double>(4,4),std::complex<double>(3,3),std::complex<double>(1,1),std::complex<double>(0,0)};
    std::complex<double>  fft_p[n];
    /**
    std::cout << "origin" << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << p[i] << std::endl;
    }
    **/
    fft(fft_p,p,n);
    std::cout << "fft" << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << fft_p[i] ;
    }
    std::cout << ' '<< std::endl;

    std::complex<double>  p2[n]{std::complex<double>(0,0),std::complex<double>(1,1),std::complex<double>(3,3),std::complex<double>(4,4),
                               std::complex<double>(4,4),std::complex<double>(3,3),std::complex<double>(1,1),std::complex<double>(0,0),};
    std::complex<double>  fft_p2[n];
    pfft(fft_p2,p2,n,8);
    
    std::cout << "pfft" << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << fft_p2[i];
    }
    std::cout << ' '<< std::endl;
    return 0;

}