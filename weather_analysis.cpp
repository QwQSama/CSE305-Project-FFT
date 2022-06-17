#include "fft.cpp"
using namespace std;
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <numeric>
#include <algorithm>


double error(complex<double>* p1, complex<double>* p2, int n){
    double total_error = 0;
    for (int i = 0; i < n; i ++){
        total_error += pow(p1[i].real() - p2[i].real(), 2);
    }
    return sqrt(total_error);
}

void compress(complex<double>* p, int n){
    vector<complex<double>> vect(n);
    for (int i = 0; i < n; i ++){
        vect[i] = p[i];
    } // build vector so we can use builtin functions

    vector<int> idxs(n);
    iota(idxs.begin(), idxs.end(), 0);
    stable_sort(idxs.begin(), idxs.end(), [&vect](size_t i1, size_t i2) {return vect[i1].real() < vect[i2].real();});

    for (int i = 0; i < n; i++){
        if (std::count(idxs.begin(), idxs.end(), i) == 0){
            p[i] = 0;
        }
    }

}

void retrieve_data(vector<vector<double>> &temp_content, vector<vector<string>> &time_content){
    string temp_name, date_name;
    temp_name = "new_york_temp.csv";
    date_name = "datetime.csv";

    vector<double> row;
    string line, word;

    fstream file (temp_name, ios::in);
    int i = 0;
    if(file.is_open())
    {
        while(getline(file, line))
        {
            row.clear();

            stringstream str(line);

            if (i != 0){
                while(getline(str, word, ','))
                    row.push_back(stod(word));
                temp_content.push_back(row);

            } else {
                i ++;
            }

        }
    }
    else
        cout<<"Could not open the file\n";

    vector<string> srow;
    line = "";
    word = "";

    fstream file2 (date_name, ios::in);
    i = 0;
    if(file2.is_open())
    {
        while(getline(file2, line))
        {
            srow.clear();

            stringstream str2(line);


            if (i != 0){
                while(getline(str2, word, ',')) {
                    srow.push_back(word);
                }

                time_content.push_back(srow);
            } else {
                i ++;
            }

        }
    }
    else
        cout<<"Could not open the file\n";



    cout << "temperature: "<< temp_content.size() << endl;
    cout << "temperature: "<< temp_content[0].size() << endl;
    cout << "temperature: "<< temp_content[0][0] << ' ' << temp_content[0][1] << endl;

    cout << "date: "<< time_content.size() << endl;
    cout << "date: "<< time_content[0].size() << endl;
    cout << "date: "<< time_content[0][0] << ' ' << time_content[0][1] << endl;

}

int main()
{
    vector<vector<double>> temp_content;
    vector<vector<string>> time_content;

    // retrieve temperature data into a vector and corresponding time values in a separate vector
    retrieve_data(temp_content, time_content);

    // This algorithm requires a power of 2 for n
    // we choose 8192 as their are 8,760 hours in a year and this corresponds to almost a full yearly cycle.
    unsigned long long n = 8192;

    complex<double> p[n], fft_p[n], pfft_p[n];
    for (int z = 0; z < n;z ++){
        p[z] = complex<double> (temp_content[z][1]); // build a complex array required for the algorithm
    }

    fft(fft_p, p, n); // run fft on this array. fft_p is our output.
    pfft(pfft_p, p, n, 4); // run pfft on this array with 4 threads. pfft_p is our output.


    // here we compress fft_p and pfft_p to keep only large values.
    // we first keep the original file to compare
    complex<double> og_fft[n];
    for (int z = 0; z < n;z ++){
        og_fft[z] = fft_p[z];
    }

    compress(fft_p, 20); // keep 20 largest values
    compress(pfft_p, 20); // same


    // now we write the results to a file
    ofstream fft_file("fft.csv");
    ofstream pfft_file("pfft.csv");
    fft_file << "fft values\n";
    for(int i=0; i < n; i++)
    {
        fft_file << fft_p[i].real() << endl;
    }
    pfft_file << "pfft values\n";
    for(int i=0; i<n; i++)
    {
        pfft_file << pfft_p[i].real() << endl;
    }


    // now let's perform inverse fft and retreive original data
    std::complex<double> approx[n];
    std::complex<double> p_approx[n];
    std::complex<double> og_approx[n];

    inverse_dft(og_approx, og_fft, n);
    inverse_dft(approx, fft_p, n);
    inverse_dft(p_approx, pfft_p, n);

    ofstream approx_file("approx.csv");
    approx_file << "approximation\n";
    for(int i=0; i<n; i++)
    {
        approx_file << approx[i].real() << endl;
    }
    ofstream p_approx_file("approx.csv");
    p_approx_file << "approximation\n";
    for(int i=0; i<n; i++)
    {
        p_approx_file << p_approx[i].real() << endl;
    }

    double og_error = error(og_approx, p, n);
    double fft_error = error(og_approx, p, n);
    double pfft_error = error(og_approx, p, n);

    cout << "error for fft then ifft data: " << og_error << endl; // should be small
    cout << "error for compressed fft data: " << fft_error << endl;
    cout << "error for compressed parallel fft data: " << pfft_error << endl;









    return 0;
}