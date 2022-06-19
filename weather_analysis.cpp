#include "fft.cpp"
using namespace std;
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <numeric>
#include <algorithm>


pair<double,double> error(complex<double>* p1, complex<double>* p2, int n){
    double total_error = 0;
    double average_error = 0;
    for (int i = 0; i < n; i ++){
        total_error += pow(p1[i].real() - p2[i].real(), 2);
        average_error += abs(p1[i].real() - p2[i].real());
    }
    pair<double,double> Pair;
    Pair.first = sqrt(total_error);
    Pair.second = average_error / n;
    return Pair;
}

void compress(complex<double>* p, int n, int num_indexes){
    vector<complex<double>> vect(n);
    for (int i = 0; i < n; i ++){
        vect[i] = p[i];
    } // build vector so we can use builtin functions

    vector<int> idxs(num_indexes);
    iota(idxs.begin(), idxs.end(), 0);
    stable_sort(idxs.begin(), idxs.end(), [&vect](size_t i1, size_t i2) {return norm(vect[i1]) < norm(vect[i2]);});

    for (int i = 0; i < n; i++){
        if (std::count(idxs.begin(), idxs.end(), i) != 1){
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

int main3()
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
    auto start = std::chrono::steady_clock::now();
    fft(fft_p, p, n); // run fft on this array. fft_p is our output.
    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for fft is " << elapsed << " microseconds" << std::endl;

    start = std::chrono::steady_clock::now();
    pfft(pfft_p, p, n, 4); // run pfft on this array with 4 threads. pfft_p is our output.
    finish = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for parallel fft is " << elapsed << " microseconds" << std::endl;









    // here we compress fft_p and pfft_p to keep only large values.
    // we first keep the original file to compare
    complex<double> og_fft[n];
    for (int z = 0; z < n;z ++){
        og_fft[z] = fft_p[z];
    }
    int num_indexes = 1000;
    compress(fft_p, n, num_indexes); // keep 20 largest values
    compress(pfft_p, n, num_indexes); // same


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

    start = std::chrono::steady_clock::now();
    ifft(approx, fft_p, n);
    finish = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for ifft is " << elapsed << " microseconds" << std::endl;

    start = std::chrono::steady_clock::now();
    ifft_parallel(p_approx, pfft_p, n, 4);
    finish = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for parallel ifft is " << elapsed << " microseconds" << std::endl;

    ifft(og_approx, og_fft, n);


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

    pair<double, double> og_error = error(og_approx, p, n);
    pair<double, double> fft_error = error(approx, p, n);
    pair<double, double> pfft_error = error(p_approx, p, n);


    cout << "First comparisons. Squared error and average error ~~~~~~~~~~~~~~~~~~~~~ " << endl;
    cout << "Exact fft then ifft data: " << og_error.first << ", " << og_error.second << endl; // should be small
    cout << "compressed fft data: " << fft_error.first << ", " << fft_error.second << endl;
    cout << "compressed parallel fft data: " << pfft_error.first << ", " << pfft_error.second << endl;



    num_indexes = 100;
    complex<double> inverse_100[n], inverse_50[n], inverse_20[n];
    compress(pfft_p, n, num_indexes);
    num_indexes = 50;
    ifft_parallel(inverse_100, pfft_p, n, 4);
    compress(pfft_p, n, num_indexes);
    ifft_parallel(inverse_50, pfft_p, n, 4);
    num_indexes = 20;
    compress(pfft_p, n, num_indexes);
    ifft_parallel(inverse_20, pfft_p, n, 4);

    pair<double, double> error_1000 = pfft_error;
    pair<double, double> error_100 = error(inverse_100, p, n);
    pair<double, double> error_50 = error(inverse_50, p, n);
    pair<double, double> error_20 = error(inverse_20, p, n);



    cout << "Varying size of compression. Squared error and average error ~~~~~~~~~~~~~~" << endl;
    cout << "n = 1000: " << error_1000.first << ", " << error_1000.second << endl;
    cout << "error for n = 100: " << error_100.first << ", " << error_100.second << endl;
    cout << "error for n = 50: " << error_50.first << ", " << error_50.second << endl;
    cout << "error for n = 20: " << error_20.first << ", " << error_20.second << endl;

    cout << "Estimated temperature by hour. Correct vs approximation." << endl;
    cout << "t = 1000. aka: " << time_content[1000][1] << ". " << p[1000].real() << ", " << approx[1000].real() << endl;
    cout << "t = 2000. aka: " << time_content[2000][1] << ". " << p[2000].real() << ", " << approx[2000].real() << endl;
    cout << "t = 3000. aka: " << time_content[3000][1] << ". " << p[3000].real() << ", " << approx[3000].real() << endl;
    cout << "t = 4000. aka: " << time_content[4000][1] << ". " << p[4000].real() << ", " << approx[4000].real() << endl;
    cout << "t = 5000. aka: " << time_content[5000][1] << ". " << p[5000].real() << ", " << approx[5000].real() << endl;
    cout << "t = 6000. aka: " << time_content[6000][1] << ". " << p[6000].real() << ", " << approx[6000].real() << endl;
    cout << "t = 7000. aka: " << time_content[7000][1] << ". " << p[7000].real() << ", " << approx[7000].real() << endl;
    cout << "t = 8000. aka: " << time_content[8000][1] << ". " << p[8000].real() << ", " << approx[8000].real() << endl;

















    return 0;
}