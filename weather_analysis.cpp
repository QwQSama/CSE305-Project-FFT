#include "fft.cpp"
using namespace std;
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

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

    retrieve_data(temp_content, time_content);


//    while (highest_power_2 * 2 <= n) {
//        highest_power_2 *= 2;
//    }
    unsigned long long n = 8192;

    complex<double> p[n], fft_p[n];
    for (int z = 0; z < n;z ++){
        p[z] = complex<double> (temp_content[z][1]);
    }

    fft(fft_p, p, n);

    for (int z = 0; z < n; z++){
        if (std::abs(fft_p[z]) < 50.) {
            fft_p[z] = 0.;
        }
    }
    ofstream myfile("fft_csv.csv");
    int vsize = 8192;
    myfile << "fft values\n";
    for(int n=0; n<vsize; n++)
    {
        myfile << fft_p[n].real() << endl;
    }

    std::complex<double> approx[8192];
    inverse_dft(approx, fft_p, 8192);

    ofstream myfile2("approx.csv");
    int vsize2 = 8192;
    myfile2 << "approximation\n";
    for(int n=0; n<vsize2; n++)
    {
        myfile2 << approx[n].real() << endl;
    }










    return 0;
}