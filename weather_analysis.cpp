#include "fft.cpp"
using namespace std;
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>



int load_and_test()
{
    string temp_name, date_name;
    temp_name = "new_york_temp.csv";
    date_name = "datetime.csv";

    vector<vector<double>> temp_content;
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

    vector<vector<string>> time_content;
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


    unsigned long long highest_power_2 = 1, n = temp_content.size();
//    while (highest_power_2 * 2 <= n) {
//        highest_power_2 *= 2;
//    }
    highest_power_2 = 8192;

    complex<double> p[highest_power_2], fft_p[highest_power_2];
    for (int z = 0; z < highest_power_2; z ++){
        p[z] = complex<double> (temp_content[z][1]);
    }

    fft(fft_p, p, highest_power_2);







    return 0;
}
