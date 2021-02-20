#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>      //sig. figures
#include <algorithm>    //sort()
#include <dirent.h>     //get file names
#include <sys/types.h>

using namespace std;

// dichiarazione della struttura per ciascun campione
struct sample_struct {
    // variabili stabili della struttura
    string filename, filepath;
    vector<double> data;
    double distance;
    // costruttore
    sample_struct(string name){
        filename = name;
        filepath = "./data/"+name;
        ReadFile();
    }

    // acquisizione dei dati
    void ReadFile(){
        ifstream input_file(filepath);
        if (!input_file.is_open()){
            cout<< "Error opening the file";
            return;
        }
        double value;
        input_file >> distance;
        while (input_file >> value)
            data.push_back(value);
        input_file.close();
    }

    // variabili non costanti della struttura
    // vengono calcolate quando servono senza salvarle, in modo che siano sempre aggiornate
    double Max(){
        double max = data[0];
        for (auto c : data)
            if (c > max)
                max = c;
        return max;
    }
    double Min(){
        double min = data[0];
        for (auto c : data)
            if (c < min)
                min = c;
        return min;
    }
    double Mean(){
        double sum = 0.0;
        for (auto c : data)
            sum += c;
        double mean = sum / data.size();
        return mean;
    }
    double Median(){
        vector<double> sort_data = data;
        sort(sort_data.begin(), sort_data.end());
        int half = sort_data.size() / 2;
        if(sort_data.size()%2)
            return sort_data[half + 1];
        return (sort_data[half] + sort_data[half + 1]) / 2;
    }

    // sommatoria del quadrato degli scarti
    double Summation(){
        double summation = 0.0;
        double m = Mean();
        for (auto c : data)
            summation += pow((c - m), 2.0);
        return summation;
    }
    // formula generica della deviazione standard
    double GetStdDev(int data_len){
        double std_dev = sqrt(Summation() / data_len);
        return std_dev;
    }
    // deviazione standard campionaria
    double StdDev(){
        double std_dev = GetStdDev(data.size());
        return std_dev;
    }
    // deviazione standard singola misura
    double StdDevCorr(){
        double std_dev_corr = GetStdDev(data.size() - 1);
        return std_dev_corr;
    }
    // deviazione standard sulla media
    double StdDevMean(){
        double std_dev_mean = StdDevCorr()/sqrt(data.size());
        return std_dev_mean;
    }

    // eliminazione dei dati secondo la regola del 3-sigma
    vector<double> Refine(){
        vector<double> new_data, removed_data;
        double three_sigma = 3.0 * StdDevCorr();
        double mean = Mean();
        for (auto c : data){
            if (c < mean - three_sigma || c > mean + three_sigma)
                removed_data.push_back(c);
            else new_data.push_back(c);
        }
        data = new_data;
        return removed_data;
    }

    // stampa delle informazioni sul campione sulla console
    void PrintData(){
        cout<< setprecision(4);
        cout<< "Data set size: "    << data.size()  << "\t\tDistance: "             << distance     <<endl;
        cout<< "Minimum value: "    << Min()        << "\t\tMaximum value: "        << Max()        <<endl;
        cout<< "Mean value: "       << Mean()       << "\t\tMedian value: "         <<Median()      <<endl;
        cout<< "Std. deviation: "   << StdDev()     << "\tCorrected std. dev.: " << StdDevCorr() 
            << "\t\tMean std. dev.: "      << StdDevMean()   <<endl <<endl;
        cout<< string(100, '-') <<endl <<endl;
        cout<< setprecision(0);
    }
};

vector<string> GetFiles();
vector<sample_struct> ElaborateData(vector<string>);
void PrintEoF();


int main(){
    vector<string> filenames = GetFiles();
    vector<sample_struct> samples = ElaborateData(filenames);
    return 0;
}

vector<string> GetFiles(){
    vector<string> filenames;
    struct dirent *entry;
    DIR *dir = opendir("./data");

    while ((entry = readdir(dir)) != NULL)
        filenames.push_back(entry->d_name);
    closedir(dir);
    filenames.pop_back();
    filenames.pop_back();
    cout<< "Reading files:" <<endl;
    for (auto c : filenames)
        cout<< c <<endl;
    cout<<"Total: " <<filenames.size() <<endl;
    PrintEoF();
    return filenames;
}

vector<sample_struct> ElaborateData(vector<string> filenames){
    vector<sample_struct> data_out;
    for (auto c : filenames){
        sample_struct sample = sample_struct(c);
        cout<< sample.filename <<endl <<endl;
        vector<double> removed_data;
        do {
            cout<< setprecision(4);
            removed_data = sample.Refine();
            if (removed_data.size()){
                cout<< "Removed data:" <<endl;
                for (auto d : removed_data)
                    cout<< d << '\t';
                cout<< endl;
            }
            else cout<< "All data is compatible" <<endl <<endl;
        } while (removed_data.size());
        sample.PrintData();
    }
    return data_out;
}

void PrintEoF(){
    cout<<endl << string(100, '=') <<endl <<endl;
}