#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>      //sig. figures
#include <algorithm>    //sort()
#include <dirent.h>     //get file names
#include <sys/types.h>
#include <limits>
using namespace std;

// variabili globali
const string ifile = "./data/period.txt";
const double dt = 0.02;

// stampa su console
template <typename T> void Cell(T text){
    cout<< setw(15) << right << text;
}
void PrintEoF(char sep){
    cout<<endl << string(45, sep) <<endl <<endl;
}

// prototipi funzioni matematiche
double Mean(vector<double>);
double Summation(vector<double>);
double GetStdDev(vector<double>, int);
double CoeffCorr(vector<double>, vector<double>);

// struttura set di data
struct sample_struct {
    double freq;
    vector<double> thetas;
    double limit;
    vector<double> periods;
    sample_struct(string line){
        istringstream input_line(line);
        input_line >> freq;
        double theta;
        while (input_line >> theta)
            thetas.push_back(theta * 2.0 * M_PI);
        cout<< "Frequency = " << freq << "Hz" <<endl<<endl;
        GetLimit();
        Periods();
        PrintEoF('#');
    }

    void GetLimit(){
        vector<double> squared;
        for(double c : thetas)
            squared.push_back(c * c);
        double mean = Mean(squared);
        limit = sqrt(mean) / 2.0;
    }

    void Periods(){
        cout<< "Periods" <<endl; Cell("Time"); Cell("Corr."); Cell("Length"); Cell("Position"); cout<<endl;
        int i = 0;
        while(i < thetas.size()){
            vector<double> ys, xs;
            while(abs(thetas[i]) < limit){
                ys.push_back(thetas[i]);
                xs.push_back(i * dt);
                i++;
            }
            i++;
            if (ys.size() > 1){
                double prev = 0.0;
                if (!periods.empty())
                    prev = periods.back();
                double time = Interpol(xs, ys) - prev;
                periods.push_back(time);
                double corr = CoeffCorr(xs, ys);
                int length = ys.size();
                Cell(time); Cell(corr); Cell(length); Cell(i); cout<<endl;
            }
        }
    }

    double Interpol(vector<double> xs, vector<double> ys){
        double xone = 0.0, xtwo = 0.0, yone = 0.0, xy = 0.0;
        int data_size = xs.size();
        for (int i = 0; i < data_size; i++){
            xone += xs[i];
            xtwo += pow(xs[i], 2.0);
            yone += ys[i];
            xy += xs[i] * ys[i];
        }
        double delta = (data_size * xtwo) - (xone * xone);
        double coeff = ((data_size * xy) - (xone * yone)) / delta;
        double rcept = ((xtwo * yone) - (xone * xy)) / delta;
        /*k_line_sig = length_sigma * sqrt(data_size / delta);
        q_line_sig = length_sigma * sqrt(xtwo / delta);
        double num = 0.0;
        for (int i = 0; i < data_size; i++)
            num += pow ((lengths[i] - (k_line * d_forces[i]) - q_line), 2.0);
        sigma_post = sqrt(num / (data_size - 2.0));
        k_line_sig_p = sigma_post * sqrt(data_size / delta);
        q_line_sig_p = sigma_post * sqrt(xtwo / delta); */
        double interpol = -1.0 * rcept / coeff;
        return interpol;
    }
};

// prototipi funzioni main
vector<sample_struct> ReadFile();

int main(){
    vector<sample_struct> samples = ReadFile();
    
    return 0;
}

// funzioni main
vector<sample_struct> ReadFile(){
    vector<sample_struct> samples;
    ifstream input_file(ifile);
    if (!input_file.is_open()){
        cout<< "Error opening the file" <<endl;
        exit(1);
    }
    string line;
    while(getline(input_file, line)){
        sample_struct sample = sample_struct(line);
        samples.push_back(sample);
    }
    return samples;
}

// funzioni matematiche
double Mean(vector<double> data){
    double sum = 0.0;
    for (auto c : data)
        sum += c;
    double mean = sum / data.size();
    return mean;
}

double Summation(vector<double> data){
    double summation = 0.0;
    double m = Mean(data);
    for (auto c : data)
        summation += pow((c - m), 2.0);
    return summation;
}

double GetStdDev(vector<double> data, int data_len){
    double std_dev = sqrt(Summation(data) / data_len);
    return std_dev;
}

double CoeffCorr(vector<double> xs, vector<double> ys){
    vector<double> xys;
    int data_len = xs.size();
    for (int i = 0; i < data_len; i++)
        xys.push_back(xs[i] * ys[i]);
    double xs_mean = Mean(xs);
    double ys_mean = Mean(ys);
    double xys_mean = Mean(xys);
    double var_x = GetStdDev(xs, data_len);
    double var_y = GetStdDev(ys, data_len);
    double corr_coeff = (xys_mean - xs_mean * ys_mean) / (var_x * var_y);
    return corr_coeff;
}