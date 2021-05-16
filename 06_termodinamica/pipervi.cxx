#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>      //sig. figures
#include <algorithm>    //sort()
#include <dirent.h>     //get file names
#include <sys/types.h>
#include <limits>

using namespace std;

// variabili globali
const string idir = "./data";
const double p_conv = 10000.0, v_conv = 0.000001, t_conv = 273.15;

// stampa su console
template <typename T> void Cell(T text){
    cout<< setw(15) << right << text;
}
void PrintEoF(char sep){
    cout<<endl << string(75, sep) <<endl <<endl;
}

// prototipi funzioni matematiche
double Mean(vector<double>);
double WeightedMean(vector<double>, vector<double>);
double MeanSigma(vector<double>);
double Summation(vector<double>);
double GetStdDev(vector<double>, int);
double CoeffCorr(vector<double>, vector<double>);
void ThreeSigma(vector<double>&, vector<double>&);
array<double, 4> Interpol(vector<double>, vector<double>);

//strutture
struct sample_struct{
    int temp;
    vector<double> ps, vs, ts;
    sample_struct(string filename){
        ReadFile(filename);
        GetLines();
    }

    void ReadFile(string filename){
        ifstream ifile(filename);
        if (!ifile.is_open()){
            cout<< "Error opening the file" <<endl;
            exit(1);
        }
        double pi, vi, ti;
        while ((ifile >> pi) && (ifile >> vi) && (ifile >> ti)){
            ps.push_back(1.0 / (pi * p_conv));
            vs.push_back(vi * v_conv);
            ts.push_back(ti + t_conv);
        }
        temp = round(ts[0] - t_conv);
        ifile.close();
    }

    void GetLines(){
        cout<< "Temperature: " << temp <<endl<<endl;
        array<double, 4> line = Interpol(ps, vs);
        for (double c : line){
            Cell(c);
        }
        cout<<endl;
        line = Interpol(vs, ps);
        for (double c : line){
            Cell(c);
        }
        PrintEoF('-');
    }
};

//prototipi funzioni main
vector<string> GetFiles(string);

int main(){
    vector<string> filenames = GetFiles(idir);        // lettura delle cartelle con i dati
    vector<sample_struct> samples;
    for (string c : filenames){
        sample_struct sample = sample_struct(c);
        samples.push_back(sample);
    }
    return 0;
}


vector<string> GetFiles(string wdir){
    vector<string> filenames;
    const char *wdirname = wdir.c_str();
    struct dirent *entry;
    DIR *dir = opendir(wdirname);

    while ((entry = readdir(dir)) != NULL){
        string filename = wdir + "/" + entry->d_name;
        filenames.push_back(filename);
    }
    closedir(dir);
    sort(filenames.begin(), filenames.end());
    reverse(filenames.begin(), filenames.end());
    filenames.pop_back();
    filenames.pop_back();
    reverse(filenames.begin(), filenames.end());
    cout<< "Reading directory: " << wdir << '/' <<endl;
    for (string file : filenames)
        cout<< file <<endl;
    cout<<"Total files found: " <<filenames.size() <<endl;
    PrintEoF('#');
    return filenames;
}

// funzioni matematiche
double Mean(vector<double> data){
    double sum = 0.0;
    for (double c : data)
        sum += c;
    double mean = sum / data.size();
    return mean;
}

double WeightedMean(vector<double> values, vector<double> sigmas){
    double num = 0.0, den = 0.0, sig = 0.0;
    int N = values.size();
        for (int i = 0; i < N; i++){
            double sq_sig = pow(sigmas[i], 2.0);
            num += values[i] / sq_sig;
            den += 1.0 / sq_sig;
            sig += sq_sig;
        }
        double mean = num / den;
        return mean;
}

double MeanSigma(vector<double> sigmas){
    double sig = 0.0;
    for(double sigma : sigmas)
        sig += sigma * sigma;
    double MeanSigma = sqrt(sig) / sigmas.size();
    return MeanSigma;
}

double Summation(vector<double> data){
    double summation = 0.0;
    double m = Mean(data);
    for (double c : data)
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

void ThreeSigma(vector<double> &data, vector<double> &sigmas){
    vector<double> old_data, new_data(data);
    vector<double> old_sigmas, new_sigmas(sigmas);
    do {
        old_data = new_data;
        old_sigmas = new_sigmas;
        new_data.clear();
        new_sigmas.clear();
        double three_sigma = 3.0 * GetStdDev(old_data, (old_data.size() - 1));
        double mean = Mean(old_data);
        for (int i = 0; i < old_data.size(); i++)
            if(((mean - three_sigma) < old_data[i]) && (old_data[i] < (mean + three_sigma))){
                new_data.push_back(old_data[i]);
                new_sigmas.push_back(old_sigmas[i]);
            }
    } while (new_data.size() < old_data.size());
    data = new_data;
    sigmas = new_sigmas;
}

array<double, 4> Interpol(vector<double> xs, vector<double> ys){
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
    double num = 0.0;
    for (int i = 0; i < data_size; i++)
        num += pow ((ys[i] - (coeff * xs[i]) - rcept), 2.0);
    double sigma_post = sqrt(num / (data_size - 2.0));
    double coeff_s = sigma_post * sqrt(data_size / delta);
    double rcept_s = sigma_post * sqrt(xtwo / delta);
    array<double, 4> out = {coeff, coeff_s, rcept, rcept_s};
    return out;
}