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

void PrintEoF(){
    cout<<endl << string(90, '=') <<endl <<endl;
}
template <typename T> void Cell(T text){
    cout<< setw(15) << right << text;
}

double Mean(vector<double>);
double WeightedMean(vector<double>, vector<double>);
double MeanSigma(vector<double>);
double Summation(vector<double>);
double GetStdDev(vector<double>, int);
double CoeffCorr(vector<double>, vector<double>);

struct est_struct{
    double k, k_sig;
    double length, length_sig;
    double diam, diam_sig;
    double young, young_sig;
    est_struct(string line){
        ReadData(line);
        YoungModulus();
        PrintData();
    }
    void ReadData(string line){
        istringstream input_line(line);
        double k_str, k_str_sig, k_shr, k_shr_sig;
        input_line >> k_str >> k_str_sig;
        input_line >> k_shr >> k_shr_sig;
        input_line >> length >> length_sig;
        input_line >> diam >> diam_sig;
        vector<double> ks {k_str, k_shr};
        vector<double> ks_sig {k_str_sig, k_shr_sig};
        k = WeightedMean(ks, ks_sig);
        k_sig = MeanSigma(ks_sig);
    }
    void YoungModulus(){
        young = (4.0 * length) / (M_PI * k * diam * diam);
        double a = (length_sig) / (k * diam * diam);
        double b = (length * k_sig) / (k * k * diam * diam);
        double c = (length * diam_sig) / (2.0 * k * diam * diam * diam);
        young_sig = (4.0 / M_PI) * sqrt((a * a) + (b * b) + (c * c));
    }
    void PrintData(){
        Cell(k); Cell(k_sig); Cell(young); Cell(young_sig); Cell(length); Cell(diam);
        cout<<endl;
    }
};


vector<string> GetFiles(string);
vector<est_struct> ReadFile(string);
array<double,2> DataAnalysis(vector<est_struct>);
    double CorrCheck(vector<est_struct>);
    array<double,2> MeanYoung(vector<est_struct>);
void Summary(vector<double>, vector<double>);

int main(){
    vector<string> filenames = GetFiles(idir);        // lettura delle cartelle con i dati
    vector<double> youngs, youngs_sig; 
    for (string filename : filenames){
        vector<est_struct> estens = ReadFile(filename);
        array<double,2> out = DataAnalysis(estens);
        youngs.push_back(out[0]);
        youngs_sig.push_back(out[1]);        
    }
    Summary(youngs, youngs_sig);
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
    PrintEoF();
    return filenames;
}

vector<est_struct> ReadFile(string filename){
        cout<< filename <<endl<<endl;
        Cell("k"); Cell("k.sigma"); Cell("Y.modulus"); Cell("Y.mod.sig."); Cell("Length"); Cell("Diam.");
        cout<<endl;
        vector<est_struct> data;
        ifstream input_file(filename);
        if (!input_file.is_open()){
            cout<< "Error opening the file" <<endl;
            exit(1);
        }
        string line;
        while(getline(input_file, line)){
            est_struct est = est_struct(line);
            data.push_back(est);
        }
        cout<<endl << string(90, '-') <<endl <<endl;
        return data;
    }

array<double,2> DataAnalysis(vector<est_struct> estens){
    double corr = CorrCheck(estens);
    cout<< corr <<endl;
    array<double,2> mean_young = MeanYoung(estens);
    cout<< "Mean young modulus: " << mean_young[0] << "\t+-" << mean_young[1] <<endl;
    PrintEoF();
    return mean_young;
}

void Summary(vector<double> youngs, vector<double> youngs_sig){
    double young_f = WeightedMean(youngs, youngs_sig);
    double young_f_sig = MeanSigma(youngs_sig);
    cout<< "Final young modulus: " << young_f << "\t+- " << young_f_sig <<endl;
    PrintEoF();
}

double CorrCheck(vector<est_struct> estens){
    vector<double> ks, ys;
    if (estens[0].diam == estens[1].diam){
        cout<< "Constant diameter: ";
        for (est_struct est : estens){
            ks.push_back(est.k);
            ys.push_back(est.length);
        }
    }
    else if (estens[0].length == estens[1].length){
        cout<< "Constant length: ";
        for (est_struct est : estens){
            ks.push_back(est.k);
            ys.push_back(1.0 / pow(est.diam, 2.0));
        }
    }
    else {
        cout<< "Error: no match was found"<<endl;
        PrintEoF();
        return nan("");
    }
    double corr = CoeffCorr(ks, ys);
    return corr;
}

array<double,2> MeanYoung(vector<est_struct> estens){
    array<double,2> out;
    vector<double> youngs, youngs_sig;
    for(est_struct est : estens){
        youngs.push_back(est.young);
        youngs_sig.push_back(est.young_sig);
    }
    out[0] = WeightedMean(youngs, youngs_sig);
    out[1] = MeanSigma(youngs_sig);
    return out;
}

double Mean(vector<double> data){
    double sum = 0.0;
    for (auto c : data)
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