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
const string odir = "./results";
const string ofile = "./output.txt";
const double p_conv = 98100.0, v_conv = 0.000001, t_conv = 273.15;
const double t_amb = 20.0 + t_conv, t_amb_s = 0.1;
const double R = 8.314462618;

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
array<double, 5> Interpol(vector<double>, vector<double>);
vector<double> LineDist(vector<double>, vector<double>, double, double);
string GetFileName(string);

//strutture
struct sample_struct{
    double temp, t_avg, t_avg_s;
    vector<double> ps, vs, ts;
    double a1, a1_s, b1, b1_s, post1, a2, a2_s, b2, b2_s, post2, corr;
    vector<double> rests_1, rests_2;
    double vol1, vol1_s, vol2, vol2_s, vol0;
    vector<double> moles;
    double mol, mol_s;
    sample_struct(string filename){
        ReadFile(filename);
        GetLines();
        GetRests();
        GetVolumes();
        GetMoles();
        PrintConsole();
        WriteFile(filename);
    }

    void ReadFile(string filename){
        ifstream ifile(filename);
        if (!ifile.is_open()){
            cout<< "Error opening the file" <<endl;
            exit(1);
        }
        double pi, vi, ti;
        while ((ifile >> pi) && (ifile >> vi) && (ifile >> ti)){
            if((4.8 < vi) && (vi < 23.0)){
                ps.push_back(1.0 / (pi * p_conv));
                vs.push_back(vi * v_conv);
                ts.push_back(ti + t_conv);
            }
        }
        t_avg = Mean(ts);
        t_avg_s = GetStdDev(ts, ts.size() - 1);
        temp = t_avg - t_conv;
        ifile.close();
    }

    void GetLines(){
        array<double, 5> line;
        corr = CoeffCorr(ps, vs);
        line = Interpol(ps, vs);
        a1 = line[0]; a1_s = line[1];
        b1 = line[2]; b1_s = line[3];
        post1 = line[4];
        line = Interpol(vs, ps);
        a2 = line[0]; a2_s = line[1];
        b2 = line[2]; b2_s = line[3];
        post2 = line[4];
    }

    void GetRests(){
        rests_1 = LineDist(ps, vs, a1, b1);
        rests_2 = LineDist(vs, ps, a2, b2);
    }

    void GetVolumes(){
        vol1 = -1.0 * b1;
        vol1_s = b1_s;
        vol2 = b2 / a2;
        double cov_ab = (Mean(vs) / vs.size()) * pow((post2 / GetStdDev(vs, vs.size())), 2.0);
        vol2_s = abs(vol2) * sqrt(pow(a2_s / a2, 2) + pow(b2_s / b2, 2) - (2.0 * cov_ab / (a2 * b2)));
    }

    void GetMoles(){
        double v_corr = vol1 * t_amb / t_avg;     //temporaneo
        for(int i = 0; i < vs.size(); i++){
            double mole = (vs[i] + v_corr) / (R * ts[i] * ps[i]);
            moles.push_back(mole);
        }
        mol = Mean(moles);
        mol_s = GetStdDev(moles, moles.size() - 1);
    }

    void PrintConsole(){
        cout<< "Temperature: " << temp <<endl;
        cout<< "Correlation: " << corr <<endl;
        PrintEoF('-');

    }

    void WriteFile(string filepath){
        ofstream fileout(ofile, ofstream::app);
        fileout<< temp <<'\t' << t_avg_s <<'\t' << corr
            <<'\t' << a1 <<'\t' << a1_s <<'\t' << b1 <<'\t' << b1_s <<'\t' << post1
            <<'\t' << a2 <<'\t' << a2_s <<'\t' << b2 <<'\t' << b2_s <<'\t' << post2
            <<'\t' << vol1 <<'\t' << vol1_s
            <<'\t' << vol2 <<'\t' << vol2_s
            <<'\t' << mol <<'\t' << mol_s <<endl;
        fileout.close();
        ofstream output(odir + GetFileName(filepath));
        if (!output.is_open()){
            cout<< "Permission denied" <<endl;
            exit(1);
        }
        for (int i = 0; i < rests_1.size(); i++){
            output<< rests_1[i] <<'\t' << rests_2[i] <<'\t' << moles[i] <<endl;
        }
        output.close();
    }
};

//prototipi funzioni main
vector<string> GetFiles(string);
void OverWriteFile();

int main(){
    OverWriteFile();
    vector<string> filenames = GetFiles(idir);        // lettura delle cartelle con i dati
    vector<sample_struct> samples;
    for (string c : filenames){
        sample_struct sample = sample_struct(c);
        samples.push_back(sample);
    }
    return 0;
}

void OverWriteFile(){
    ofstream overwrite(ofile, ofstream::trunc);
    if (!overwrite.is_open()){
        cout<< "Permission denied" <<endl;
        exit(1);
    }
    overwrite.close();
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

array<double, 5> Interpol(vector<double> xs, vector<double> ys){
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
    array<double, 5> out = {coeff, coeff_s, rcept, rcept_s, sigma_post};
    return out;
}

vector<double> LineDist(vector<double> xs, vector<double> ys, double coeff, double rcept){
    vector<double> dist;
    for (int i = 0; i < xs.size(); i++){
        dist.push_back(ys[i] - xs[i] * coeff - rcept);
    }
    return dist;
}

string GetFileName(string filepath){
    return filepath.substr(filepath.rfind('/'), filepath.size() - 1);
}