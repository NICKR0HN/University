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
const string idir = "./data";
const string odir = "./output";
const double f_factor = 4.0 * 9.806 * 0.001;
const double l_factor = 0.000001;

//funzioni utili per la stampa su console
void PrintEoF(){
    cout<<endl << string(72, '=') <<endl <<endl;
}
template <typename T> void Cell(T text){
    cout<< setw(14) << right << text;
}

double Compatibility(double, double, double, double);
double WeightedMean(vector<double>, vector<double>);
double MeanSigma(vector<double>);

struct sample_struct{
    string filename;
    vector<double> forces, lengths;
    double force_sigma, length_sigma;
    double k_line, q_line, k_line_sig, q_line_sig;
    double sigma_post, k_line_sig_p, q_line_sig_p;
    vector<double> k_vector, k_sigma_vector;
    double k_mean, k_mean_sigma, compatible;
    sample_struct(string filepath){
        filename = filepath;
        ReadFile();
        Interpol();
        KVector();
        compatible = Compatibility(k_line, k_mean, k_line_sig_p, k_mean_sigma);
        PrintData();
    }

    void ReadFile(){
        ifstream input_file(filename);
        if (!input_file.is_open()){
            cout<< "Error opening the file" <<endl;
            exit(1);
        }
        double force_sigma_r, length_sigma_r;
        input_file >> force_sigma_r;
        force_sigma = force_sigma_r * f_factor;
        input_file >> length_sigma_r;
        length_sigma = length_sigma_r * l_factor;
        double force_r, force, length_r, length;
        while ((input_file >> force_r) && (input_file >> length_r)){
            force = force_r * f_factor;
            forces.push_back(force);
            length = length_r * l_factor;
            lengths.push_back(length);
        }
        input_file.close();
    }

    void Interpol(){
        double xone = 0.0, xtwo = 0.0, yone = 0.0, xy = 0.0;
        int N = forces.size();
        for (int i = 0; i < N; i++){
            xone += forces[i];
            xtwo += pow(forces[i], 2.0);
            yone += lengths[i];
            xy += forces[i] * lengths[i];
        }
        double delta = (N * xtwo) - (xone * xone);
        k_line = ((N * xy) - (xone * yone)) / delta;
        q_line = ((xtwo * yone) - (xone * xy)) / delta;
        k_line_sig = length_sigma * sqrt(N / delta);
        q_line_sig = length_sigma * sqrt(xtwo / delta);
        double num = 0.0;
        for (int i = 0; i < N; i++)
            num += pow ((lengths[i] - (k_line * forces[i]) - q_line), 2.0);
        sigma_post = sqrt(num / (N - 2.0));
        k_line_sig_p = sigma_post * sqrt(N / delta);
        q_line_sig_p = sigma_post * sqrt(xtwo / delta);
    }
    
    void KVector(){
        double df, df_sigma, dx, dx_sigma, k, k_sigma;
        int k_size = floor((forces.size() / 2));
        dx_sigma = length_sigma * sqrt(2.0);
        df_sigma = force_sigma * sqrt(2.0);
        for (int i = 0; i < k_size; i++){
            dx = lengths[2 * i + 1] - lengths[2 * i];
            df = forces[2 * i + 1] - forces[2 * i];
            k = dx / df;
            k_vector.push_back(k);
            k_sigma = abs(k) * sqrt(pow((dx_sigma / dx), 2.0) + pow((df_sigma / df), 2.0));
            k_sigma_vector.push_back(k_sigma);
        }
        k_mean = WeightedMean(k_vector, k_sigma_vector);
        k_mean_sigma = MeanSigma(k_sigma_vector);
    }

    void PrintData(){
        cout<< filename <<endl<<endl;
        PrintTable();
        PrintEoF();
    }
    
    void PrintTable(){
        cout<< setprecision(5);
        cout<< setw(15) << left << ' ';             Cell("k (m/N)");    Cell("sigma (m/N)");    Cell("q (m)");  Cell("sigma (m)");  cout<<endl;
        cout<< setw(15) << left << "Interpol.";     Cell(k_line);       Cell(k_line_sig);       Cell(q_line);   Cell(q_line_sig);   cout<<endl;
        cout<< setw(15) << left << "Post. sigma";   Cell(k_line);       Cell(k_line_sig_p);     Cell(q_line);   Cell(q_line_sig_p); cout<<endl;
        cout<< setw(15) << left << "Mean k";        Cell(k_mean);       Cell(k_mean_sigma);     cout<<endl;
        cout<< setw(15) << left << "Compatibility"; Cell(compatible);   cout<<endl;
    }
};

//prototipi funzioni
vector<string> GetFiles(string);
vector<sample_struct> ElaborateData(vector<string>);
void Compare(vector<sample_struct>);

int main(){
    vector<string> foldernames = GetFiles(idir);        // lettura delle cartelle con i dati
    for (string folder : foldernames){
        vector<string> filenames = GetFiles(folder);
        vector<sample_struct> samples = ElaborateData(filenames);
        Compare(samples);
    }
    return 0;
}

// ottiene i nomi di tutti i file presenti nella cartella di input
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

vector<sample_struct> ElaborateData(vector<string> filenames){
    vector<sample_struct> samples;
    for (string file : filenames){
        sample_struct sample = sample_struct(file);
        samples.push_back(sample);
    }
    return samples;
}

void Compare(vector<sample_struct> samples){
    double comp_lines, comp_means;
    comp_lines = Compatibility(samples[0].k_line, samples[1].k_line, samples[0].k_line_sig_p, samples[1].k_line_sig_p);
    comp_means = Compatibility(samples[0].k_mean, samples[1].k_mean, samples[0].k_mean_sigma, samples[1].k_mean_sigma);
    cout << setw(20) << left << "Line compatibility";   Cell(comp_lines); cout<<endl;
    cout << setw(20) << left << "Mean compatibility";   Cell(comp_means); cout<<endl<<endl;
    PrintEoF();
}

double Compatibility(double x, double y, double x_s, double y_s){
    double comp, num, den;
    num = abs(x - y);
    den = sqrt((x_s * x_s) + (y_s * y_s));
    comp = num / den;
    return comp;
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
    double MeanSigma = sig / sigmas.size();
    return MeanSigma;
}