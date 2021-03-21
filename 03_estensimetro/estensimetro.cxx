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

//funzioni algebriche
double Compatibility(double, double, double, double);
double WeightedMean(vector<double>, vector<double>);
double MeanSigma(vector<double>);
double Mean(vector<double>);
double Summation(vector<double>);
double GetStdDev(vector<double>, int);

struct sample_struct{
    // variabili della struttura
    string filename;
    int data_size = 0;
    vector<double> forces, d_forces, lengths;
    double force_sigma, length_sigma;
    double k_line, q_line, k_line_sig, q_line_sig;
    double corr_coeff, sigma_post, k_line_sig_p, q_line_sig_p;
    array<double,2> k_mean, k_mean_sigma;
    array<int,2> times_k;
    double comp_a_0, comp_a_1, comp_0_1;
    double k_final, k_final_sigma;

    // costruttore
    sample_struct(string filepath){
        filename = filepath;
        ReadFile();
        for (double force : forces){
            double d_force = force - forces[0];
            d_forces.push_back(d_force);
        }
        Interpol();
        CoeffCorr();
        for (int i = 0; i < 2; i++)
            KVector(i);
        comp_a_0 = Compatibility(k_line, k_mean[0], k_line_sig_p, k_mean_sigma[0]);
        comp_a_1 = Compatibility(k_line, k_mean[1], k_line_sig_p, k_mean_sigma[1]);
        comp_0_1 = Compatibility(k_mean[0], k_mean[1], k_mean_sigma[0], k_mean_sigma[1]);
        FinalK();
        PrintData();
    }

    //funzioni della struttura
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
            data_size++;
        }
        input_file.close();
    }
    
    void Interpol(){
        double xone = 0.0, xtwo = 0.0, yone = 0.0, xy = 0.0;
        for (int i = 0; i < data_size; i++){
            xone += d_forces[i];
            xtwo += pow(d_forces[i], 2.0);
            yone += lengths[i];
            xy += d_forces[i] * lengths[i];
        }
        double delta = (data_size * xtwo) - (xone * xone);
        k_line = ((data_size * xy) - (xone * yone)) / delta;
        q_line = ((xtwo * yone) - (xone * xy)) / delta;
        k_line_sig = length_sigma * sqrt(data_size / delta);
        q_line_sig = length_sigma * sqrt(xtwo / delta);
        double num = 0.0;
        for (int i = 0; i < data_size; i++)
            num += pow ((lengths[i] - (k_line * d_forces[i]) - q_line), 2.0);
        sigma_post = sqrt(num / (data_size - 2.0));
        k_line_sig_p = sigma_post * sqrt(data_size / delta);
        q_line_sig_p = sigma_post * sqrt(xtwo / delta);
    }

    void CoeffCorr(){
        vector<double> xys;
        for (int i = 0; i < data_size; i++)
            xys.push_back(d_forces[i] * lengths[i]);
        double xs_mean = Mean(d_forces);
        double ys_mean = Mean(lengths);
        double xys_mean = Mean(xys);
        double var_x = GetStdDev(d_forces, data_size);
        double var_y = GetStdDev(lengths, data_size);
        corr_coeff = (xys_mean - xs_mean * ys_mean) / (var_x * var_y);
    }
    
    void KVector(int odd){
        double df, df_sigma, dx, dx_sigma, k, k_sigma;
        int count = 0;
        vector<double> k_vector, k_sigma_vector;
        dx_sigma = length_sigma * sqrt(2.0);
        df_sigma = force_sigma * sqrt(2.0);
        for (int i = odd; i < (data_size - 1); i = i + 2){
            dx = lengths[i + 1] - lengths[i];
            df = d_forces[i + 1] - d_forces[i];
            k = dx / df;
            k_vector.push_back(k);
            k_sigma = abs(k) * sqrt(pow((dx_sigma / dx), 2.0) + pow((df_sigma / df), 2.0));
            k_sigma_vector.push_back(k_sigma);
            count++;
        }
        k_mean[odd] = Mean(k_vector);
        k_mean_sigma[odd] = MeanSigma(k_sigma_vector);
        times_k[odd] = count;
    }

    void FinalK(){
        vector<double> ks {k_line}, ks_sigma {k_line_sig_p};
        if (comp_a_0 < comp_a_1){
            ks.push_back(k_mean[0]);
            ks_sigma.push_back(k_mean_sigma[0]);
        }
        else {
            ks.push_back(k_mean[1]);
            ks_sigma.push_back(k_mean_sigma[1]);
        }
        k_final = WeightedMean(ks, ks_sigma);
        k_final_sigma = MeanSigma(ks_sigma);
    }

    void PrintData(){
        cout<< filename <<endl<<endl;
        PrintTable();
        PrintEoF();
    }
    
    void PrintTable(){
        cout<< setprecision(5);
        cout<< setw(13) << left << "Data size = "   << data_size;   Cell("k (m/N)");    Cell("sigma (m/N)");    Cell("q (m)");  Cell("sigma (m)");  cout<<endl;
        cout<< setw(15) << left << "Interpol.";                     Cell(k_line);       Cell(k_line_sig);       Cell(q_line);   Cell(q_line_sig);   cout<<endl;
        cout<< setw(15) << left << "Post. sigma";                   Cell(k_line);       Cell(k_line_sig_p);     Cell(q_line);   Cell(q_line_sig_p); cout<<endl;
        cout<< setw(15) << left << "Mean k1";                       Cell(k_mean[0]);    Cell(k_mean_sigma[0]);  Cell(times_k[0]);                   cout<<endl;
        cout<< setw(15) << left << "Mean k2";                       Cell(k_mean[1]);    Cell(k_mean_sigma[1]);  Cell(times_k[1]);                   cout<<endl;
        cout<< setw(15) << left << "Final k";                       Cell(k_final);      Cell(k_final_sigma);    cout<<endl<<endl;

        cout<< setw(15) << left << "Post. sigma";                   Cell(sigma_post);   cout<<endl;
        cout<< setprecision(6);
        cout<< setw(15) << left << "Corr. coeff.";                  Cell(corr_coeff);   cout<<endl;
        cout<< setprecision(5);
        cout<< setw(15) << left << "Comp. line-k1";                 Cell(comp_a_0);     cout<<endl;
        cout<< setw(15) << left << "Comp. line-k2";                 Cell(comp_a_1);     cout<<endl;
        cout<< setw(15) << left << "Comp. k1-k2";                   Cell(comp_0_1);     cout<<endl;
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
    double comp_lines = Compatibility(samples[0].k_line, samples[1].k_line, samples[0].k_line_sig_p, samples[1].k_line_sig_p);
    cout<< setw(25) << left << "Line compatibility";   Cell(comp_lines); cout<<endl;
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++){
            double comp_means = Compatibility(samples[0].k_mean[i], samples[1].k_mean[j], samples[0].k_mean_sigma[i], samples[1].k_mean_sigma[j]);
            cout<< "Mean compatibility " << (i + 1) << " - " << (j + 1);  Cell(comp_means);  cout<<endl;
        }
    double comp_k_final = Compatibility(samples[0].k_final, samples[1].k_final, samples[0].k_final_sigma, samples[1].k_final_sigma);
    cout<< setw(25) << left << "Final k compatibility"; Cell(comp_k_final); cout<<endl;
    cout<<endl << string(72, '-') <<endl <<endl;
    vector<double> ks_line {samples[0].k_line, samples[1].k_line};
    vector<double> ks_line_sigma {samples[0].k_line_sig_p, samples[1].k_line_sig_p};
    double k_line_final = WeightedMean(ks_line, ks_line_sigma);
    double k_line_final_sigma = MeanSigma(ks_line_sigma);
    vector<double> ks_final {samples[0].k_final, samples[1].k_final};
    vector<double> ks_final_sigma {samples[0].k_final_sigma, samples[1].k_final_sigma};
    double k_final_final = WeightedMean(ks_final, ks_final_sigma);
    double k_final_final_sigma = MeanSigma(ks_final_sigma);
    cout<< setprecision(5);
        cout<< setw(15) << left << ' ';             Cell("k (m/N)");        Cell("sigma (m/N)");        cout<<endl;
        cout<< setw(15) << left << "Mean lines";    Cell(k_line_final);     Cell(k_line_final_sigma);   cout<<endl;
        cout<< setw(15) << left << "Mean finals";   Cell(k_final_final);    Cell(k_final_final_sigma);  cout<<endl;
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
    double MeanSigma = sqrt(sig) / sigmas.size();
    return MeanSigma;
}

// media aritmetica
double Mean(vector<double> data){
    double sum = 0.0;
    for (auto c : data)
        sum += c;
    double mean = sum / data.size();
    return mean;
}
// sommatoria del quadrato degli scarti
double Summation(vector<double> data){
    double summation = 0.0;
    double m = Mean(data);
    for (auto c : data)
        summation += pow((c - m), 2.0);
    return summation;
}
// formula generica della deviazione standard
double GetStdDev(vector<double> data, int data_len){
    double std_dev = sqrt(Summation(data) / data_len);
    return std_dev;
}