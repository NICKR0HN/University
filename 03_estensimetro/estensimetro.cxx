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
const double x_factor = 4.0 * 9.806 * 0.001;
const double y_factor = 0.000001;

void PrintEoF(){
    cout<<endl << string(40, '=') <<endl <<endl;
}

struct sample_struct{
    string filename;
    vector<double> forces, lengths;
    double force_sigma, length_sigma;
    double k_line, q_line, k_line_sig, q_line_sig;
    vector<double> k_vector, k_sigma_vector;
    double k_mean, k_mean_sigma;
    sample_struct(string filepath){
        filename = filepath;
        ReadFile();
        Interpol();
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
        force_sigma = force_sigma_r * x_factor;
        input_file >> length_sigma_r;
        length_sigma = length_sigma_r * y_factor;
        double force_r, force, length_r, length;
        while ((input_file >> force_r) && (input_file >> length_r)){
            force = force_r * x_factor;
            forces.push_back(force);
            length = length_r * y_factor;
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
        double num = 0.0, den = 0.0, sig = 0.0;
        for (int i = 0; i < k_size; i++){
            double sq_sig = pow(k_sigma_vector[i], 2.0);
            num += k_vector[i] / sq_sig;
            den += 1.0 / sq_sig;
            sig += sq_sig;
        }
        k_mean = num / den;
        k_mean_sigma = sqrt(sig) / k_size;
    }

    void PrintData(){
        cout<< filename <<endl<<endl;
        PrintInterpol();
        cout<< string(40, '-') <<endl <<endl;
        PrintKVector();
        PrintEoF();
    }
    void PrintInterpol(){
        cout<< setprecision(4);
        cout<< "y = kx + q"             <<endl;
        cout<< "k = " << k_line     << " m/N\t±"    << k_line_sig   << " m/N"   <<endl;
        cout<< "q = " << q_line     << " m\t±"      << q_line_sig   << " m"     <<endl;
    }
    void PrintKVector(){
        cout<< setprecision(4);
        cout<< "Mean k = " << k_mean    << "m/N\t±" << k_mean_sigma << " m/N"   <<endl;
    }
};

//prototipi funzioni
vector<string> GetFiles(string);
vector<sample_struct> ElaborateData(vector<string>);

int main(){

    vector<string> foldernames = GetFiles(idir);        // lettura delle cartelle con i dati
    for (string folder : foldernames){
        vector<string> filenames = GetFiles(folder);
        vector<sample_struct> samples = ElaborateData(filenames);
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