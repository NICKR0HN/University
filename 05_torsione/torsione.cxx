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
    cout<<endl << string(60, sep) <<endl <<endl;
}

// prototipi funzioni matematiche
double Mean(vector<double>);
double Summation(vector<double>);
double GetStdDev(vector<double>, int);
double CoeffCorr(vector<double>, vector<double>);

// struttura set di data
struct sample_struct {
    double freq;
    vector<double> forces, angles;
    vector<double> f_periods, a_periods;
    sample_struct(string forces_in, string angles_in){
        istringstream forces_line(forces_in), angles_line(angles_in);
        forces_line >> freq;
        cout<< "Frequency = " << freq << "Hz" <<endl<<endl;
        double temp;
        angles_line >> temp;
        if(freq != temp){
            cout << "Missing line" <<endl;
            exit(1);
        }
        double force, angle;
        while ((forces_line >> force) && (angles_line >> angle)){
            if ((force == 0) && (!forces.empty()) && (forces.back() == 0)){
                cout<< "Motor stopped" <<endl;
                break;
            }
            forces.push_back(force * 2.0 * M_PI);
            angles.push_back(angle * 2.0 * M_PI);
        }
        f_periods = Periods(forces);
        a_periods = Periods(angles);
        PrintEoF('#');
    }

    double GetLimit(vector<double> data){
        vector<double> squared;
        for(double c : data)
            squared.push_back(c * c);
        double mean = Mean(squared);
        double limit = sqrt(mean) / 2.0;
        return limit;
    }

    vector<double> Periods(vector<double> data){
        double limit = GetLimit(data);
        vector<double> times, periods;
        cout<< "Periods" <<endl; Cell("Period"); Cell("Corr."); Cell("Length"); Cell("Position"); cout<<endl;
        int i = 0;
        while(i < data.size()){
            vector<double> ys, xs;
            while((i < data.size()) && (abs(data[i]) < limit)){
                ys.push_back(data[i]);
                xs.push_back(i * dt);
                i++;
            }
            i++;
            if (ys.size() > 1){
                double time = Interpol(xs, ys);
                double period = 0.0;
                if (!times.empty()){
                    period = time - times.back();
                    periods.push_back(period);
                }
                times.push_back(time);
                double corr = CoeffCorr(xs, ys);
                int length = ys.size();
                Cell(period); Cell(corr); Cell(length); Cell(i); cout<<endl;
            }
        }
        PrintEoF('-');
        return periods;
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
    string forces, angles;
    while(getline(input_file, forces) && getline(input_file, angles)){
        sample_struct sample = sample_struct(forces, angles);
        samples.push_back(sample);
    }
    return samples;
}

// funzioni matematiche
double Mean(vector<double> data){
    double sum = 0.0;
    for (double c : data)
        sum += c;
    double mean = sum / data.size();
    return mean;
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