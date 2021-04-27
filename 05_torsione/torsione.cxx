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
const string ofile = "./data/results.txt";
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
vector<double> ThreeSigma(vector<double>);

// struttura set di data
struct sample_struct {
    double freq, real_p;
    vector<double> forces, angles;
    vector<double> f_periods, a_periods;
    array<double, 2> f_avg_p, a_avg_p;
    sample_struct(string forces_in, string angles_in){
        ReadData(forces_in, angles_in);
        real_p = 1.0 / freq;
        cout<< "Forcing" <<endl;
        f_periods = Periods(forces);
        f_avg_p = AvgPeriod(f_periods);
        PrintEoF('-');
        cout<< "Pendulum" <<endl;
        a_periods = Periods(angles);
        a_avg_p = AvgPeriod(a_periods);
        PrintEoF('#');
        WriteFile();
    }

    void ReadData(string forces_in, string angles_in){
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
    }

    double GetLimit(vector<double> data){
        vector<double> squared;
        for(double c : data)
            squared.push_back(c * c);
        double mean = Mean(squared);
        double limit = sqrt(mean) / 1.5;
        return limit;
    }

    vector<double> GetIntervals(vector<double> data){
        double limit = GetLimit(data);
        vector<double> times, intervals;
        // Cell("Interval"); Cell("Corr."); Cell("Length"); Cell("Position"); cout<<endl;
        int i = 0;
        while(i < data.size()){
            vector<double> ys, xs;
            while((i < data.size()) && (abs(data[i]) < limit)){
                ys.push_back(data[i]);
                xs.push_back(i * dt);
                i++;
            }
            if (ys.size() > 1){
                double time = Interpol(xs, ys);
                double interval = 0.0;
                if (!times.empty()){
                    interval = time - times.back();
                    intervals.push_back(interval);
                }
                times.push_back(time);
                double corr = CoeffCorr(xs, ys);
                int length = ys.size();
                // Cell(interval); Cell(corr); Cell(length); Cell(i); cout<<endl;
            }
            i++;
        }
        cout<< "Intersections found: " << intervals.size() <<endl;
        return intervals;
    }

    vector<double> Periods(vector<double> data){
        vector<double> intervals = GetIntervals(data);
        vector<double> periods;
        for (int i = 0; i < intervals.size(); i = i+2)
            periods.push_back(intervals[i]);
        periods = ThreeSigma(periods);
        cout<< "Refined: " << periods.size() <<endl;
        return periods;
    }

    array<double, 2> AvgPeriod(vector<double> periods){
        double avg_p = 2.0 * Mean(periods);
        double std_dev = 2.0 * GetStdDev(periods, (periods.size() - 1));
        cout<< "Exp. period = \t"   << real_p   <<endl;
        cout<< "Avg. period = \t"   << avg_p    <<endl;
        cout<< "Std. deviation=\t"  << std_dev  <<endl;
        array<double, 2> out = {avg_p, std_dev};
        return out;
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
        double interpol = -1.0 * rcept / coeff;
        return interpol;
    }

    void WriteFile(){
        ofstream output(ofile, ofstream::app);
        if (!output.is_open())
            exit(1);
        output<< freq <<'\t' << f_avg_p[0] <<'\t' << f_avg_p[1] <<'\t' << a_avg_p[0] <<'\t' << a_avg_p[1] <<endl;
        output.close();
    }
};

// prototipi funzioni main
vector<sample_struct> ReadFile();

int main(){
    ofstream overwrite(ofile, ofstream::trunc);
    if (!overwrite.is_open()){
        cout<< "Permission denied" <<endl;
        exit(1);
    }
    overwrite.close();
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

vector<double> ThreeSigma(vector<double> data){
    vector<double> new_data;
    while (new_data.size() < data.size()){
        if (!new_data.empty()){
            data = new_data;
            new_data.clear();
        }
        double three_sigma = 3.0 * GetStdDev(data, (data.size() - 1));
        double mean = Mean(data);
        for (double c : data)
            if(((mean - three_sigma) < c) && (c < (mean + three_sigma)))
                new_data.push_back(c);
    }
    return new_data;
}