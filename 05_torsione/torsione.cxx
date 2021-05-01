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
const string ifile = "./data/1_period.txt";
const string ofile = "./data/1_results.txt";
const double dt = 0.02;

// stampa su console
template <typename T> void Cell(T text){
    cout<< setw(15) << right << text;
}
void PrintEoF(char sep){
    cout<<endl << string(75, sep) <<endl <<endl;
}

// prototipi funzioni matematiche
double Min(vector<double>);
double Max(vector<double>);
double Mean(vector<double>);
double WeightedMean(vector<double>, vector<double>);
double MeanSigma(vector<double>);
double Summation(vector<double>);
double GetStdDev(vector<double>, int);
double CoeffCorr(vector<double>, vector<double>);
void ThreeSigma(vector<double>&, vector<double>&);
void DelFirstLast(vector<double>&);

// struttura dei dati di una singola sinusoide
struct dataset_struct {
    double per_real, limit;
    vector<double> data;
    vector<double> times, times_s, periods, per_sigmas, correlations, per_counter;
    double per_avg, per_stddev, per_sigma, puls, puls_sigma;
    vector<double> tops, tops_s, lows, lows_s, max_counter;
    double top, top_stddev, top_sigma;
    double low, low_stddev, low_sigma;
    double amp, amp_sigma;
    dataset_struct(){};
    dataset_struct(double period_in, vector<double> data_in){
        per_real = period_in;
        data = data_in;
        GetLimit();
        GetIntersections();
        GetPeriods();
        GetMaxes();
        GetAmplitude();
        PrintData();
    }

    void GetLimit(){
        vector<double> squared;
        for(double c : data)
            squared.push_back(c * c);
        double mean = Mean(squared);
        limit = sqrt(mean) / 1.5;
    }

    // intersezioni con l'asse x
    void GetIntersections(){
        int i = 0;
        while(i < data.size()){
            vector<double> ys, xs;
            while((i < data.size()) && (abs(data[i]) < limit)){
                ys.push_back(data[i]);
                xs.push_back(i * dt);
                i++;
            }
            if (ys.size() > 1){
                array<double, 2> line = Interpol(xs, ys);
                double time = line[0];
                times.push_back(time);
                double sigma = line[1];
                times_s.push_back(sigma);
                double corr = CoeffCorr(xs, ys);
                correlations.push_back(abs(corr));
                int count = ys.size();
                per_counter.push_back(count);
            }
            i++;
        }
        DelFirstLast(times);
        DelFirstLast(correlations);
        DelFirstLast(per_counter);
        DelFirstLast(times_s);
    }

    void GetPeriods(){
        for (int i = 1; i < times.size(); i = i+2){
            periods.push_back(times[i] - times[i-1]);
            double sigma = sqrt(pow(times_s[i], 2.0) + pow(times_s[i-1], 2.0));
            per_sigmas.push_back(sigma);
        }
        ThreeSigma(periods, per_sigmas);
        per_avg = 2.0 * WeightedMean(periods, per_sigmas);
        per_stddev = 2.0 * GetStdDev(periods, (periods.size() - 1));
        per_sigma = 2.0 * MeanSigma(per_sigmas);
        puls = 2.0 * M_PI / per_avg;
        puls_sigma = 2.0 * M_PI * per_sigma / pow(puls, 2.0);
    }

    // massimi e minimi
    void GetMaxes(){
        int i = 1;
        while(i < data.size()){
            vector<double> xs, ys;
            while((i < data.size()) && ((data[i-1] * data[i]) > 0)){
                if(abs(data[i]) > limit){
                    xs.push_back(i * dt);
                    ys.push_back(data[i]);
                }
                i++;
            }
            if(ys.size() > 3){
                array<double, 2> parabola = Parabola(xs, ys);
                double max = parabola[0];
                double sigma = parabola[1];
                int count = ys.size();
                if(max > 0){
                    tops.push_back(max);
                    tops_s.push_back(sigma);
                }
                else{
                    lows.push_back(max);
                    lows_s.push_back(sigma);
                }
                max_counter.push_back(count);
            }
            i++;
        }
        DelFirstLast(tops);
        DelFirstLast(lows);
        DelFirstLast(max_counter);
        DelFirstLast(tops_s);
        DelFirstLast(lows_s);
    }

    void GetAmplitude(){
        ThreeSigma(tops, tops_s);
        ThreeSigma(lows, lows_s);
        top = WeightedMean(tops, tops_s);
        low = WeightedMean(lows, lows_s);
        top_stddev = GetStdDev(tops, tops.size() - 1);
        low_stddev = GetStdDev(lows, lows.size() - 1);
        top_sigma = MeanSigma(tops_s);
        low_sigma = MeanSigma(lows_s);
        vector<double> vals = {top, abs(low)}, sigs = {top_sigma, low_sigma};
        amp = WeightedMean(vals, sigs);
        amp_sigma = MeanSigma(sigs);
    }

    void PrintData(){
        cout<< "Expected period = " << per_real <<endl;
        cout<< "Limit = " << limit <<endl;
        cout<< "Intersections = " << times.size() <<endl;
        cout<< "Periods = " << periods.size() <<endl;
        cout<< "Tops = " << tops.size() <<endl;
        cout<< "Lows = " << lows.size() <<endl;
        Cell(' '); Cell("Average"); Cell("Minimum"); Cell("Maximum"); Cell("Std. dev."); cout<<endl;
        cout<< setw(15) << left << "Period"; Cell(per_avg); Cell(2.0 * Min(periods)); Cell(2.0 * Max(periods)); Cell(per_stddev); cout<<endl;
        cout<< setw(15) << left << "Correlation"; Cell(Mean(correlations)); Cell(Min(correlations)); Cell(Max(correlations)); cout<<endl;
        cout<< setw(15) << left << "Number points"; Cell(Mean(per_counter)); Cell(Min(per_counter)); Cell(Max(per_counter)); cout<<endl;
        cout<< setw(15) << left << "Top amplitude"; Cell(top); Cell(Min(tops)); Cell(Max(tops)); Cell(top_stddev); cout<<endl;
        cout<< setw(15) << left << "Low amplitude"; Cell(low); Cell(Max(lows)); Cell(Min(lows)); Cell(low_stddev); cout<<endl;
        cout<< setw(15) << left << "Number points"; Cell(Mean(max_counter)); Cell(Min(max_counter)); Cell(Max(max_counter)); cout<<endl;
    }

    array<double, 2> Interpol(vector<double> xs, vector<double> ys){
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
        double num = 0.0;
        for (int i = 0; i < data_size; i++)
            num += pow ((ys[i] - (coeff * xs[i]) - rcept), 2.0);
        double sigma_post = sqrt(num / (data_size - 2.0));
        array<double, 2> out = {interpol, sigma_post};
        return out;
    }

    // interpolazione parabolica
    array<double, 2> Parabola(vector<double> xs, vector<double> ys){
        double x1 = 0.0, x2 = 0.0, x3 = 0.0, x4 = 0.0, y1 = 0.0, x1y1 = 0.0, x2y1 = 0.0;
        int n = xs.size();
        for(int i = 0; i < n; i++){
            x1 += xs[i];
            x2 += pow(xs[i], 2.0);          //  [x4] [x3] [x2] | [x2y1]
            x3 += pow(xs[i], 3.0);          //  [x3] [x2] [x1] | [x1y1]
            x4 += pow(xs[i], 4.0);          //  [x2] [x1] [n]  |  [y1]
            y1 += ys[i];
            x1y1 += xs[i] * ys[i];
            x2y1 += pow(xs[i], 2.0) * ys[i];
        }
        // determinante delle sottomatrici
        double lt = (n * x2) - (x1 * x1);   //  [lt] [mt] [rt]
        double lm = (n * x3) - (x1 * x2);   //  [lm] [mm] [rm]
        double lb = (x1 * x3) - (x2 * x2);  //  [lb] [mb] [rb]
        double mt = lm;
        double mm = (n * x4) - (x2 * x2);
        double mb = (x1 * x4) - (x2 * x3);
        double rt = lb;
        double rm = mb;
        double rb = (x2 * x4) - (x3 * x3);
        // metodo di Cramer
        double delta = (x4 * lt) - (x3 * lm) + (x2 * lb);
        // sviluppo di LaPlace sulle colonne della variabile corrispondente
        double a = ((x2y1 * lt) - (x1y1 * lm) + (y1 * lb)) / delta;
        double b = (- (x2y1 * mt) + (x1y1 * mm) - (y1 * mb)) / delta;
        double c = ((x2y1 * rt) - (x1y1 * rm) + (y1 * rb)) / delta;
        double vertix = (c - ((b * b) / (4.0 * a)));
        double num = 0.0;
        for (int i = 0; i < n; i++){
            num += pow(ys[i] - (a * pow(xs[i], 2.0) + b * xs[i] + c), 2.0);
        }
        double sigma_post = sqrt(num / (n - 3));
        array<double, 2> out = {vertix, sigma_post};
        return out;
    }
};

// struttura per ciascuna frequenza
struct sample_struct {
    double freq;
    dataset_struct forces, angles;
    sample_struct(string forces_in, string angles_in){
        ReadLine(forces_in, angles_in);
        WriteFile();
    }

    void ReadLine(string forces_in, string angles_in){
        istringstream forces_line(forces_in), angles_line(angles_in);
        forces_line >> freq;
        cout<< "FREQUENCY = " << freq << "Hz" <<endl;
        PrintEoF('-');
        double temp;
        angles_line >> temp;
        if(freq != temp){
            cout << "Missing line" <<endl;
            exit(1);
        }
        double real_p = 1.0 / freq;
        double force, angle;
        vector<double> forces_data, angles_data;
        while ((forces_line >> force) && (angles_line >> angle)){
            if ((force == 0) && (!forces_data.empty()) && (forces_data.back() == 0)){
                cout<< "Motor stopped" <<endl;
                break;
            }
            forces_data.push_back(force * 2.0 * M_PI);
            angles_data.push_back(angle * 2.0 * M_PI);
        }
        cout<< "FORCING" <<endl<<endl;
        forces = dataset_struct(real_p, forces_data);
        PrintEoF('/');
        cout<< "PENDULUM" <<endl<<endl;
        angles = dataset_struct(real_p, angles_data);
        PrintEoF('#');
    }

    void WriteFile(){
        ofstream write(ofile, ofstream::app);
        if (!write.is_open()){
            cout<< "Permission denied" <<endl;
            exit(1);
        }
        write<< freq <<'\t' << forces.per_avg <<'\t' << forces.per_sigma <<'\t' << forces.amp <<'\t' << forces.amp_sigma <<'\t';
        write<< angles.per_avg <<'\t' << angles.per_sigma <<'\t' << angles.top <<'\t' << angles.top_sigma <<'\t'
             << angles.low <<'\t' << angles.low_sigma <<'\t' << angles.amp <<'\t' << angles.amp_sigma <<'\t' << angles.puls <<'\t' <<angles.puls_sigma <<endl;
        write.close();
    }
};

// prototipi funzioni main
void OverWriteFile();
vector<sample_struct> ReadFile();

int main(){
    OverWriteFile();
    vector<sample_struct> samples = ReadFile();
    return 0;
}

// funzioni main
void OverWriteFile(){
    ofstream overwrite(ofile, ofstream::trunc);
    if (!overwrite.is_open()){
        cout<< "Permission denied" <<endl;
        exit(1);
    }
    overwrite.close();
}

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
double Min(vector<double> data){
    double min = data[0];
        for (auto c : data)
            if (c < min)
                min = c;
        return min;
};

double Max(vector<double> data){
    double max = data[0];
        for (double c : data)
            if (c > max)
                max = c;
        return max;
};

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

void DelFirstLast(vector<double> &data){
    data.erase(data.begin());
    data.pop_back();
}