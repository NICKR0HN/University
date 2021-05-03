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
const string ifile = "./data/2_smorzato.txt";
const string ofile = "./data/2_results.txt";
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
template <typename T> void DelFirstLast(vector<T>&);

struct curve_struct{
    double limit;
    double rcept, rcept_s, x_max, y_max, y_max_s, x_min, y_min, y_min_s;
    int line_points, parabola_points;
    curve_struct(vector<double> xs, vector<double> ys){
        GetLimit(ys);
        Intercept(xs, ys);
        GetMaxMin(xs, ys);
    }
    
    void GetLimit(vector<double> data){
        vector<double> squared;
        for(double c : data)
            squared.push_back(c * c);
        double mean = Mean(squared);
        limit = sqrt(mean) / 1.5;
    }

    void Intercept(vector<double> xs, vector<double> ys){
        vector<double> xs_line, ys_line;
        int i = 0;
        double max_temp = Max(ys);
        double min_temp = Min(ys);
        while (ys[i] != max_temp)
            i++;
        while ((ys[i] != min_temp) && (i < xs.size())){
            if(abs(ys[i]) < limit){
                xs_line.push_back(xs[i]);
                ys_line.push_back(ys[i]);
            }
            i++;
        }
        array<double, 2> line = Interpol(xs_line, ys_line);
        rcept = line[0];
        rcept_s = line[1];
        line_points = xs_line.size();
    }

    void GetMaxMin(vector<double> xs, vector<double> ys){
        vector<double> xs_par, ys_par;
        int i = 0;
        while ((ys[i] >= 0) && (i < xs.size())){
            if (ys[i] > limit){
                xs_par.push_back(xs[i]);
                ys_par.push_back(ys[i]);
            }
            i++;
        }
        array<double, 3> parabola = Parabola(xs_par, ys_par);
        x_max = parabola[0];
        y_max = parabola[1];
        y_max_s = parabola[2];
        parabola_points = xs_par.size();
        xs_par.clear();
        ys_par.clear();
        while ((ys[i] <= 0) && (i < xs.size())){
            if (ys[i] < limit){
                xs_par.push_back(xs[i]);
                ys_par.push_back(ys[i]);
            }
            i++;
        }
        parabola = Parabola(xs_par, ys_par);
        x_min = parabola[0];
        y_min = parabola[1];
        y_min_s = parabola[2];
        if(parabola_points > xs_par.size())
            parabola_points = xs_par.size();
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
        double x_sigma = sigma_post / coeff;
        array<double, 2> out = {interpol, x_sigma};
        return out;
    }

    // interpolazione parabolica
    array<double, 3> Parabola(vector<double> xs, vector<double> ys){
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
        double axis = -0.5 * (b / a);
        double vertix = (c - ((b * b) / (4.0 * a)));
        double num = 0.0;
        for (int i = 0; i < n; i++){
            num += pow(ys[i] - (a * pow(xs[i], 2.0) + b * xs[i] + c), 2.0);
        }
        double sigma_post = sqrt(num / (n - 3));
        array<double, 3> out = {axis, vertix, sigma_post};
        return out;
    }
};

// struttura dei dati di una singola sinusoide
struct sample_struct {
    vector<double> data;
    vector<curve_struct> curves;
    double per_avg, per_avg_s, pulse, pulse_s;
    vector<double> periods, periods_s, periods_points;
    double max_a, max_a_s, max_b, max_b_s, max_corr;
    double min_a, min_a_s, min_b, min_b_s, min_corr;
    vector<double> amp_points;
    sample_struct(string data_in){
        istringstream data_line(data_in);
        double angle;
        while ((data_line >> angle))
            data.push_back(angle * 2.0 * M_PI);
        GetCurves();
        GetPeriods();
        GetLines();
        PrintData();
        WriteFile();
    }

    void GetCurves(){
        int n = data.size();
        int i = 0;
        while (i < n){
            vector<double> xs, ys;
            while ((i < n) && (data[i] >= 0)){
                xs.push_back(i * dt);
                ys.push_back(data[i]);
                i++;
            }
            while((i < n) && (data[i] < 0)){
                xs.push_back(i * dt);
                ys.push_back(data[i]);
                i++;
            }
            curve_struct curve = curve_struct(xs, ys);
            curves.push_back(curve);
        }
        DelFirstLast(curves);
    }

    void GetPeriods(){
        for (int i = 1; i < curves.size(); i = i+2){
            double period = curves[i].rcept - curves[i-1].rcept;
            periods.push_back(period);
            double period_s = sqrt(pow(curves[i].rcept_s, 2.0) + pow(curves[i-1].rcept_s, 2.0));
            periods_s.push_back(period_s);
        }
        for (curve_struct c : curves)
            periods_points.push_back(c.line_points);
        ThreeSigma(periods, periods_s);
        per_avg = WeightedMean(periods, periods_s);
        per_avg_s = MeanSigma(periods_s);
        pulse = 2.0 * M_PI / per_avg;
        pulse_s = 2.0 * M_PI * per_avg_s / pow(per_avg, 2.0);
    }

    void GetLines(){
        vector<double> x_maxs, y_maxs, y_maxs_s;
        vector<double> x_mins, y_mins, y_mins_s;
        for (curve_struct c : curves){
            x_maxs.push_back(c.x_max);
            y_maxs.push_back(log(c.y_max));
            y_maxs_s.push_back(c.y_max_s / c.y_max);
            x_mins.push_back(c.x_min);
            y_mins.push_back(log(abs(c.y_min)));
            y_mins_s.push_back(abs(c.y_min_s / c.y_min));
            amp_points.push_back(c.parabola_points);
        }
        max_corr = CoeffCorr(x_maxs, y_maxs);
        array<double, 4> line = InterpolS(x_maxs, y_maxs, y_maxs_s);
        max_a = line[0];
        max_a_s = line[1];
        max_b = line[2];
        max_b_s = line[3];
        min_corr = CoeffCorr(x_mins, y_mins);
        line = InterpolS(x_mins, y_mins, y_mins_s);
        min_a = line[0];
        min_a_s = line[1];
        min_b = line[2];
        min_b_s = line[3];
        
    }

    array<double, 4> InterpolS(vector<double> xs, vector<double> ys, vector<double> ys_s){
        double one = 0.0, xone = 0.0, xtwo = 0.0, yone = 0.0, xy = 0.0;
        for (int i = 0; i < xs.size(); i++){
            double den = 1.0 / pow(ys_s[i], 2.0);
            one += den;
            xone += xs[i] * den;
            xtwo += pow(xs[i], 2.0) * den;
            yone += ys[i] * den;
            xy += xs[i] * ys[i] * den;
        }
        double delta = (one * xtwo) - pow(xone, 2.0);
        double a = ((one * xy) - (xone * yone)) / delta;
        double b =  ((xtwo * yone) - (xone * xy)) / delta;
        double a_s = sqrt(one / delta);
        double b_s = sqrt(xtwo / delta);
        array<double, 4> out = {a, a_s, b, b_s};
        return out;
    }

    void PrintData(){
        cout<< "Period = " << per_avg << " s \t +- " << per_avg_s << " s \t (Std. Dev. = " << GetStdDev(periods, periods.size()-1) << ')' <<endl;
        cout<< " Using " << Mean(periods_points) << " points (from " << Min(periods_points) << " to " << Max(periods_points) << ')' <<endl;
        cout<< " Using " << periods.size() << " periods" <<endl;
        cout<< "Pulse = " << pulse << " Hz \t +- " << pulse_s << " Hz" <<endl;
        cout<< "Amplitudes: using " << Mean(amp_points) << " points (from " << Min(amp_points) << " to " << Max(amp_points) << ')' <<endl;
        cout<< "Maximums logarithmic correlation = " << max_corr <<endl;
        cout<< "a  = " << max_a << "\t +- " << max_a_s <<endl;
        cout<< "b  = " << max_b << "\t +- " << max_b_s <<endl;
        cout<< "Minimums logarithmic correlation = " << min_corr <<endl;
        cout<< "a  = " << min_a << "\t +- " << min_a_s <<endl;
        cout<< "b  = " << min_b << "\t +- " << min_b_s <<endl;
    }

    void WriteFile(){
        ofstream fileout(ofile, ofstream::app);
        fileout<< per_avg <<'\t' << per_avg_s <<'\t' << GetStdDev(periods, periods.size()-1) <<'\t' << Mean(periods_points) <<'\t' << pulse <<'\t'
            << pulse_s <<'\t' << Mean(amp_points) <<'\t' << max_corr <<'\t' << max_a <<'\t' << max_a_s <<'\t' << max_b <<'\t' << max_b_s
            <<'\t' << min_corr <<'\t' << min_a <<'\t' << min_a_s <<'\t' << min_b <<'\t' << min_b_s <<endl;
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
    string data;
    int i = 1;
    while(getline(input_file, data)){
        PrintEoF('#');
        cout<< "Line " << i <<endl;
        PrintEoF('-');
        sample_struct sample = sample_struct(data);
        samples.push_back(sample);
        i++;
    }
    PrintEoF('#');
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

template <typename T> void DelFirstLast(vector<T> &data){
    data.erase(data.begin());
    data.pop_back();
}