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

// prototipi di funzioni utilizzate nelle strutture
double M_Mean(vector<double>);
double M_Summation(vector<double>);
double M_GetStdDev(vector<double>, int);

// dichiarazione della struttura per ciascun campione
struct sample_struct {
    // variabili stabili della struttura
    string filepath;
    vector<double> data;
    double dist;
    // costruttore
    sample_struct(string path){
        filepath = path;
        ReadFile();
    }
    // costruttore per il fototraguardo iniziale (overload)
    sample_struct(double distance){
        filepath = "start";
        data = {0.0, 0.0};
        dist = distance;
    }

    // acquisizione dei dati
    void ReadFile(){
        ifstream input_file(filepath);
        if (!input_file.is_open()){
            cout<< "Error opening the file" <<endl;
            exit(1);
        }
        double value;
        input_file >> dist;
        while (input_file >> value)
            data.push_back(value);
        input_file.close();
    }

    // variabili non costanti della struttura
    // vengono calcolate quando servono senza salvarle, in modo che siano sempre aggiornate
    double Max(){
        double max = data[0];
        for (auto c : data)
            if (c > max)
                max = c;
        return max;
    }
    double Min(){
        double min = data[0];
        for (auto c : data)
            if (c < min)
                min = c;
        return min;
    }
    double Mean(){
        return M_Mean(data);
    }
    double Median(){
        vector<double> sort_data = data;
        sort(sort_data.begin(), sort_data.end());
        int half = sort_data.size() / 2;
        if(sort_data.size()%2)
            return sort_data[half + 1];
        return (sort_data[half] + sort_data[half + 1]) / 2.0;
    }

    // deviazione standard campionaria
    double StdDev(){
        return M_GetStdDev(data, data.size());
    }
    // deviazione standard singola misura
    double StdDevCorr(){
        return M_GetStdDev(data, (data.size() - 1));
    }
    // deviazione standard sulla media
    double StdDevMean(){
        return StdDevCorr()/sqrt(data.size());
    }

    // eliminazione dei dati secondo la regola del 3-sigma
    vector<double> Refine(){
        vector<double> new_data, removed_data;
        double three_sigma = 3.0 * StdDevCorr();
        double mean = Mean();
        for (auto c : data){
            if (c < mean - three_sigma || c > mean + three_sigma)
                removed_data.push_back(c);
            else new_data.push_back(c);
        }
        data = new_data;
        return removed_data;
    }

    // stampa delle informazioni sul campione sulla console
    void PrintData(){
        cout<< setprecision(5);
        cout<< "Data set size: "    << data.size()  << "\t\tDistance: "             << dist         << 'm' <<endl;
        cout<< "Minimum value: "    << Min()        << "s\t\tMaximum value: "       << Max()        << 's' <<endl;
        cout<< "Mean value: "       << Mean()       << "s\t\tMedian value: "        << Median()     << 's' <<endl;
        cout<< "Std. deviation: "   << StdDev()     << "s\tCorrected std. dev.: "   << StdDevCorr() 
        << "s\t\tMean std. dev.: "  << StdDevMean() << 's' <<endl <<endl;
        cout<< string(100, '-') <<endl <<endl;
    }
};

// struttura per le velocità di ciascun segmento (confronto tra file nella stessa cartella)
struct speed_struct{
    double speed, sp_sigma, time, tm_sigma, ds, ds_sigma, dt, dt_sigma, min, max;
    string ranges;
    // costruttore
    speed_struct(sample_struct n1, sample_struct n2, double d_sigma){
        min = n1.dist;
        max = n2.dist;
        ds = (n2.dist - n1.dist);
        ds_sigma = sqrt(pow(d_sigma, 2.0) + pow(d_sigma, 2.0));
        dt = (n2.Mean() - n1.Mean());
        dt_sigma = sqrt(pow(n2.StdDevMean(), 2.0) + pow(n1.StdDevMean(), 2.0));
        time = (n2.Mean() + n1.Mean()) / 2.0;
        tm_sigma = 0.5 * dt_sigma;
        speed = ds / dt;
        sp_sigma = abs(speed) * sqrt(pow((ds_sigma / ds), 2.0) + pow((dt_sigma / dt), 2.0));
    }
    // stampa su console
    void PrintData(){
        cout<< "Range:\t\t("        << min      << " - "        << max      << ") m"    <<endl;
        cout<< "Distance:\t"        << ds       << " m\t\t±"    << ds_sigma << " m"     <<endl;
        cout<< "Time:\t\t"          << dt       << " s\t±"      << dt_sigma << " s"     <<endl;
        cout<< "Average speed:\t"   << speed    << " m/s\t±"    << sp_sigma << " m/s"   <<endl;
        cout<< "Int. time:\t"       << time     << " s\t±"      << tm_sigma << " s"     <<endl;
        cout<< (tm_sigma / time)    << '\t'     << (sp_sigma / speed)       <<endl;
        cout<< string(100, '-') <<endl <<endl;
    }
};

// struttura per ciascuna retta di interpolazione (unione dei dati in ogni cartella)
struct interpol_struct{
    string dirname;
    double acc, acc_sigma, sp0, sp0_sigma, grav, grav_sigma, coeff_corr;
    double sp_sigma_p, acc_p, acc_sigma_p, sp0_p, sp0_sigma_p, grav_p, grav_sigma_p;
    // costruttore
    interpol_struct(vector<speed_struct> speeds, string foldername, double sin_a, double a_sigma){
        dirname = foldername;
        AccQ(speeds);
        grav = acc / sin_a;
        grav_sigma = GSigma(acc, acc_sigma, sin_a, a_sigma);
        SpSigmaPost(speeds);
        AccQ_P(speeds);
        grav_p = acc_p / sin_a;
        grav_sigma_p = GSigma(acc_p, acc_sigma_p, sin_a, a_sigma);
        CoeffCorr(speeds);
    }
    // funzioni per il calcolo delle variabili
    void AccQ(vector<speed_struct> speeds){
        double one = 0.0, xone = 0.0, xtwo = 0.0, yone = 0.0, xy = 0.0;
        for (auto c : speeds){
            double den = 1.0 / pow(c.sp_sigma, 2.0);
            one += den;
            xone += c.time * den;
            xtwo += c.time * c.time * den;
            yone += c.speed * den;
            xy += c.time * c.speed * den;
        }
        double delta = (one * xtwo) - pow(xone, 2.0);
        acc = ((one * xy) - (xone * yone)) / delta;
        sp0 =  ((xtwo * yone) - (xone * xy)) / delta;
        acc_sigma = sqrt(one / delta);
        sp0_sigma = sqrt(xtwo / delta);
    }
    double GSigma(double acc_f, double acc_s, double sin_a, double a_sigma){
        double a = pow((acc_s / sin_a), 2.0);
        double b = (1.0 - pow(sin_a, 2.0)) * pow((acc_f * a_sigma / sin_a), 2.0);
        return sqrt(a + b);
    }
    void CoeffCorr(vector<speed_struct> speeds){
        vector<double> xs, ys, xys;
        int data_len = speeds.size();
        for (auto c : speeds){
            xs.push_back(c.time);
            ys.push_back(c.speed);
            xys.push_back(c.time * c.speed);
        }
        double xs_mean = M_Mean(xs);
        double ys_mean = M_Mean(ys);
        double xys_mean = M_Mean(xys);
        double var_x = M_GetStdDev(xs, data_len);
        double var_y = M_GetStdDev(ys, data_len);
        coeff_corr = (xys_mean - xs_mean * ys_mean) / (var_x * var_y);
    }
    void SpSigmaPost(vector<speed_struct> speeds){
        double num = 0.0;
        for (auto c : speeds)
            num += pow ((c.speed - (acc * c.time) - sp0), 2.0);
        double den = speeds.size() - 2.0;
        sp_sigma_p = sqrt(num / den);
    }
    void AccQ_P(vector<speed_struct> speeds){
        double xone = 0.0, xtwo = 0.0, yone = 0.0, xy = 0.0;
        int N = speeds.size();
        for (auto c : speeds){
            xone += c.time;
            xtwo += c.time * c.time;
            yone += c.speed;
            xy += c.time * c.speed;
        }
        double delta = (N * xtwo) - (xone * xone);
        acc_p = ((N * xy) - (xone * yone)) / delta;
        sp0_p = ((xtwo * yone) - (xone * xy)) / delta;
        acc_sigma_p = sp_sigma_p * sqrt(N / delta);
        sp0_sigma_p = sp_sigma_p * sqrt(xtwo / delta);
    }
    // stampa su console
    void PrintData(){
        cout<< setprecision(4);
        cout<< "Dataset: "  << dirname  <<endl;
        cout<< "y = ax + b"             <<endl;
        cout<< "a = " << acc    << " m/s²\t±"   << acc_sigma    << " m/s²"  <<endl;
        cout<< "b = " << sp0    << " m/s\t±"    << sp0_sigma    << " m/s"   <<endl;
        cout<< "g = " << grav   << " m/s²\t±"   << grav_sigma   << " m/s²"  <<endl;
        cout<< "Correletion coefficient = "     << coeff_corr   <<endl<<endl;
        cout<< "Post. speed sigma = "           << sp_sigma_p   << " m/s"   <<endl;
        cout<< "a = " << acc_p  << " m/s²\t±"   << acc_sigma_p  << " m/s²"  <<endl;
        cout<< "b = " << sp0_p  << " m/s\t±"    << sp0_sigma_p  << " m/s"   <<endl;
        cout<< "g = " << grav_p << " m/s²\t±"   << grav_sigma_p << " m/s²"  <<endl<<endl;
        cout<< string(100, '-') <<endl <<endl;
    }
};

// prototipi di funzioni usate nel main
vector<string> GetFiles(string);
vector<sample_struct> ElaborateData(vector<string>);
vector<speed_struct> SpeedCalc(vector<sample_struct>, double);
void RemDataPrint(vector<double>);
void SpeedPrint(vector<speed_struct>);
void SpeedFileOut(vector<speed_struct>, string);
void PrintEoF();

// variabili globali per la gestione delle cartelle su cui lavora il programma
const string idir = "./data";
const string odir = "./output";

int main(){
    double d_sigma, sin_a, a_sigma;
    cout<< "Insert the standard deviation for the distances (in m): ";
    cin >> d_sigma;     // 0.0002
    cout<< "Insert the sine of the angle of inclination: ";
    cin >> sin_a;       // 0.01308959557
    cout<< "Insert the standard deviation for the angle (in rad): ";
    cin >> a_sigma;     // 0.00002908882087
    vector<string> foldernames = GetFiles(idir);        // lettura delle cartelle con i dati
    vector<interpol_struct> lines;
    for (auto c : foldernames){
        vector<string> filenames = GetFiles(c);         // lettura dei file in ciascuna cartella
        vector<sample_struct> samples = ElaborateData(filenames);
        vector<speed_struct> speeds = SpeedCalc(samples, d_sigma);
        SpeedPrint(speeds);
        SpeedFileOut(speeds, c);
        interpol_struct line = interpol_struct(speeds, c, sin_a, a_sigma);
        lines.push_back(line);
    }
    for (auto c : lines){
        c.PrintData();
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
    for (auto c : filenames)
        cout<< c <<endl;
    cout<<"Total files found: " <<filenames.size() <<endl;
    PrintEoF();
    return filenames;
}

//legge i file e salva i dati dentro le strutture
vector<sample_struct> ElaborateData(vector<string> filenames){
    vector<sample_struct> samples;
    double dist0;
    cout<< "Insert the first photogate position (in m): ";
    cin>> dist0;
    sample_struct sample0 = sample_struct(dist0);
    samples.push_back(sample0);
    PrintEoF();
    for (auto c : filenames){
        sample_struct sample = sample_struct(c);
        cout<< sample.filepath <<endl <<endl;
        vector<double> removed_data;
        do {
            removed_data = sample.Refine();
            RemDataPrint(removed_data);
        } while (removed_data.size());
        sample.PrintData();
        samples.push_back(sample);
    }
    return samples;
}

// stampa dei dati rimossi per la regola del 3-sigma
void RemDataPrint(vector<double> removed_data){
    if (removed_data.size()){
        cout<< "Removed data:" <<endl;
        cout<< setprecision(5);
        for (auto d : removed_data)
            cout<< d << '\t';
        cout<< endl;
    }
    else cout<< "All data is compatible" <<endl <<endl;
}

// crea una struttura velocità per ciascun segmento
vector<speed_struct> SpeedCalc(vector<sample_struct> samples, double d_sigma){
    vector<speed_struct> speeds;
    PrintEoF();
    for (int i = 1; i < samples.size(); i++){
        speed_struct speed = speed_struct(samples[(i - 1)], samples[i], d_sigma);
        speeds.push_back(speed);
    }
    return speeds;
}

void SpeedPrint(vector<speed_struct> speeds){
    for (auto c : speeds)
        c.PrintData();
}

// crea i file di output del programma (uno per ogni cartella)
void SpeedFileOut(vector<speed_struct> speeds, string foldername){
    istringstream folderin(foldername);
    string filename;
    while (getline(folderin, filename, '/'));
    string filepath = odir + "/" + filename + ".txt";
    ofstream ofile(filepath);
    if (!ofile.is_open()){
        cout<< "Error: can't write the file" <<endl;
        return;
    }
    for (auto c : speeds)
        ofile<< c.time << '\t' << c.tm_sigma << '\t' << c.speed << '\t' << c.sp_sigma <<endl;
    ofile.close();
    cout<< "File successfully created" <<endl;
    PrintEoF();
}

void PrintEoF(){
    cout<<endl << string(100, '=') <<endl <<endl;
}

// funzioni utilizzate dentro le strutture

double M_Mean(vector<double> data){
    double sum = 0.0;
    for (auto c : data)
        sum += c;
    double mean = sum / data.size();
    return mean;
}
// sommatoria del quadrato degli scarti
double M_Summation(vector<double> data){
    double summation = 0.0;
    double m = M_Mean(data);
    for (auto c : data)
        summation += pow((c - m), 2.0);
    return summation;
}
// formula generica della deviazione standard
double M_GetStdDev(vector<double> data, int data_len){
    double std_dev = sqrt(M_Summation(data) / data_len);
    return std_dev;
}