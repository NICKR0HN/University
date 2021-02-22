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
        double sum = 0.0;
        for (auto c : data)
            sum += c;
        double mean = sum / data.size();
        return mean;
    }
    double Median(){
        vector<double> sort_data = data;
        sort(sort_data.begin(), sort_data.end());
        int half = sort_data.size() / 2;
        if(sort_data.size()%2)
            return sort_data[half + 1];
        return (sort_data[half] + sort_data[half + 1]) / 2.0;
    }

    // sommatoria del quadrato degli scarti
    double Summation(){
        double summation = 0.0;
        double m = Mean();
        for (auto c : data)
            summation += pow((c - m), 2.0);
        return summation;
    }
    // formula generica della deviazione standard
    double GetStdDev(int data_len){
        double std_dev = sqrt(Summation() / data_len);
        return std_dev;
    }
    // deviazione standard campionaria
    double StdDev(){
        double std_dev = GetStdDev(data.size());
        return std_dev;
    }
    // deviazione standard singola misura
    double StdDevCorr(){
        double std_dev_corr = GetStdDev(data.size() - 1);
        return std_dev_corr;
    }
    // deviazione standard sulla media
    double StdDevMean(){
        double std_dev_mean = StdDevCorr()/sqrt(data.size());
        return std_dev_mean;
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
            << "s\t\tMean std. dev.: " << StdDevMean() << 's' <<endl <<endl;
        cout<< string(100, '-') <<endl <<endl;
    }
};

//struttura per immagazzinare i dati relativi alle velocità dei singoli segmenti
struct speed_struct{
    double speed, sp_sigma, time, tm_sigma, ds, ds_sigma, dt, dt_sigma, min, max;
    string ranges;
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
    void PrintData(){
        cout<< "Range:\t\t("        << min      << " - "        << max      << ") m"    <<endl;
        cout<< "Distance:\t"        << ds       << " m\t\t±"    << ds_sigma << " m"     <<endl;
        cout<< "Time:\t\t"          << dt       << " s\t±"      << dt_sigma << " s"     <<endl;
        cout<< "Average speed:\t"   << speed    << " m/s\t±"    << sp_sigma << " m/s"   <<endl;
        cout<< "Int. time:\t"       << time     << " s\t±"      << tm_sigma << " s"     <<endl <<endl;
        cout<< string(100, '-') <<endl <<endl;
    }
};

struct interpol_struct{
    string dirname;
    double acc = 0.0, acc_sigma, q = 0.0, q_sigma;
    interpol_struct(vector<speed_struct> speeds, string foldername){
        dirname = foldername;
        double h = 0.0, k = 0.0, n = speeds.size();
        for (auto c : speeds){
            h += pow(c.time, 2.0);
            k += c.time;
        }
        double delta = n * h + k;
        for (auto c : speeds){
            acc += (n * c.time - k) * c.speed / delta;
        }
        // formula acc_sigma
        for (auto c : speeds){
            q += (h - (k * c.time)) * c.speed / delta;
        }
        // formula q_sigma
    }
    void PrintData(){
        cout<< "Dataset: "  << dirname  <<endl;
        cout<< "y = ax + b"             <<endl;
        cout<< "a = " << acc    << " m/s²\t±"   << acc_sigma    << " m/s²"  <<endl;
        cout<< "b = " << q      << " m/s\t±"    << q_sigma      << " m/s"   <<endl <<endl;
        cout<< string(100, '-') <<endl <<endl;
    }
};

vector<string> GetFiles(string);
vector<sample_struct> ElaborateData(vector<string>);
vector<speed_struct> SpeedCalc(vector<sample_struct>);
void RemDataPrint(vector<double>);
void SpeedPrint(vector<speed_struct>);
void SpeedFileOut(vector<speed_struct>, string);
void PrintEoF();

const string idir = "./data";
const string odir = "./output";

int main(){
    vector<string> foldernames = GetFiles(idir);
    vector<interpol_struct> lines;
    for (auto c : foldernames){
        vector<string> filenames = GetFiles(c);
        vector<sample_struct> samples = ElaborateData(filenames);
        vector<speed_struct> speeds = SpeedCalc(samples);
        // SpeedPrint(speeds);
        SpeedFileOut(speeds, c);
        interpol_struct line = interpol_struct(speeds, c);
        lines.push_back(line);
    }
    for (auto c : lines){
        c.PrintData();
    }
    return 0;
}

//ottiene i nomi di tutti i file presenti nella cartella "data"
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
    filenames.pop_back();
    filenames.pop_back();
    sort(filenames.begin(), filenames.end());
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
            // RemDataPrint(removed_data);
        } while (removed_data.size());
        // sample.PrintData();
        samples.push_back(sample);
    }
    return samples;
}

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

//crea una struttura velocità per ciascun segmento
vector<speed_struct> SpeedCalc(vector<sample_struct> samples){
    vector<speed_struct> speeds;
    double d_sigma;
    cout<< "Insert the standard deviation for the distances (in m): ";
    cin>> d_sigma;
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