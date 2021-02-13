#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <algorithm> 
using namespace std;

//dichiarazione struttura delle classi dell'istogramma
struct class_struct{
    double min, max, centroid;
    int freq = 0;
    double gauss;
    //costruttore della struttura, usata quando viene inizializzata una nuova struttura
    class_struct(double minimum, double maximum, vector<double> data, double mean, double sigma, double delta, int n_classes){
        min = minimum;
        max = maximum;
        centroid = (max + min) / 2.0;
        for(auto c : data)
            if((min <= c) && (c < max))
                freq++;
        double den = sqrt(2.0 * M_PI) * sigma;
        double gaussian = exp(Exponent(mean, sigma)) / den;
        gauss = gaussian * delta * n_classes;
    }
    double Exponent(double mean, double sigma){
        double num = pow((centroid - mean), 2.0);
        double den = 2.0 * pow(sigma, 2.0);
        double exponent = (-1.0) * (num / den);
        return exponent;
    }
};

// dichiarazione della struttura per ciascun campione
struct sample_struct {
    // variabili stabili della struttura
    string filename;
    vector<double> data;
    int n_classes = 20;
    vector<class_struct> classes;

    // costruttore
    sample_struct(string filepath){
        filename = filepath;
    }

    // acquisizione dei dati
    void ReadFile(){
        ifstream input_file(filename);
        if (!input_file.is_open())
            throw;
        double value;
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
        return (sort_data[half] + sort_data[half + 1]) / 2;
    }
    // ampiezza delle classi
    double Delta(){
        double delta = (Max() - Min()) / n_classes;
        return delta;
    }
    // sommatoria del quadrato degli scarti
    double Summation(){
        double summation = 0.0;
        double m = Mean();
        for(auto c : data)
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

    // calcolo della deviazione standard normalizzata
    // viene attivata solo se è indicato 'true' in FileAnalysis() nel main()
    // sostituisce i dati con i relativi scarti normalizzati
    void NormDev(){
        vector<double> old_data = data;
        double sigma = StdDevCorr();
        double mean = Mean();
        data.clear();
        for (auto c : old_data){
            double norm_dev = (c - mean) / sigma;
            data.push_back(norm_dev);
        }
    }

    // creazione delle classi per costruire l'istogramma
    vector<class_struct> Classes(){
        vector<class_struct> classes;
        double delta = Delta();
        double mean = Mean();
        double sigma = StdDevCorr();
        double min = Min() - delta;
        double max = Min();
        double top = Max() + (delta / 2);
        while(min < top){
            class_struct cl = class_struct(min, max, data, mean, sigma, delta, n_classes);
            classes.push_back(cl);
            min = max;
            max = min + delta;
        }
        classes[classes.size()-2].freq += classes[classes.size()-1].freq;
        classes[classes.size()-1].freq = 0;
        return classes;
    }

    // stampa delle informazioni sul campione sulla console
    void PrintData(){
        cout<< setprecision(4);
        cout<< filename <<endl <<endl;
        cout<< "Data set size: "    << data.size()  << "\t\tNumber of classes: "    << n_classes    
            << "\t\t\tSize of each class: "<< Delta()        <<endl;
        cout<< "Minimum value: "    << Min()        << "\t\tMaximum value: "        << Max()        <<endl;
        cout<< "Mean value: "       << Mean()       << "\t\tMedian value: "         <<Median()      <<endl;
        cout<< "Std. deviation: "   << StdDev()     << "\t\tCorrected std. dev.: " << StdDevCorr() 
            << "\t\tMean std. dev.: "      << StdDevMean()   <<endl <<endl;
        cout<< string(100, '-') <<endl <<endl;
        cout<< setprecision(0);
    }
    // stampa dell'istogramma sulla console
    void PrintGraph(){
        ios::fmtflags restore = cout.flags();
        cout.setf(ios::fixed, ios::floatfield);
        cout<< setprecision(3) << setw(5);

        cout<< "Hystogram" <<endl <<endl;
        vector<class_struct> classes = Classes();
        int counter = 0;
        for (auto c : classes){
            cout<< c.min << " - " << c.max << '\t' << string(c.freq, '#') <<endl;
            counter += c.freq;
        }
        cout<<endl << "Counted: " << counter <<endl <<endl;

        cout.flags(restore);
    }
    // scrittura su file dei valori delle classi
    void WriteFile(string stage, bool use_norm_dev = 0){
        string filename_out = filename.substr(0, (filename.length() - 4)) + "_hystogram_" + stage + ".txt";
        ofstream output_file (filename_out);
        if (!output_file.is_open())
            throw;
        int norm = 1;
        if (use_norm_dev) norm = data.size();
        vector<class_struct> classes = Classes();
        for (auto c : classes)
            output_file << c.min << '\t' << 1.0 * c.freq / norm << '\t' << c.centroid << '\t' << c.gauss <<endl;
        output_file.close();
    }

    // eliminazione dei dati secondo la regola del 3-sigma
    vector<double> Refine(){
        vector<double> new_data, removed_data;
        double three_sigma = 3.0 * StdDevCorr();
        double mean = Mean();
        for (auto c : data){
            if(c < mean - three_sigma || c > mean + three_sigma)
                removed_data.push_back(c);
            else new_data.push_back(c);
        }
        data = new_data;
        return removed_data;
    }
};

// prototipo delle funzioni
sample_struct FileAnalysis(string filepath, bool use_norm_dev = false);
void PrintRemovedData(vector<double> vect);
void PrintEndofFile();

int main(){
    PrintEndofFile();
    sample_struct sample_1a = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1a.txt");
    sample_struct sample_1b = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1b.txt");
    sample_struct sample_1ab_normdev = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1ab.txt", true);
    sample_struct sample_1c = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1c.txt");
    sample_struct sample_4a = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione4a.txt");
    sample_struct sample_4b = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione4b.txt");
    sample_struct sample_4c = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione4c.txt");
    sample_struct sample_4ab_normdev = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione4ab.txt", true);
    return 0;
}

// step applicati per ogni set di dati
sample_struct FileAnalysis(string filepath, bool use_norm_dev){
    sample_struct sample = sample_struct(filepath);
    sample.ReadFile();
    if (use_norm_dev)
        sample.NormDev();
    sample.PrintData();
    sample.PrintGraph();
    sample.WriteFile("raw", use_norm_dev);
    vector<double> removed_data = sample.Refine();
    while (removed_data.size() && !use_norm_dev){
        PrintEndofFile();
        PrintRemovedData(removed_data);
        sample.PrintData();
        sample.PrintGraph();
        sample.WriteFile("refined");
        removed_data = sample.Refine();
    }
    cout<< "All data is compatible" <<endl;
    PrintEndofFile();
    return sample;
}

//stampa dei dati scartati con i 3-sigma
void PrintRemovedData(vector<double> vect){
    cout<< "Removed data:"<<endl;
    cout<< setprecision(4);
    for (auto c : vect)
        cout<< c <<'\t';
    cout<<endl;
    cout<< setprecision(0);
}

void PrintEndofFile(){
    cout<<endl << string(100, '=') <<endl <<endl;
}