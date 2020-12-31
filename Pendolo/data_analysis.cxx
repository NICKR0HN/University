#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <algorithm> 
using namespace std;

struct class_struct{
    double min, max, centroid;
    int freq = 0;
    double gauss;
    class_struct(double minimum, double maximum, vector<double> data, double mean, double sigma){
        min = minimum;
        max = maximum;
        centroid = (max + min) / 2.0;
        for(auto c : data)
            if((min <= c) && (c < max))
                freq++;
        double den = sqrt(2.0 * M_PI) * sigma;
        gauss = exp(Exponent(mean, sigma)) / den;
    }
    double Exponent(double mean, double sigma){
        double num = pow((centroid - mean), 2.0);
        double den = 2.0 * pow(sigma, 2.0);
        double exponent = (-1.0) * (num / den);
        return exponent;
    }
};

struct sample_struct {
    string filename;
    vector<double> data;
    int n_classes = 20;
    vector<class_struct> classes;

    sample_struct(string filepath){
        filename = filepath;
    }

    void ReadFile(){
        ifstream input_file(filename);
        if (!input_file.is_open())
            throw;
        double value;
        while (input_file >> value)
            data.push_back(value);
        input_file.close();
    }
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
    double Delta(){
        double delta = (Max() - Min()) / n_classes;
        return delta;
    }
    double Summation(){
        double summation = 0.0;
        double m = Mean();
        for(auto c : data)
            summation += pow((c - m), 2.0);
        return summation;
    }
    double GetStdDev(int data_len){
        double std_dev = sqrt(Summation() / data_len);
        return std_dev;
    }
    double StdDev(){
        double std_dev = GetStdDev(data.size());
        return std_dev;
    }
    double StdDevCorr(){
        double std_dev_corr = GetStdDev(data.size() - 1);
        return std_dev_corr;
    }
    double StdDevMean(){
        double std_dev_mean = StdDevCorr()/sqrt(data.size());
        return std_dev_mean;
    }

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

    vector<class_struct> Classes(){
        vector<class_struct> classes;
        double delta = Delta();
        double mean = Mean();
        double sigma = StdDevCorr();
        double min = Min() - delta;
        double max = Min();
        double top = Max() + (delta / 2);
        while(min < top){
            class_struct cl = class_struct(min, max, data, mean, sigma);
            classes.push_back(cl);
            min = max;
            max = min + delta;
        }
        classes[classes.size()-2].freq += classes[classes.size()-1].freq;
        classes[classes.size()-1].freq = 0;
        return classes;
    }

    void PrintData(){
        cout<< setprecision(4);
        cout << filename <<endl <<endl;
        cout<< "Data set size: "    << data.size()  << "\t\tNumber of classes: "    << n_classes    << "\t\t\tSize of each class: "<< Delta()        <<endl;
        cout<< "Minimum value: "    << Min()        << "\t\tMaximum value: "        << Max()        <<endl;
        cout<< "Mean value: "       << Mean()       << "\t\tMedian value: "         <<Median()      <<endl;
        cout<< "Std. deviation: "   << StdDev()     << "\t\tCorrected std. dev.: " << StdDevCorr() << "\t\tMean std. dev.: "      << StdDevMean()   <<endl <<endl;
        cout<< string(100, '-') <<endl <<endl;
        cout<< setprecision(0);
    }
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

sample_struct FileAnalysis(string filepath, bool use_norm_dev = false);
void PrintRemovedData(vector<double> vect);
void PrintEndofFile();

int main(){
    // string file_path;
    // while(cin >> file_path){
    //     sample_struct sample_1A = init_sample(file_path);
    // }
    PrintEndofFile();
    // sample_struct sample_1a = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1a.txt");
    // sample_struct sample_1b = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1b.txt");
    sample_struct sample_1ab = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1ab.txt", true);
    // sample_struct sample_1c = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1c.txt");
    // sample_struct sample_4a = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione4a.txt");
    // sample_struct sample_4b = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione4b.txt");
    // sample_struct sample_4c = FileAnalysis("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione4c.txt");
    return 0;
}

sample_struct FileAnalysis(string filepath, bool use_norm_dev){
    sample_struct sample = sample_struct(filepath);
    sample.ReadFile();
    if (use_norm_dev)
        sample.NormDev();
    sample.PrintData();
    sample.PrintGraph();
    sample.WriteFile("raw", use_norm_dev);
    vector<double> removed_data = sample.Refine();
    while (removed_data.size()){
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