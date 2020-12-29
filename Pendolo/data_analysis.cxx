#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <stdexcept>
#include <iomanip>
using namespace std;

struct class_struct{
    double min, max, centroid;
    int freq = 0;
    double gauss;
    class_struct(double minimum, double maximum, vector<double> data){
        min = minimum;
        max = maximum;
        centroid = (max + min) / 2.0;
        for(auto c : data)
            if((min <= c) && (c < max))
                freq++;
    }
};

struct sample_struct {
    string filename;
    vector<double> data;
    int n_classes = 15;
    vector<class_struct> classes;

    sample_struct(string filepath){
        filename = filepath;
        ReadFile();
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
    vector<class_struct> Classes(){
        vector<class_struct> classes;
        double delta = Delta();
        for(int i = -1; i <= n_classes; i++){
            double min = Min() + delta * i;
            double max = min + delta;
            class_struct cl = class_struct(min, max, data);
            classes.push_back(cl);
        }
        classes[classes.size()-2].freq += classes[classes.size()-1].freq;
        classes[classes.size()-1].freq = 0;
        return classes;
    }

    void PrintData(){
        cout<< setprecision(4);
        cout<<endl << filename <<endl <<endl;
        cout<< "Data set size: "        << data.size()  << "\t\tNumber of classes: "   << n_classes    << "\t\t\tSize of each class: "<< Delta() <<endl;
        cout<< "Minimum value: "        << Min()        << "\t\tMaximum value: "       << Max()        << "\t\t\tMean value: "        << Mean() <<endl;
        cout<< "Std. deviation: "       << StdDev()     << "\t\tCorrected std. dev.: " << StdDevCorr() << "\t\tMean std. dev.: "      << StdDevMean() <<endl <<endl;
        cout<< string(100, '=') <<endl <<endl;
        cout<< setprecision(0);
    }
    void PrintGraph(){
        ios::fmtflags restore = cout.flags();
        cout.setf(ios::fixed, ios::floatfield);
        cout<<setprecision(3)<<setw(5);

        cout<< "Hystogram" <<endl <<endl;
        vector<class_struct> classes = Classes();
        for (auto c : classes){
            cout<< c.min << " - " << c.max << '\t' << string(c.freq, '#') <<endl;
        }
        cout<< string(100, '=') <<endl <<endl;

        cout.flags(restore);
    }
};

void FilePrintData();
void FilePrintGraph();

int main(){
    // string file_path;
    // while(cin >> file_path){
    //     sample_struct sample_1A = init_sample(file_path);
    // }
    sample_struct sample_1A = sample_struct("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1a.txt");
    sample_1A.PrintData();
    sample_1A.PrintGraph();
    return 0;
}
