#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <stdexcept>
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
        for(int i = -1; i <= data.size(); i++){
            double min = Min() + delta * i;
            double max = min + delta;
            class_struct(min, max, data);
        }
        classes[classes.size()-2].freq += classes[classes.size()-1].freq;
        classes[classes.size()-1].freq = 0;
        return classes;
    }
};

sample_struct InitSample(string filename);
vector<class_struct> CreateClasses(double min, double delta, int n_classes);
void CompileClasses(vector<class_struct> &classes, vector<double> data);
void ConsolePrintData();
void ConsolePrintGraph();
void FilePrintData();
void FilePrintGraph();

int main(){
    // string file_path;
    // while(cin >> file_path){
    //     sample_struct sample_1A = init_sample(file_path);
    // }
    sample_struct sample_1A = sample_struct("C:\\Users\\Admin\\projects\\University\\Pendolo\\Data\\campione1a.txt");
    cout<<sample_1A.Delta();
    return 0;
}
