#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
using namespace std;

struct class_struct{
    double min, max;
    int freq = 0;
};

struct sample_struct {
    string file_name;
    vector<double> data;
    int data_len;
    double max, min, mean, delta;
    double std_dev, std_dev_corr, std_dev_mean;
    int n_classes = 15;
    vector<class_struct> classes;
};

sample_struct init_sample(string filename);
double mean(vector<double> data, int data_len);
double GetSummation(vector<double> data, double mean);
vector<class_struct> CreateClasses(double min, double delta, int n_classes);
//void CompileClasses()

int main(){
    sample_struct sample_1A = init_sample("C:\\Users\\Admin\\projects\\Pendolo\\campione1a.txt");
    return 0;
}


sample_struct init_sample(string filename){
    sample_struct sample;
    sample.file_name = filename;
    ifstream input_file(filename);
    if (!input_file.is_open()){
        throw "Error: cannot read the file";
    }
    int data_len = 0;
    double value, max = 0.0, min = 10.0;
    while (input_file >> value){
        sample.data.push_back(value);
        if (value > max) max = value;
        if (value < min) min = value;
        data_len++;
    }
    input_file.close();

    sample.data_len = data_len;
    sample.max = max;
    sample.min = min;
    sample.delta = (max - min) / sample.n_classes;
    sample.mean = mean(sample.data, data_len);
    double summation = GetSummation(sample.data, sample.mean);
    sample.std_dev = sqrt(summation / data_len);
    sample.std_dev_corr = sqrt(summation / (data_len - 1.0));
    sample.std_dev_mean = sample.std_dev_corr / sqrt(data_len);
    sample.classes = CreateClasses(min, sample.delta, sample.n_classes);
    return sample;
}

double mean(vector<double> data, int data_len){
    double sum = 0.0;
    for (auto c : data){
        sum += c;
    }
    double mean = sum / data_len;
    return mean;
}
 
double GetSummation(vector<double> data, double mean){
    double summation = 0.0;
    for(auto c : data){
        summation += pow((c - mean), 2.0);
    }
    return summation;
}

vector<class_struct> CreateClasses(double min, double delta, int n_classes){
    vector<class_struct> classes;
    class_struct myClass;
    myClass.min = min - delta;
    myClass.max = min;
    classes.push_back(myClass);
    for (int i = 0; i <= n_classes; i++){
        myClass.min = min + delta * i;
        myClass.max = myClass.min + delta;
        classes.push_back(myClass);
    }
    return classes;
}