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
};

struct sample_struct {
    string filename;
    vector<double> data;
    int data_len;
    double max, min, mean, delta;
    double std_dev, std_dev_corr, std_dev_mean;
    int n_classes = 15;
    vector<class_struct> classes;
};

sample_struct InitSample(string filename);
vector<double> ReadFile(string filename);
double GetMax(vector<double> data);
double GetMin(vector<double> data);
double Mean(vector<double> data, int data_len);
double GetSummation(vector<double> data, double mean);
vector<class_struct> CreateClasses(double min, double delta, int n_classes);
void CompileClasses(vector<class_struct> &classes, vector<double> data);

int main(){
    // string file_path;
    // while(cin >> file_path){
    //     sample_struct sample_1A = init_sample(file_path);
    // }
    sample_struct sample_1A = InitSample("C:\\Users\\Admin\\projects\\Pendolo\\campione1a.txt");
    return 0;
}


sample_struct InitSample(string filename){
    sample_struct sample;
    sample.filename = filename;
    sample.data = ReadFile(sample.filename);
    sample.data_len = sample.data.size();
    sample.max = GetMax(sample.data);
    sample.min = GetMin(sample.data);
    sample.delta = (sample.max - sample.min) / sample.n_classes;
    sample.mean = Mean(sample.data, sample.data_len);
    double summation = GetSummation(sample.data, sample.mean);
    sample.std_dev = sqrt(summation / sample.data_len);
    sample.std_dev_corr = sqrt(summation / (sample.data_len - 1.0));
    sample.std_dev_mean = sample.std_dev_corr / sqrt(sample.data_len);
    sample.classes = CreateClasses(sample.min, sample.delta, sample.n_classes);
    CompileClasses(sample.classes, sample.data);
    return sample;
}

vector<double> ReadFile(string filename){
    vector<double> data;
    ifstream input_file(filename);
    if (!input_file.is_open()){
        cout<<"Error: cannot read the file"<<endl;
        throw;
    }
    double value;
    while (input_file >> value){
        data.push_back(value);
    }
    input_file.close();
    return data;
}

double GetMax(vector<double> data){
    double max = data[0];
    for (auto c : data){
        if (c > max){
            max = c;
        }
    }
    return max;
}

double GetMin(vector<double> data){
    double min = data[0];
    for (auto c : data){
        if (c < min){
            min = c;
        }
    }
    return min;
}

double Mean(vector<double> data, int data_len){
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
    myClass.centroid = (myClass.min + myClass.max)/2.0;
    classes.push_back(myClass);
    for (int i = 0; i <= n_classes; i++){
        myClass.min = min + delta * i;
        myClass.max = myClass.min + delta;
        myClass.centroid = (myClass.min + myClass.max)/2.0;
        classes.push_back(myClass);
    }
    return classes;
}

void CompileClasses(vector<class_struct> &classes, vector<double> data){
    for(auto c : classes){
        for(auto d : data){
            if((c.min <= d) && (d < c.max)){
                c.freq++;
            }
        }
    }
    classes[classes.size()-2].freq += classes[classes.size()-1].freq;
    classes[classes.size()-1].freq = 0;
}