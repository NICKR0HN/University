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
        double bas = (centroid - mean) / sigma;
        double exponent = (-0.5) * pow(bas, 2.0);
        return exponent;
    }
};

// dichiarazione della struttura per ciascun campione
struct sample_struct {
    // variabili stabili della struttura
    string filepath;
    vector<double> data;
    vector<class_struct> classes;

    // costruttore
    sample_struct(string filename){
        filepath = "./data/" + filename;
    }

    // acquisizione dei dati
    void ReadFile(){
        ifstream input_file(filepath);
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
        return ((sort_data[half] + sort_data[half + 1]) / 2.0);
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
};


