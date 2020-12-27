#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
using namespace std;

struct campione
{
	vector<double> dati;				//vettore con i dati
	vector<int> classi;					//vettore con le frequenze delle misure (n. misure per classe di frequenza)
	double max, min;					//massimo e minimo del campione, da calcolare
	int C;								//numero delle classi di frequenza
	double delta;						//larghezza delle classi di frequenza
	double media;						//valore medio dei dati
	string nome;						//nome del file da cui proviene il campione
	double dev_standard_singola_misura; //deviazione standard della singola misura
	double dev_standard_campionaria;	//deviazione standard campionaria
	double dev_standard_media;			//deviazione standard sulla media (misure statisticamente indipendenti)
	//per centroide
	vector<double> punti_ascissa;
};

int N = 50;								//dimensione del campione
const int numero_classi_frequenza = 15; //numero delle classi di frequenza

void init_campione(campione &c, string filename);
void print_campione(campione);
void set_classi(campione &c);
void classi_out(campione);
void print_classi(campione);

double media_campione(campione);
double dev_standard_camp(campione);
double dev_standard_sing_mis(campione);
double dev_standard_med(campione);

void centroide_bin(campione);
void centroide_out(campione);

int main()
{

	//Di seguito sono dichiarate le variabili di tipo campione con cui si identificano i vari set di dati

	campione annachiara1; //campione 1A
	campione tommaso1;	  //campione 1B
	campione annachiara4; //campione 4A
	campione tommaso4;	  //campione 4B
	campione ga1;		  //campione 1C
	campione ga4;		  //campione 4C
	campione ab1;		  //campione 1AB
	campione ab4;		  //campione 4AB

	string filedati1a = "campione1a.txt";	//campione 1A
	string filedati1t = "campione1t.txt";	//campione 1B
	string filedati4a = "campione4a.txt";	//campione 4A
	string filedati4t = "campione4t.txt";	//campione 4B
	string filedati1c = "campione1c.txt";	//campione 1C
	string filedati4c = "campione4c.txt";	//campione 4C
	string filedati1ab = "campione1ab.txt"; //campione 1AB
	string filedati4ab = "campione4ab.txt"; //campione 4AB

	init_campione(annachiara1, filedati1a);
	//print_campione(annachiara1);
	classi_out(annachiara1);
	centroide_bin(annachiara1);
	centroide_out(annachiara1);

	/*

	init_campione(tommaso1, filedati1t);
	print_campione(tommaso1);
	classi_out(tommaso1);


	init_campione(annachiara4, filedati4a);
	print_campione(annachiara4);
	classi_out(annachiara4);


	init_campione(tommaso4, filedati4t);
	print_campione(tommaso4);
	classi_out(tommaso4);



	N=100;


	init_campione(ga1, filedati1c);
	print_campione(ga1);
	classi_out(ga1);

	init_campione(ga4, filedati4c);
	print_campione(ga4);
	classi_out(ga4);
	*/
	/*
	init_campione(ab1, filedati1ab);
	print_campione(ab1);
	classi_out(ab1);


	init_campione(ab4, filedati4ab);
	print_campione(ab4);
	classi_out(ab4);
	*/

	return 0;
}

//Inizializzazione del campione di dati
void init_campione(campione &c, string filename)
{
	c.nome = filename; //Imposto c.nome = al nome del file da cui prendo i dati
	ifstream ifile(filename);
	if (ifile.is_open())
	{
		for (int i = 0; i < N; i++)
		{ //inizializza il vettore inserendo i dati leggendoli da file
			double temp;
			ifile >> temp;
			c.dati.push_back(temp);
		}
		ifile.close();
	}

	//Calcolo dei valori massimo e minimo misurati
	c.max = 0.;
	c.min = 10.;
	for (int j = 0; j < N; j++)
	{ //calcola massimo e minimo del campione
		if (c.dati.at(j) > c.max)
			c.max = c.dati.at(j);
		if (c.dati.at(j) < c.min)
			c.min = c.dati.at(j);
	}

	c.C = numero_classi_frequenza; //numero delle classi di frequenza, qui viene fissato uguale per tutti i campioni

	c.delta = (c.max - c.min) / c.C; //calcolo ampiezza delle classi di frequenza

	c.media = media_campione(c);
	c.dev_standard_singola_misura = dev_standard_sing_mis(c);
	c.dev_standard_campionaria = dev_standard_camp(c);
	c.dev_standard_media = dev_standard_med(c);

	set_classi(c);
}

void print_campione(campione c)
{
	cout << "**********************" << endl
		 << endl
		 << "CAMPIONE	" << c.nome << endl
		 << endl;
	for (auto a : c.dati)
	{
		cout << a << "\t";
	}
	cout << endl;
	cout << "Il massimo e il minimo valgono:	max = " << c.max << "\t min = " << c.min << endl;
	cout << "Delta vale: " << c.delta << endl;
	cout << "La media vale:	" << c.media << endl;
	cout << "La deviazione standard campionaria vale:	" << c.dev_standard_campionaria << endl;
	cout << "La deviazione standard sulla singola misura vale:	" << c.dev_standard_singola_misura << endl;
	cout << "La deviazione standard sulla media vale:	" << c.dev_standard_media << endl;

	print_classi(c);
}

void set_classi(campione &c)
{
	for (int j = 0; j < numero_classi_frequenza; j++)
		c.classi.push_back(0.0);

	for (int i = 0; i < numero_classi_frequenza; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if ((c.dati.at(j) >= (c.min + i * c.delta)) && (c.dati.at(j) < (c.min + (i + 1) * c.delta)))
				c.classi.at(i)++;
		}
	}
	c.classi.at(c.C - 1)++; //Aggiungo 1 all'ultima classe di frequenza perchè nella condizione if ho messo < e non <= (altrimenti avrei alterato tutte le altre classi di frequenza)

	/*
cout<<endl<<"Le classi di frequenza con i relativi valori valgono: "<<endl;
	for(int j = 0; j < c.C; j++)
	cout<<"Classe "<< (j+1) <<":\t da  " << c.min+j*c.delta << "\t a  "<<c.min + (j+1)*c.delta <<"\t"<< c.classi.at(j) <<endl;
cout<<endl;
*/
}

void print_classi(campione c)
{
	cout << endl
		 << "Le classi di frequenza con i relativi valori valgono: " << endl;
	for (int j = 0; j < c.C; j++)
		cout << "Classe " << (j + 1) << ":\t da  " << c.min + j * c.delta << "\t a  " << c.min + (j + 1) * c.delta << "\t" << c.classi.at(j) << endl;

	cout << endl;
}

void classi_out(campione c)
{
	string fileout = c.nome.substr(0, c.nome.length() - 4) + "_out.txt";
	//fileout + "_out";
	ofstream ofile(fileout);
	if (ofile.is_open())
	{
		for (int j = 0; j < c.C; j++)
			ofile << c.min + j * c.delta << "\t" << c.classi.at(j) << endl;
	}
	ofile.close();
}

double media_campione(campione c)
{
	double somma = 0., media;
	for (int i = 0; i < N; i++)
	{
		somma += c.dati.at(i);
	}
	media = somma / N;
	return media;
}

double dev_standard_camp(campione c)
{
	double sigma, sommatoria = 0.;
	for (int i = 0; i < N; i++)
	{
		sommatoria += pow((c.dati.at(i) - c.media), 2.0);
	}
	sigma = sqrt(sommatoria / (N));
	return sigma;
}

double dev_standard_sing_mis(campione c)
{
	double S, sommatoria = 0.;
	for (int i = 0; i < N; i++)
	{
		sommatoria += pow((c.dati.at(i) - c.media), 2.0);
	}
	S = sqrt(sommatoria / (N - 1));
	return S;
}

double dev_standard_med(campione c)
{
	double sigma_x;
	sigma_x = c.dev_standard_singola_misura / (sqrt(N));

	return sigma_x;
}

void centroide_bin(campione c)
{
	double ascissa;
	double ordinata;

	for (int i = 0; i < numero_classi_frequenza; i++)
	{
		ascissa = c.min + (i * (c.delta)) + ((c.delta) / 2);
		c.punti_ascissa.push_back(ascissa); //scrivo ascissa centroide
	}
	cout << endl;
	for (auto a : c.punti_ascissa)
	{
		cout << a << endl;
	}
	for (auto a : c.classi)
	{
		cout << a << endl;
	}
}

void centroide_out(campione c)
{
	string fileout = "centroide.txt";
	ofstream ofilec(fileout);
	if (ofilec.is_open())
	{
		for (int j = 0; j < numero_classi_frequenza; j++)
			ofilec << c.punti_ascissa[j] << "\t" << c.classi[j] << endl;
	}
	ofilec.close();
}
