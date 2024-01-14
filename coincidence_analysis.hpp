////////////////////////////////////////////////////////////////////////////
//
//	Program do analizy plikow zawierajacych cale sekwencje pomiarow przeprowadzonych przy uzyciu ukladu
//  Caen DT5730 oraz programu Compass.
//  Jednakze rdzen sortowania eventow i dopasowywania funckji jest ten dla wszystkich pomiarow
//  w ktorych uzywany byl uklad firmy Caen
//  
//
//
//  ver. 2.0 BETA: Wprowadzenie dzielenia wektora koincydencji na czesci
//
//	Autor: Przemyslaw Sekowski
//
//	email: przemyslaw.sekowski@fuw.edu.pl
//
///////////////////////////////////////////////////////////////////////////
#include <TFile.h>

#include <TTree.h>

#include <TH1F.h>

#include <TH2F.h>

#include <TF1.h>

#include <TMath.h>

#include <iostream>

#include <vector>

#include <algorithm>

#include <fstream>

#include <cmath>

#include <cstring>

using namespace std;

typedef unsigned short UShort_t;

UShort_t energia;

UShort_t energia_other;

ULong64_t czas;

UShort_t channel;

TFile* inputFile = nullptr; // plik źródłowy

TTree* inputTree = nullptr; // drzewo

TFile* f_output = nullptr; // plik wyjściowy



Double_t gausswithlinearbkg(Double_t * xarg, Double_t * par) {
    Double_t x = xarg[0], result = 0.;
    result = par[0] / (par[2] * TMath::Sqrt(2 * TMath::Pi())) * TMath::Exp((-TMath::Power((x - par[1]), 2)) / (2 * TMath::Power(par[2], 2))) + par[3] * x + par[4];
    return result;
}


Int_t liczba_pomiarow;

Int_t liczba_par_det, start_det, stop_det, liczba_det;

vector < ULong64_t > czas_koincydencji = {25000, 20000, 30000};

Double_t ile_sigma = 1;

Int_t numberofbins = 300; //liczba binow w widme energergetycznym; zmniejszajac jej wartosc otrzymujemy lepsza statystykew binach ale gorsza rozdzielczosc (mniej pewna wartosc dopasowanej centroidy)

Float_t minimum = 0, maksimum = 1000; //zakres widma energetycznego wyrazone w kanalach

Double_t zliczenia_amplituda = maksimum / numberofbins;

vector < ULong64_t > wektor_czasu;

Int_t n_entries;

Int_t stan, pozycja;

vector < vector < Double_t >> zakres_energii;
vector < vector < Double_t >> zliczenia;
vector < vector < Double_t >> blad_zliczenia;
vector < vector < Double_t >> blad_zliczenia_fixed;
vector < vector < vector < ULong64_t >>> wektor_timestamp;
vector < vector < vector < Int_t >>> wektor_entry;
vector < ULong64_t > czasy;
vector < Int_t > pozycje;

TH1F*** h;
TH2F*** h_2d;
TH1F** total_h;
TH2F** total_h_2d;
TH1F** h_delta_time;

TF1 * pre_dopasowanie[6];


bool within511kevRange(vector < vector < Double_t >> &zakres_energii, UShort_t E, UShort_t Ch){
    bool within_range = true;
    if (E > zakres_energii[Ch][0] && E < zakres_energii[Ch][1]) within_range = false;
    return within_range;
}


ULong64_t searchClosest(std::vector < ULong64_t > & sorted_array, ULong64_t x) {
    auto iter_geq = std::lower_bound(sorted_array.begin(), sorted_array.end(), x);
    if (iter_geq == sorted_array.begin()) return 0;
    ULong64_t a = * (iter_geq - 1);
    ULong64_t b = * (iter_geq);
    if (fabs(x - a) < fabs(x - b)) return iter_geq - sorted_array.begin() - 1;
    return iter_geq - sorted_array.begin();
}


Int_t measurementPoint(vector < ULong64_t > &wektor_czasu, ULong64_t czas) {
    Int_t pozycja = 0;
    for (Int_t i = 0; i < wektor_czasu.size(); i++) {
        pozycja = i;
        if (czas < wektor_czasu[i]) break;
    }
    return pozycja;
}


void prepareVariables(Int_t para_detektorow, Int_t lp){
    char name[100];
    char title[100];
    liczba_pomiarow	= lp;
    if (para_detektorow == -1) {
        liczba_par_det = 3;
        liczba_det = 6;
        start_det = 0;
        stop_det = 6;
    } else {
        liczba_par_det = 1;
        liczba_det = 2;
        start_det = 2 * para_detektorow;
        stop_det = start_det + 2;
    };
    zakres_energii = {{175, 135, 2400, 1700, 1400, 3700}, //Dolny zakres det0 det1 det2...
                    {260, 200, 3050, 2300, 2000, 4400}}; // Górny zakres det0 det1 det2...
    zliczenia.resize(liczba_det, vector <Double_t> (liczba_pomiarow));
    blad_zliczenia.resize(liczba_det, vector <Double_t> (liczba_pomiarow));
    blad_zliczenia_fixed.resize(liczba_det, vector <Double_t>(liczba_pomiarow));
    h = new TH1F**[liczba_det];
    h_2d = new TH2F**[liczba_par_det];
    total_h = new TH1F*[liczba_det];
    h_delta_time = new TH1F*[liczba_par_det];
    for (Int_t i = 0; i < liczba_det; i++) {
        h[i] = new TH1F*[liczba_pomiarow];
        if (i % 2 == 0) h_2d[i] = new TH2F*[liczba_pomiarow];
        for (Int_t m = 0; m < liczba_pomiarow; m++) {
            sprintf(name, "spek_%d_det_%d", m + 1, i);
            sprintf(title, "Spektrum pozycji h%d det %d", m + 1, i);
            h[i][m] = new TH1F(name, title, numberofbins, minimum, maksimum);
            sprintf(name, "spek_2D_%d_det_%d", m + 1, i);
            sprintf(title, "Spektrum 2D pozycji h%d det %d", m + 1, i);
            if (i % 2 == 0) h_2d[i / 2][m] = new TH2F(name, title, 100, 0, 1000, 150, 0, 1000);
        }
        cout << "Przygotowano wszystkie histogramy h " << i << endl;
        sprintf(name, "spek_total_det_%i", i);
        sprintf(title, "Spektrum calkowite det %i", i);
        total_h[i] = new TH1F(name, title, 4000, minimum, maksimum);
        sprintf(name, "spek_delta_det_%i", i);
        sprintf(title, "Spektrum delta time det %i", i);
        h_delta_time[i] = new TH1F(name, title, 5e3, 0, 1e5);
    }
    for (Int_t i = 0; i < liczba_par_det; i++) {
        sprintf(name, "spek_total_2d_det_%d", i);
        sprintf(title, "Spektrum calkowite 2D det %d", i);
        total_h_2d[i] = new TH2F(name, title, 400, minimum, maksimum, 400, minimum, maksimum);
    }
};


void openAndSetupFiles(const char* sourceFileName, const char* outputFileName,
                TFile*& inputFile, TTree*& inputTree,
                TFile*& outputFile,
                UShort_t& energiaVar,
                ULong64_t& czasVar, UShort_t& channelVar){
    inputFile = new TFile(sourceFileName);
    inputTree = dynamic_cast<TTree*>(inputFile->Get("Data_R"));
    if (!inputTree) {
        std::cerr << "Błąd: Nie udało się uzyskać drzewa 'Data_R' z pliku źródłowego." << std::endl;
        return;
    }
    // Przypisanie zmiennych do gałęzi drzewa
    inputTree->SetBranchAddress("Energy", &energiaVar);
    inputTree->SetBranchAddress("Timestamp", &czasVar);
    inputTree->SetBranchAddress("Channel", &channelVar);
    // Utworzenie pliku wyjściowego
    outputFile = new TFile(outputFileName, "RECREATE");
    n_entries = ((Int_t) inputTree -> GetEntries());
};


void openAndSetupFiles(const char* sourceFileName, const char* outputFileName, Int_t custom_n_entries,
                TFile*& inputFile, TTree*& inputTree,
                TFile*& outputFile,
                UShort_t& energiaVar,
                ULong64_t& czasVar, UShort_t& channelVar){
    inputFile = new TFile(sourceFileName);
    inputTree = dynamic_cast<TTree*>(inputFile->Get("Data_R"));
    if (!inputTree) {
        std::cerr << "Błąd: Nie udało się uzyskać drzewa 'Data_R' z pliku źródłowego." << std::endl;
        return;
    }
    // Przypisanie zmiennych do gałęzi drzewa
    inputTree->SetBranchAddress("Energy", &energiaVar);
    inputTree->SetBranchAddress("Timestamp", &czasVar);
    inputTree->SetBranchAddress("Channel", &channelVar);
    // Utworzenie pliku wyjściowego
    outputFile = new TFile(outputFileName, "RECREATE");
    n_entries = custom_n_entries;
};


void timeVectorComputing(const char* fileName) {
    std::ifstream in (fileName);
    std::vector < ULong64_t > vecOfStr; // Sprawdza czy plik jest ok
    std::string str; // Czyta kolejne linijki dopoki plik sie nie skonczy
    while (std::getline( in , str)) { // Jesli string jest niezerowej dlugosci to jest zapisywany do wektora czasu
        
        if (str.size() > 0) vecOfStr.push_back(stod(str));
    } in .close(); //Zamyka plik
    wektor_czasu = vecOfStr;
};


void fillTotalVectors(bool limit_511_kev){
    for (Int_t i = 0; i < n_entries; i++) {
        inputTree -> GetEntry(i);
        // if (i < 100) cout <<"Czas dla "<< i << ": " << czas << endl;
        if (channel > 5) continue;
        if (channel < start_det || channel >= stop_det) continue;
        channel %= liczba_det;
        if (i%1000000 == 0) cout << "obrabiam " << i << endl;
        total_h[channel] -> Fill(energia);
    }
    for (Int_t i = 0; i < liczba_det; i++) {
        total_h[i] -> Write();
    }
} 


void preFitting(){
    for (Int_t det = 0; det < liczba_det; det++) {
        pre_dopasowanie[det] = new TF1("dopasowanie", gausswithlinearbkg, zakres_energii[0][start_det+det], zakres_energii[1][start_det+det], 5);
        pre_dopasowanie[det] -> SetParameters(20000, (zakres_energii[0][start_det+det] + zakres_energii[1][start_det+det]) / 2, 25, -0.00001, 10);
        pre_dopasowanie[det] -> SetParNames("Amplituda", "Srednia", "Sigma", "A", "B");
        pre_dopasowanie[det] -> SetParLimits(0, 20 * zliczenia_amplituda, 200000 * zliczenia_amplituda);
        pre_dopasowanie[det] -> SetParLimits(1, zakres_energii[0][start_det+det], zakres_energii[1][start_det+det]);
        pre_dopasowanie[det] -> SetParLimits(2, 25, 100);
        pre_dopasowanie[det] -> SetParLimits(3, -0.1, 0.1);
        pre_dopasowanie[det] -> SetParLimits(4, -10000, 10000);
        total_h[det] -> Fit("dopasowanie", "LMQR0", "", zakres_energii[0][start_det+det], zakres_energii[1][start_det+det]);
    }
    for (int i = 0 ; i<liczba_det; i++) {
        cout << total_h[i] -> GetEntries() << endl;
        cout << zakres_energii[0][i] << " " << zakres_energii[1][i] << endl;
        cout << "centroida " << pre_dopasowanie[i] -> GetParameter(1) << " sigma "<< pre_dopasowanie[i] -> GetParameter(2) << endl;;
    }
}


void coincidenceTimeVectors(Int_t beginning, Int_t end) {
    if (beginning!=0) {
        czasy.clear();
        pozycje.clear();
        wektor_timestamp.clear();
        wektor_entry.clear();
    };
    czasy.resize(end - beginning);
    pozycje.resize(end - beginning);
    wektor_timestamp.resize(liczba_par_det, vector < vector < ULong64_t >> (liczba_pomiarow, czasy));
    wektor_entry.resize(liczba_par_det, vector < vector < Int_t >> (liczba_pomiarow, pozycje));
    cout << "start petli w coincidenceTimeVectors" << endl;
    for (Int_t i = beginning; i < end; i++) {
        inputTree -> GetEntry(i);
        if (i%100==0) cout << i << endl;
        if (i < 100) cout <<"Czas dla "<< i << ": " << czas << endl;
        stan = measurementPoint(wektor_czasu, czas); //stan tarczy: wartosc nieparzysta - obrot, parzysta - pomiar
        cout << "Stan: " << stan << endl;
        pozycja = stan / 2;
        cout << "Pozycja: " << pozycja << endl;
        if ((channel > 5) && ((stan % 2) == 0)) continue;
        if ((channel % 2) == 0) continue;
        if ((channel < start_det) || (channel >= stop_det)) continue;
        channel %= liczba_det;
        wektor_timestamp[(channel - 1) / 2][pozycja].push_back(czas); // wypelniany jest wektor z czasem do koincydencji
        wektor_entry[(channel - 1) / 2][pozycja].push_back(i); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
    }
    cout << "stworzyla sie macierz koincydencji entry od "<< beginning <<" do "<< end << endl;
}









