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


Int_t liczba_pomiarow = 2;

Int_t liczba_par_det = 3;
Int_t liczba_det = 6;
Int_t start_det = 0;
Int_t stop_det = 6;

vector < ULong64_t > czas_koincydencji = {25000, 20000, 30000};

Double_t ile_sigma = 1;

Int_t numberofbins = 300; //liczba binow w widme energergetycznym; zmniejszajac jej wartosc otrzymujemy lepsza statystykew binach ale gorsza rozdzielczosc (mniej pewna wartosc dopasowanej centroidy)

Float_t minimum = 0, maksimum = 1000; //zakres widma energetycznego wyrazone w kanalach

Double_t zliczenia_amplituda = maksimum / numberofbins;

vector < ULong64_t > wektor_czasu;

ULong64_t closest_time_high, closest_time_low, delta_time_high, delta_time_low, delta_time;

ULong64_t n_entries;

Int_t stan, pozycja;

ULong64_t chunk = 0;
ULong64_t begin_chunk, end_chunk, chunks, length_of_chunk;
ULong64_t limit_list = 1000000;
ULong64_t bufor = 100000;
ULong64_t n_entried_entries=0; 

vector < vector < Double_t >> zakres_energii = {{175, 135, 2400, 1700, 1400, 3700}, //Dolny zakres det0 det1 det2...
                                                {260, 200, 3050, 2300, 2000, 4400}}; // Górny zakres det0 det1 det2...
vector < vector < Double_t >> zliczenia(liczba_det, vector <Double_t> (liczba_pomiarow));
vector < vector < Double_t >> blad_zliczenia(liczba_det, vector <Double_t> (liczba_pomiarow));
vector < vector < Double_t >> blad_zliczenia_fixed(liczba_det, vector <Double_t>(liczba_pomiarow));
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


auto h_time = new TH1F("spek_time", "Widmo czasowe", 1e5, 0, 6e5);
auto h_step_time = new TH1F("spek_step_time", "Widmo krokow", 1e5, 0, 6e5);
auto h_rot_time = new TH1F("spek_rot_time", "Widmo momentow obrotu", 1e5, 0, 6e5);



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


void openAndSetupFiles(const char* sourceFileName, const char* outputFileName,
                ULong64_t custom_n_entries,
                TFile*& inputFile, TTree*& inputTree, TFile*& outputFile,
                UShort_t& energiaVar, ULong64_t& czasVar, UShort_t& channelVar){
    inputFile = new TFile(sourceFileName);
    inputTree = dynamic_cast<TTree*>(inputFile->Get("Data_R"));
    if (!inputTree) {
        std::cerr << "Błąd: Nie udało się uzyskać drzewa 'Data_R' z pliku źródłowego." << std::endl;
        return;
    }
    cout << "Przypisanie zmiennych do gałęzi drzewa"<< endl;
    inputTree->SetBranchAddress("Energy", &energiaVar);
    inputTree->SetBranchAddress("Timestamp", &czasVar);
    inputTree->SetBranchAddress("Channel", &channelVar);
    cout << "Utworzenie pliku wyjściowego"<< endl;
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

void timeVectorComputing() {
    inputTree -> GetEntry(n_entries-1);
    std::vector < ULong64_t > vecOfStr; // Sprawdza czy plik jest ok
    vecOfStr.push_back(czas);
    wektor_czasu = vecOfStr;
};


void fillTotalVectors(bool limit_511_kev){
    for (Int_t i = 0; i < n_entries; i++) {
        inputTree -> GetEntry(i);
        if (channel == 7) h_step_time -> Fill(czas / 1e12); // widmo krokow silnika
        if (channel > 5) continue;
        if (channel < start_det || channel >= stop_det) continue;
        channel %= liczba_det;
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


void coincidenceTimeVectors(ULong64_t beginning, ULong64_t ending) {
    if (beginning!=0) {
        czasy.clear();
        pozycje.clear();
        wektor_timestamp.clear();
        wektor_entry.clear();
    };
    czasy.resize(ending - beginning);
    pozycje.resize(ending - beginning);
    wektor_timestamp.resize(liczba_par_det, vector < vector < ULong64_t >> (liczba_pomiarow, czasy));
    wektor_entry.resize(liczba_par_det, vector < vector < Int_t >> (liczba_pomiarow, pozycje));
        for (ULong64_t i = beginning; i < ending; i++) {
            inputTree -> GetEntry(i);
            stan = measurementPoint(wektor_czasu, czas); //stan tarczy: wartosc nieparzysta - obrot, parzysta - pomiar
            pozycja = stan / 2;
            if (i%1000000 == 0) cout<<i<<endl;
            if ((channel > 5) && ((stan % 2) == 0)) continue;
            if ((channel % 2) == 0) continue;
            if ((channel < start_det) || (channel >= stop_det)) continue;
            channel %= liczba_det;
            wektor_timestamp[(channel - 1) / 2][pozycja].push_back(czas); // wypelniany jest wektor z czasem do koincydencji
            wektor_entry[(channel - 1) / 2][pozycja].push_back(i); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
        }
}


bool w_zakresie_elipsy(UShort_t energia_dol, UShort_t energia_gora, UShort_t channel,TF1 * pre_dopasowanie[], Double_t ile_sigma) {
    UShort_t srodek_1 = pre_dopasowanie[channel-1] -> GetParameter(1), srodek_2 = pre_dopasowanie[channel] -> GetParameter(1);
    UShort_t sigma_1 = pre_dopasowanie[channel-1] -> GetParameter(2), sigma_2 = pre_dopasowanie[channel] -> GetParameter(2);
    if (TMath::Power((energia_dol - srodek_1) / (ile_sigma*sigma_1), 2) + TMath::Power((energia_gora - srodek_2) 
        / (ile_sigma*sigma_2), 2) > 1) return true;
    return false;
}


void getChunksTimeVector(ULong64_t *length_of_chunk, ULong64_t *chunks, ULong64_t n_entries, ULong64_t limit_list){
    if (n_entries < limit_list){
        *chunks = 1;
        *length_of_chunk = n_entries;
        return;
    }
    *chunks = (n_entries/limit_list)+1;
    *length_of_chunk = n_entries/(*chunks);
}


void getLimitsNewTimeVector(ULong64_t *chunk, ULong64_t *begin_chunk, ULong64_t *end_chunk, ULong64_t n_entries, ULong64_t length_of_chunk, ULong64_t chunks, ULong64_t bufor){
    cout << "chunks" << chunks << "chunk" << *chunk << endl;
    if (chunks == 1){
        *begin_chunk = 0;
        *end_chunk = n_entries;
        return;
    }
    if (*chunk == 0){
        *begin_chunk = static_cast<ULong64_t>(0);
        *end_chunk = min((*chunk + 1) * length_of_chunk + bufor, n_entries);
        *chunk = *chunk + 1;
        return;
    }
    if (*chunk < (chunks)){
        *begin_chunk = max(*chunk * length_of_chunk - bufor, static_cast<ULong64_t>(0));
        *end_chunk = min((*chunk + 1) * length_of_chunk + bufor, n_entries);
        *chunk = *chunk + 1;
        return;
    }
}


void findCoincidence(){
    getChunksTimeVector(&length_of_chunk, &chunks, n_entries, limit_list);
    getLimitsNewTimeVector(&chunk, &begin_chunk, &end_chunk, n_entries, length_of_chunk, chunks, bufor);
    coincidenceTimeVectors(begin_chunk, end_chunk);
    cout << "Zakres danych to " << begin_chunk << " " << end_chunk << endl;
    for (ULong64_t i = 0; i < n_entries; i++) {
        if ((end_chunk - i) < bufor && i < (n_entries - bufor)) {
            getLimitsNewTimeVector(&chunk, &begin_chunk, &end_chunk, n_entries, length_of_chunk, chunks, bufor);
            coincidenceTimeVectors(begin_chunk, end_chunk);
            }; 
        if (chunk == chunks) {
            break;
            }; 
		inputTree -> GetEntry(i);
		stan = measurementPoint(wektor_czasu, czas); //stan tarczy: wartosc nieparzysta - obrot, parzysta - pomiar
		pozycja = stan / 2;
		// if (w_zakresie_511kev(energia, channel)) continue;
		if (channel > 5) continue;
		if (channel % 2 == 1) continue; // jesli event jest zebrany na kanale innym niz do ktorych byly podlaczone detektory, kod przechodzi do kolejnego eventu
		if (channel < start_det || channel >= stop_det) continue;
		channel %= liczba_det;
		auto closest_higher_time_index = searchClosest(wektor_timestamp[(channel) / 2][pozycja], czas);
		auto closest_lower_time_index = (closest_higher_time_index - 1) * ((closest_higher_time_index - 1) > 0);
		if (closest_higher_time_index > wektor_timestamp[(channel) / 2][pozycja].size() - 1 ||
			closest_lower_time_index > wektor_timestamp[(channel) / 2][pozycja].size() - 1 ||
			closest_higher_time_index > wektor_entry[(channel) / 2][pozycja].size() - 1 ||
			closest_lower_time_index > wektor_entry[(channel) / 2][pozycja].size() - 1) continue;
		closest_time_high = wektor_timestamp[(channel) / 2][pozycja][closest_higher_time_index];
		closest_time_low = wektor_timestamp[(channel) / 2][pozycja][closest_lower_time_index];
		delta_time_high = max(closest_time_high, czas) - min(closest_time_high, czas);
		delta_time_low = max(closest_time_low, czas) - min(closest_time_low, czas);
		auto closest_time_index = closest_higher_time_index;
		delta_time = delta_time_high;
		if (delta_time_high > delta_time_low) {
			closest_time_index = closest_lower_time_index;
			delta_time = delta_time_low;
		}	
		h_delta_time[(channel) / 2] -> Fill(delta_time);
		if (i%1000000 == 0) cout<<i<<endl;
		if (delta_time < czas_koincydencji[(channel) / 2]) {
			energia_other = energia;
			inputTree -> GetEntry(wektor_entry[(channel) / 2][pozycja][closest_time_index]);
			if (channel < start_det || channel >= stop_det) continue;
			channel %= liczba_det;	
			if (w_zakresie_elipsy(energia_other, energia, channel, pre_dopasowanie, ile_sigma)) continue;
			h[channel-1][pozycja] -> Fill(energia_other); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
			h[channel][pozycja] -> Fill(energia); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
			h_2d[(channel) / 2][pozycja] -> Fill(energia_other, energia);
			total_h_2d[(channel) / 2] -> Fill(energia_other, energia);
			n_entried_entries++; //zwiekszana jest liczba oznaczajaca eventy ktore przeszly analize
		} else continue;
	}
	cout<<"Eksport danych"<<endl;
    cout << "Weszlo" << n_entried_entries << endl;
    wektor_czasu.insert(wektor_czasu.end(), czas); // na koniec wektora czasu dorzucany jest koniec pomiaru.
    for (Int_t i = 0; i < liczba_det; i++){
        for (Int_t k = 0; i < liczba_pomiarow; i++){
            h[i][k] -> Write(); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
            h_2d[i / 2][k] -> Write();
            }
        total_h_2d[i / 2] -> Write();
        }
}

void countsExtractionToVector(){
    for (Int_t det = 0; det < liczba_det; det++) {
		for (Int_t i = 0; i < liczba_pomiarow; i++) {
			zliczenia[det][i] = (h[det][i]->IntegralAndError(h[det][i]->FindBin(zakres_energii[0][det]), h[det][i]->FindBin(zakres_energii[1][det]), blad_zliczenia[det][i]));
		}
	}
}


void countsProcessing(Double_t activity, vector < ULong64_t > &timeVector){
	for (Int_t det = 0; det < liczba_det; det++) {
		for (Int_t i = 0; i < liczba_pomiarow; i++) {
			zliczenia[det][i] = 1e12*zliczenia[det][i]/((timeVector[2*i+1] - timeVector[2*i])*activity);
			blad_zliczenia[det][i] = 1e12*blad_zliczenia[det][i]/((timeVector[2*i+1] - timeVector[2*i])*activity);
        }
	}
}


void countsProcessing(vector < ULong64_t > &timeVector){
	for (Int_t det = 0; det < liczba_det; det++) {
		for (Int_t i = 0; i < liczba_pomiarow; i++) {
			zliczenia[det][i] = 1e12*zliczenia[det][i]/((timeVector[2*i+1] - timeVector[2*i]));
			blad_zliczenia[det][i] = 1e12*blad_zliczenia[det][i]/((timeVector[2*i+1] - timeVector[2*i]));
        }
	}
}


void countsProcessing(){
	for (Int_t det = 0; det < liczba_det; det++) {
		for (Int_t i = 0; i < liczba_pomiarow; i++) {
			zliczenia[det][i] = zliczenia[det][i];
			blad_zliczenia[det][i] = blad_zliczenia[det][i];
        }
	}
}


void resultsPrint(){
	for (Int_t det = 0; det < liczba_det; det++) {
		cout << "Liczba zliczen w pomiarze - detektor "<< det << endl;
		for (Int_t i = 0; i < liczba_pomiarow; i++) {
			cout << "\t\t" << zliczenia[det][i];
			cout << "\t\t" << blad_zliczenia[det][i] << endl;
		}
	}
}


void resultsExtraction(const char* extractionFileName){
    char plik[200];
	for (Int_t det = 0; det < liczba_det; det++) {
		sprintf(plik, "%s_%d",extractionFileName, det);
		ofstream aktywnosci(plik);
		for (Int_t i = 0; i < liczba_pomiarow; i++) {
			aktywnosci << zliczenia[det][i] << "\t\t" << blad_zliczenia[det][i] << endl;
		}
		aktywnosci.close();
	}
}