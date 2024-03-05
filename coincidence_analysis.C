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

#include "coincidence_analysis.hpp"

using namespace std;


void coincidence_analysis() {

    // cout << "Liczba analizowanych pomiarow " << liczba_pomiarow << endl;

    // cout << "Przygotowanie wektorow" << endl;

    char name[20];
	char title[100];
    h = new TH1F**[liczba_det];
    h_2d = new TH2F**[liczba_par_det];
    total_h = new TH1F*[liczba_det];
    h_delta_time = new TH1F*[liczba_par_det];
    total_h_2d = new TH2F*[liczba_par_det];
	for (Int_t i = 0; i < liczba_det; i++) {
        if (i % 2 == 0) h_2d[i / 2] = new TH2F*[liczba_pomiarow];
        h[i] = new TH1F*[liczba_pomiarow];
		for (Int_t m = 0; m < liczba_pomiarow; m++) {
			sprintf(name, "spek_%d_det_%d", m + 1, i);
			sprintf(title, "Spektrum pozycji h%d det %d", m + 1, i);
			h[i][m] = new TH1F(name, title, numberofbins, minimum, maksimum);
			sprintf(name, "spek_2D_%d_det_%d", m + 1, i);
			sprintf(title, "Spektrum 2D pozycji h%d det %d", m + 1, i);
			if (i % 2 == 0) h_2d[i / 2][m] = new TH2F(name, title, 100, 0, 1000, 150,  0, 1000);
		}
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

	// cout << "Przygotowanie plikow" << endl;

    openAndSetupFiles("DataR_run-001.root", "new_Analysis_DataR_run-001.root", inputFile, inputTree, f_output, energia, czas, channel);

	// cout << "Otwarto pliki \n Tworzenie wektora obrotow" << endl;

    timeVectorComputing(); 
    // timeVectorComputing("czas_do_analizy.txt"); 

	// cout << "Wektor obrotow policzony \n Przygotowanie histogramow totalnych dla n_entries " << n_entries << endl;

    fillTotalVectors(false);

	// cout << "Wypelniono histogramy h_total \n Przygotowanie do przedswstepnych dopasowan " << n_entries << endl;

    preFitting();

	// cout << "Przeprowadzono wstepne dopasowanie \n Rozpoczecie analizy koincydencyjnej w chunks " << n_entries << endl;

    findCoincidence();

	// cout << "Analiza koincydencyjna zakonczona. Wprowadzono elementow " << n_entried_entries  << "\n Dalsza analiza danych" << n_entries << endl;

	countsExtractionToVector();

	// cout << "Wyeksportowano liczbę zliczeń do wektorów"  << "\n Dalsza analiza danych" << endl;


	//Poniższa funkcja może przyjmować 2, 1 lub żadnego parametru. 
	// - countsProcessing(aktywnosc, wektor_czasu) wtedy obliczona zostaje wydajność pomiaru
	// - countsProcessing(wektor_czasu) wtedy obliczona zostaje liczba zliczeń na jednostkę czasu
	// - countsProcessing() wtedy zostaje liczba zliczeń
	// countsProcessing(38000, wektor_czasu);
	// countsProcessing(wektor_czasu);
	countsProcessing();

	// cout << "Opracowano wymagany wynik."  << "\n Zwrot danych" << endl;

	resultsPrint(); // Wyprintowanie wyników, można pomijać

	resultsExtraction("results.txt"); // Funkcja ta przyjmuje nazwe koncowego pliku do ktorego wyeksportuje wszystkie dane.

    f_output->Write();
}