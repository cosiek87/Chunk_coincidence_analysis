////////////////////////////////////////////////////////////////////////////
//
//	Program do analizy plikow zawierajacych cale sekwencje pomiarow przeprowadzonych przy uzyciu ukladu
//  Caen DT5730 oraz programu Compass. Skrypt ten w glownej mierze sluzy tylko do pomiarow z 09.2021
//  Jednakze rdzen sortowania eventow i dopasowywania funckji jest ten dla wszystkich pomiarow
//  w ktorych uzywany byl uklad firmy Caen
//
//	Autor: Przemyslaw Sekowski
//
//	email: przemyslaw.sekowski@fuw.edu.pl
//
///////////////////////////////////////////////////////////////////////////

#include <iostream>

using namespace std;

//Glowna czesc skryptu
void coincidence_analysis_time_export() {
    char plik[200];
    ULong64_t czas;
    UShort_t channel;

    Char_t nazwa[200] = "DataR_run-001.root"; //plik zrodlowy z danymi

    auto * file_input = new TFile(nazwa); //tworzone sa tu plik wyjsciowy i wszystkie histogramy, a takze do zmiennych energia, czas i channel przypisywane sa wartosci eventow
    auto * tree_input = (TTree * ) file_input -> Get("Data_R");

    tree_input -> SetBranchAddress("Timestamp", & czas);
    tree_input -> SetBranchAddress("Channel", & channel);

    ULong64_t nentries = (ULong64_t) tree_input -> GetEntries();

    sprintf(plik, "kroki_%s.txt", nazwa);
    ofstream kroki(plik);

    for (ULong64_t i = 0; i < nentries; i++) {
        tree_input -> GetEntry(i);
        if (channel == 7) kroki << czas << endl;
    }
    kroki.close();
	}