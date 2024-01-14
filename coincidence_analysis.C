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

    Int_t para_detektorow = 0;
    Int_t liczba_pomiarow = 1;

	cout << "Rozpoczeto program" << endl;

    openAndSetupFiles("DataR_run-001.root", "new_Analysis_DataR_run-001.root",20000000, inputFile, inputTree, f_output, energia, czas, channel);

	cout << "Otwarto pliki" << endl;

    timeVectorComputing("czas_do_kalibracji.txt"); // wektor bedzie dostepny pod PA.wektor_czasu

	cout << "N_entries" << n_entries << endl;

    fillTotalVectors(false);

	cout << "Wypelniono podstawowe histogramy" << endl;

    preFitting();

	cout << "Przeprowadzono wstepne dopasowania" << endl;

    findCoincidence();

	cout << "Znaleziono wszystkie koincydencjie" << endl;

    f_output->Write();
}