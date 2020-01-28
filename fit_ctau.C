#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <TString.h>
#include <TH1D.h>
#include <TClonesArray.h>

using namespace std;

TF1 fExpAndHill(std::string, Double_t, Double_t, char, char);

template<typename TypeNumeric>
TypeNumeric myHeaviside(TypeNumeric x, TypeNumeric xShift=0, TypeNumeric yCentral=1) {
    return x == xShift ? yCentral : (x > xShift ? 1 : 0);
}

/*
char* strConcat(const char *s1, const char *s2) {
    // Adapted from David Heffernan's work at https://stackoverflow.com/a/8465083
    char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}
*/

TF1 fExpAndHill(std::string tfName, Double_t xmin, Double_t xmax, char signHorizontal=1, char signVerticle=1) {
    // myHeaviside<Double_t>((x - [1])*signHorizontal)*[2]*exp(([1] - x)*signHorizontal/[4]) 
    // + myHeaviside<Double_t>(([1] - x)*signHorizontal)*[2]*(1 + (x - [1])*signHorizontal/[3])
    std::string* formula = new string();
    std::string* strXMinusArg1 = "(x - [1])";
    std::string* strXMinusArg1TimesSignHorizontal = signHorizontal > 0 ? strcat("-", strXMinusArg1) : strXMinusArg1;
    std::string* strXMinusArg1TimesSignHorizontalNeg = signHorizontal < 0 ? strcat("-", strXMinusArg1) : strXMinusArg1;
    if (signVerticle < 0) strcat(formula, (string*)"-");
    strcat(formula, (string*)"[2]*(");
    strcat(formula, (string*)"myHeaviside<Double_t>(");
    strcat(formula, strXMinusArg1TimesSignHorizontal);
    strcat(formula, (string*)")");
    strcat(formula, "*exp(");
    strcat(formula, strXMinusArg1TimesSignHorizontalNeg);
    strcat(formula, (string*)"/[4])");
    strcat(formula, (string*)"+");
    strcat(formula, (string*)"myHeaviside<Double_t>(");
    strcat(formula, strXMinusArg1TimesSignHorizontalNeg);
    strcat(formula, (string*)")");
    strcat(formula, (string*)"*(");
    strcat(formula, strXMinusArg1TimesSignHorizontal);
    strcat(formula, (string*)"/[3])");
    strcat(formula, (string*)")");
    std::cout << *formula << std::endl;
    return new TF1(tfName, *formula, xmin, xmax);

}

