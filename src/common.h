/* Copyright (C) 2010 Romain Dubessy */
#ifndef COMMON_H
#define COMMON_H
#include <string>
#include <map>
using namespace std;
/*! \brief This type contains the options, as a map of string pairs, keys and
 * values */
typedef map<string,string> ConfigMap;
/*! \brief This method reads a string and add an entry in the map. */
bool readString(string,ConfigMap &,string &, bool);
/*! \brief This method reads the cmd lines options and trigger the program
 * execution. */
bool parseOptions(const int, char *[],ConfigMap &);
/*! \brief This method displays a convevience 'usage' screen. */
void printUsage(char []);
/*! \brief This method reads a string and remove any blanks or comments. */
string strip(string);
/*! \brief This method populates the map. */
bool parseConfig(ConfigMap &);
int getConfig(ConfigMap &,const string &,int);
double getConfig(ConfigMap &,const string &,double);
string getConfig(ConfigMap &,const string &,const string &);
#endif
/* common.h */
