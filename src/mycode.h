/* Copyright (C) 2010 Romain Dubessy */
#ifndef MYCODE_H
#define MYCODE_H
#include <iostream>
#include <sstream>
#include <fstream>                               //For ifstream / ofstream
#include <cmath>                                 //For fabs()
#include <time.h>
#include <stdlib.h>
#include "octtree.h"
#include "common.h"
/*! \file */
/* Options: data structure containing the problem's parameters. {{{ */
/*! \brief Data structure containing the program parameters.
 */
struct Options {
    int nEsp;           /*!<\brief Number of charges species. */
    double size;        /*!<\brief Size of the initial distribution. */
    int nTot;           /*!<\brief Total number of ions. */
    int interactions;   /*!<\brief Flag for the interaction algorithm. */
    int search;         /*!<\brief Flag for the stepsize algorithm. */
    int seed;           /*!<\brief Seed for the initial (random) conditions. */
    double stdo[3];     /*!<\brief Array of stiffness of the potential. */
    int *n;             /*!<\brief Ions' number in each species array. */
    double *m;          /*!<\brief Species' mass array. */
    double *eff;        /*!<\brief Array of effective stiffness. */
    double dpres;       /*!<\brief Target precision on displacement. */
    double gpres;       /*!<\brief Target precision on step size. */
    bool file;          /*!<\brief Monitoring active (or not). */
    string monitor;     /*!<\brief Name of the output monitor file. */
    string save_file;   /*!<\brief Name of the output positions file. */
    string init;        /*!<\brief Name of an input positions file. */
    bool *cDir;         /*!<\brief Array of cooling directions. */
    double *cInt;       /*!<\brief Array of cooling intensities. */
};
/* }}} */
/*! \brief This method contains the main code. */
int myFunction(ConfigMap &);
/*! \brief This method initializes the ions positions. */
bool initialization(Options &,double *);
/*! \brief This method initializes the options. */
void getOptions(Options &, ConfigMap &);
/*! \brief This method displays the ions position on the standard output. */
void display(Options &,double *,double);
/*! \brief This method write the ions position in a file. */
int write_position_file(Options &,double *,double);
/*! \brief This method display synthetic information on the cloud. */
void cloud_analysis(Options &,double *,double);
/*! \brief This method implements the gradient descent algorithm. */
double gradient_descent(Options &,double *,double *);
/*! \brief This method implements an algorithm to find the stepsize. */
int linear_search(Options &,double *,double *,double,double,double &);
/*! \brief This method computes the energy and the gradient of the system. */
int energy_grad(Options &,double *,double *,double &,double &);
/*! \brief This method computes the first and second derivatives along the
 * gradient. */
bool guess(Options &,double *,double *,double &);
/*! \brief This method permuts the ions positions to minimize the energy. */
void shuffle(Options &,double *,double *);
void test(Options &,double *);
#endif
/* mycode.h */
