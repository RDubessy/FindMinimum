/* This file is a part of findMinimum. {{{
 * Copyright (C) 2010 Romain Dubessy
 *
 * findMinimum is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *
 * findMinimum is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with findMinimum.  If not, see <http://www.gnu.org/licenses/>.
 *
 * }}} */
/*! \mainpage  findMinimum documentation.
 *
 * \section Installation
 * Download the latest version of this code and extract the archive.
 * In the main folder simply issue the command:
 * \code
 * make all
 * \endcode
 * If you want this program to be installed in your /usr/local/bin directory, 
 * use (as root):
 * \code
 * make install
 * \endcode
 * (of course this directory must exists).
 *
 * See the src/Makefile file for documentation on the various possible
 * optimizations.
 *
 * \section Usage
 * Usage information may be obtained issuing:
 * \code
 * findMinimum --usage
 * \endcode
 * The generic usage of findMinimum would be:
 * \code
 * findMinimum options [--key=value] file
 * \endcode
 * Where options is one of the arguments detailled hereafter and file is a 
 * configuration file.
 * The arguments may be passed inside the configuration file or as command line
 * options and in this case they override the contents of the 
 * file.
 * See the examples to find out how the configuration file is formatted.
 * The arguments may be any of those:
 * - algorithm::interactions, integer value, 0 for full, 1 for Barnes-Hut;
 * - algorithm::search, integer value, 0 for exact, 1 for backtracking;
 * - algorithm::monitor, string, name of the monitoring file;
 * - algorithm::dpres, double, target precision on displacement;
 * - algorithm::gpres, double target precision on step size;
 * - trap::ox, double, stiffness (in MHz) along the x axis, for Hydrogen;
 * - trap::oy, double, stiffness (in MHz) along the y axis, for Hydrogen;
 * - trap::oz, double, stiffness (in MHz) along the z axis (mass independant);
 * - ions::size, double, diameter of the initial cloud (in microns);
 * - ions::nEsp, integer, number of ion species (with different mass);
 * - ions::n1, integer, number of ions of the first specie;
 * - ions::m1, double, mass of the first specie (in proton unit mass);
 * - ions::n2, ions::m2, and so on...
 *
 * \section Output
 * The output of the program consists in three kind of data : monitoring of the
 * algorithm, ions' positions and warning.
 *
 * \subsection Monitoring
 * The monitoring output is written to the file passed by the option
 * \code--algorithm::monitor=name\endcode.
 * This file consists in three columns:
 * - the step number (integer);
 * - the current energy (double);
 * - the last step magnitude (double).
 *
 * \subsection Positions
 * The ions' positions are written on the program standard output once in a
 * while, when the relative energy variation is small enough.
 * 
 * \subsection Warnings
 * The warnings are displayed on the program standard error output.
 * Small ones are labelled by a <code>[W]</code> and more serious ones are
 * labelled by a <code>[E]</code>.
 */
#include <iostream>
#include "common.h"
#include "mycode.h"
int main(int argc, char *argv[]) {
    ConfigMap config;
    if(!parseOptions(argc,argv,config)) {        //Parse cmd line options
        std::cerr << "==> Try '" << argv[0] << " --usage'" << std::endl;
        return -1;
    }
    if(config["usage"].size()>0||config["help"].size()>0) {
        printUsage(argv[0]);
        return -1;
    }
    if(!parseConfig(config)) {
        std::cerr << "==> Try 'man " << argv[0] << "'" << std::endl;
        return -1;
    }
    return myFunction(config);
}
/* main.cpp */
