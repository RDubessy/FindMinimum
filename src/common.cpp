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
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "common.h"
/* readString: {{{ */
bool readString(string s, ConfigMap &config, string &appKey, bool verbose) {
    bool res=true;
    if(s.size()==0)                              //Empty line, skipping
        return true;
    if(s[0]=='[') {                              //Marker
        int index=s.find(']');
        switch(index) {
            case -1:
            case 1:
                cerr << "[E] Bad marker : '" << s << "' !" << endl;
                res=false;
                break;
            default:
                appKey=s.substr(1,index-1);
                if(verbose)
                    cerr << "[I] Found marker : '" << appKey << "'." << endl;
                appKey+="::";
        }
    }
    if(s[0]=='<') {                              //Option
        int indexA=s.find(',');
        int indexB=s.find('=');
        if(indexA>indexB)
            res=false;
        switch(indexA) {
            case -1:
            case 1:
                res=false;
                break;
            default:
                break;
        }
        switch(indexB) {
            case -1:
            case 2:
                res=false;
                break;
            default:
                break;
        }
        if(res) {
            string key=s.substr(1,indexA-1);
            if(verbose)
                cerr << "[I] Found option : '" << key << "'." << endl;
            string value=s.substr(indexB+1);
            if(config[appKey+key].size()>0) {
                if(verbose)
                    cerr << "[I] Previous value kept : '" << appKey+key << "="
                        << config[appKey+key] << "'." << endl;
            } else {
                config[appKey+key]=value;
                if(verbose)
                    cerr << "[I] '" << appKey+key << "=" << value << "'."
                        << endl;
            }
        } else {
            cerr << "[E] Bad option : '" << s << "' !" << endl;
        }
    }
    return res;
}
/* }}} */
/* parseOptions: {{{ */
bool parseOptions(const int argc, char *argv[], ConfigMap &config) {
    bool res=true;
    if(argc==1) {
        cerr << "[W] Need at least one argument !" << endl;
        res=false;
    }
    int index;
    for(int i=1;i<argc;i++) {
        string s(argv[i]);
        if(s.size()==1) {
            cerr << "[E] option too short !" << endl;
        } else if(s[0]=='-') {
            switch(s[1]) {
                case '-':                        //Long style option
                    s=s.substr(2);               //Remove '--'
                    index=s.find('=');
                    switch(index) {
                        case -1:
                            if(s.compare("usage")||s.compare("help"))
                                config[s]=" ";
                            break;
                        case 0:
                            cerr << "[W] Bad cmd line option : ' --" << s
                                << "' ! [ignored]" << endl;
                            res=false;
                            break;
                        default:
                            string key=s.substr(0,index);
                            string value=s.substr(index+1);
                            config[key]=value;
                    }
                    break;
                case 'v':
                    config["general::verbose"]="yes";
                    break;
                default:
                    cerr << "[W] Bad cmd line option : ' --" << s
                        << "' ! [ignored]" << endl;
                    res=false;
            }
        }
        else {                                   //Valid config file ?
            if(config["configFile"].size()>0) {  //Conflict !
                cerr << "[W] Conflicting config file names !" << endl;
                res=false;
            }
            else {
                config["configFile"]=s;
            }
        }
    }
    return res;
}
/* }}} */
/* printUsage: {{{ */
void printUsage(char argv[]) {
    cerr << "More information at www-link\n"
        << "Usage :\n"
        << "%" << argv << " options filename\n"
        << "Where 'filename' is the name of a valid config file\n"
        << "and 'options' stand for :\n"
        << "\t--usage, --help : display this screen and exit."
        << endl;
    return;
}
/* }}} */
/* strip: {{{ */
string strip(string s) {
    string res("");
    int n=s.size();
    for(int i=0;i<n;i++) {
        if(s[i]=='#')
            break;
        if(s[i]!=' ')
            res+=s[i];
    }
    return res;
}
/* }}} */
/* parseConfig: {{{ */
bool parseConfig(ConfigMap &config) {
    string appKey("");
    ifstream file(config["configFile"].c_str(),ios::in);
    if(!file.good()) {
        cerr << "[E] Error opening the configuration file : '" 
            << config["configFile"] << "' !" << endl;
        return false;
    }
    bool res=true;
    bool verbose=false;
    if(config["general::verbose"].size()>0)
        verbose=true;
    while(file) {
        string s;
        getline(file,s);
        s=strip(s);
        res=res && readString(s,config,appKey,verbose);
    }
    file.close();
    return res;
}
/* }}} */
/* getConfig: {{{ */
int getConfig(ConfigMap &config, const string &name, int def) {
    if(config[name].size()>0)
        return atoi(config[name].c_str());
    cerr << "[W] Key : '" << name.c_str() 
        << "' not found, using default value : " << def << endl;
    return def;
}
double getConfig(ConfigMap &config, const string &name, double def) {
    if(config[name].size()>0)
        return atof(config[name].c_str());
    cerr << "[W] Key : '" << name.c_str() 
        << "' not found, using default value : "
        << def << endl;
    return def;
}
string getConfig(ConfigMap &config, const string &name, const string &def) {
    if(config[name].size()>0)
        return config[name];
    cerr << "[W] Key : '" << name.c_str() 
        << "' not found, using default value : "
        << def.c_str() << endl;
    return def;
}
/* }}} */
/* common.cpp */
