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
#include "mycode.h"
/*! \file */
/*!\brief Maximum number of iterations in the linear search loop. */
#define LSTOP 100
/* Usefull constants. {{{ */
double eps0=8.854187817e-12;           /*!<\brief Vacuum permittivity. */
double e=1.602176487e-19;              /*!<\brief Elementary charge. */
double pi=3.141592653589793;           /*!<\brief Value of pi (16 digits). */
double mp=1.672621637e-27;             /*!<\brief Proton mass. */
double h=6.62606896e-34;               /*!<\brief Planck constant. */
/*!\brief Energy normalization coefficient. */
double coeff=8*pi*pi*pi*eps0*mp/(e*e)*1e-6;
/*!\brief Cooling intensity normalization coefficient. */
double coeffCool=4*pi*eps0*h/(e*e)*1e3;
/* }}} */
/* toString: translate to string an object {{{ */
/*! Usefull method to translate an object into a standard string using the
 * stringstream class.
 */
template <class T>
inline std::string toString(const T &t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
};
/* }}} */
/* mySwap: swap two objects of same type {{{ */
/*! Usefull method to swap two objects (not pointers !) and save some code
 * lines.
 */
template <class T>
inline void mySwap(T &a, T &b) {
    T tmp=a;
    a=b;
    b=tmp;
    return;
};
/* }}} */
/* shuffle: permuts the ions to try and minimize the energy. {{{ */
/*! This method tries and minimize the system energy by permuting the ions on
 * their sites.
 * For each ion it randomly choose an ion with a different mass and whenever it
 * decreases the system energy, swap the two ions positions.
 */
void shuffle(Options &options, double *ions) {
    if(options.nEsp==1)
        return;
    int max=0;
    int c=0;
    int n=options.nTot;
    for(int i=0;i<options.nEsp;i++) {
        max+=options.n[i];
        while(c<max) {
            //Randomly choose on ion of an other mass.
            int j=((rand()%(n-options.n[i]))+max)%n;
            double Eact=options.eff[3*c]*ions[3*c]*ions[3*c]
                +options.eff[3*c+1]*ions[3*c+1]*ions[3*c+1]
                +options.eff[3*c+2]*ions[3*c+2]*ions[3*c+2]
                +options.eff[3*j]*ions[3*j]*ions[3*j]
                +options.eff[3*j+1]*ions[3*j+1]*ions[3*j+1]
                +options.eff[3*j+2]*ions[3*j+2]*ions[3*j+2];
            Eact*=0.5;
            for(int k=0;k<3;k++) {
                if(options.cDir[3*j+k])
                    Eact-=options.cInt[3*j+k]*ions[3*j+k];
                if(options.cDir[3*c+k])
                    Eact-=options.cInt[3*c+k]*ions[3*c+k];
            }
            double Etest=options.eff[3*j]*ions[3*c]*ions[3*c]
                +options.eff[3*j+1]*ions[3*c+1]*ions[3*c+1]
                +options.eff[3*j+2]*ions[3*c+2]*ions[3*c+2]
                +options.eff[3*c]*ions[3*j]*ions[3*j]
                +options.eff[3*c+1]*ions[3*j+1]*ions[3*j+1]
                +options.eff[3*c+2]*ions[3*j+2]*ions[3*j+2];
            Etest*=0.5;
            for(int k=0;k<3;k++) {
                if(options.cDir[3*j+k])
                    Etest-=options.cInt[3*j+k]*ions[3*c+k];
                if(options.cDir[3*c+k])
                    Etest-=options.cInt[3*c+k]*ions[3*j+k];
            }
            if(Etest<Eact) {                     //Swap ions
                mySwap(ions[3*c],ions[3*j]);
                mySwap(ions[3*c+1],ions[3*j+1]);
                mySwap(ions[3*c+2],ions[3*j+2]);
            }
            c++;
        }
    }
    return;
}
/* }}} */ 
/*******************************************************************************
 * Implements the gradient descent method, with two alternative stepsize choice
 * algorithm : backtracking and exact line seach.
 ******************************************************************************/
/* energy_grad: computes the system energy and its gradient. {{{ */
int energy_grad(Options &options, double *ions, double *grad, double &epot,
        double &eint) {
    int n=options.nTot;
    int res=0;
    epot=0.0;
    for(int i=0;i<3*n;i++) {                 //External
        epot+=0.5*options.eff[i]*ions[i]*ions[i];
        grad[i]=options.eff[i]*ions[i];
        if(options.cDir[i]) {
            epot-=options.cInt[i]*ions[i];
            grad[i]-=options.cInt[i];
        }
    }
    eint=0.0;
    if(options.interactions==0) {
        /* Full interactions {{{ */
        for(int i=0;i<n;i++) {
            int ii=3*i;
            for(int j=i+1;j<n;j++) {
                int jj=3*j;
                double invr=(ions[ii]-ions[jj])*(ions[ii]-ions[jj])
                    +(ions[ii+1]-ions[jj+1])*(ions[ii+1]-ions[jj+1])
                    +(ions[ii+2]-ions[jj+2])*(ions[ii+2]-ions[jj+2]);
                invr=1.0/sqrt(invr);
                eint+=invr;
                invr=pow(invr,3.0);
                double tmp=invr*(ions[ii]-ions[jj]);
                grad[ii]-=tmp;
                grad[jj]+=tmp;
                ii++;
                jj++;
                tmp=invr*(ions[ii]-ions[jj]);
                grad[ii]-=tmp;
                grad[jj]+=tmp;
                ii++;
                jj++;
                tmp=invr*(ions[ii]-ions[jj]);
                grad[ii]-=tmp;
                grad[jj]+=tmp;
                ii-=2;
                res++;
            }
        }
        /* }}} */
    } else if(options.interactions==1) {
        octTree tree(n,ions,grad);
        res=tree.energy_grad(grad,eint);
    }
    return res;
}
/* }}} */
/* guess: computes the estimate for the minimum position. {{{ */
bool guess(Options &options, double *ions, double *grad, double &up) {
    int n=options.nTot;
    double fp=0.0;
    double fpp=0.0;
    for(int i=0;i<3*n;i++) {
        fp+=options.eff[i]*ions[i]*grad[i];
        fpp+=options.eff[i]*grad[i]*grad[i];
        if(options.cDir[i])
            fp-=options.cInt[i]*grad[i];
    }
    if(options.interactions==0) {
        for(int i=0;i<n;i++) {
            int ii=3*i;
            /* Full interactions {{{ */
            for(int j=i+1;j<n;j++) {
                int jj=3*j;
                double r=ions[ii]-ions[jj];
                double g=grad[jj]-grad[ii];
                double invr=r*r;
                double rg=r*g;
                double gg=g*g;
                ii++;
                jj++;
                r=ions[ii]-ions[jj];
                g=grad[jj]-grad[ii];
                invr+=r*r;
                rg+=r*g;
                gg+=g*g;
                ii++;
                jj++;
                r=ions[ii]-ions[jj];
                g=grad[jj]-grad[ii];
                invr+=r*r;
                rg+=r*g;
                gg+=g*g;
                ii-=2;
                invr=1.0/sqrt(invr);
                fp+=pow(invr,3.0)*rg;
                fpp+=pow(invr,3.0)*(3*invr*invr*rg*rg-gg);
            }
            /* }}} */
        }
    } else if(options.interactions==1) {
        octTree tree(n,ions,grad);
        tree.guess(fp,fpp);
    }
    up=fp/fpp;
    if(fpp<0)
        return false;
    return true;
}
/* }}} */
/* linear_search: computes the step size in the optimization loop. {{{ */
int linear_search(Options &options, double *ions, double *grad,
        double epot, double eint, double &res) {
    int n=options.nTot;
    double up=0.0;
    double h=0.0;
    double *np=new double[3*n];
    int c=0;
    double fp,fpp;
    fp=fpp=1.0;
    do {
        for(int i=0;i<3*n;i++) {
            np[i]=ions[i]+h*grad[i];
        }
        if(guess(options,np,grad,up)) {
            h-=up;
            c++;
        } else {                                           //Dichotomie
            cerr << "[I] Switching to dichotomie" << endl;
            h=-2.0;
            double Enew,Enewb;
            Enew=Enewb=0;
            do {
                h/=2;
                for(int i=0;i<3*n;i++)
                    np[i]=ions[i]+h*grad[i];
                energy_grad(options,np,grad,Enew,Enewb);
            } while(Enew+Enewb>epot+eint);
            break;
        }
    } while(fabs(up)>options.gpres && c<LSTOP);
    if(c==LSTOP) {
        cerr << "[I] Linear search converging slowly, "
            << "computing the gradient at current position."
            << endl;
    }
    res=0.0;
    //double tmp=0.0;
    for(int i=0;i<3*n;i++) {
        ions[i]=np[i];
        res+=grad[i]*grad[i];
        //tmp=grad[i]*grad[i];
        //if(tmp>res)
        //    res=tmp;
    }
    delete[] np;
    res=sqrt(res)*fabs(h);
    return c;
}
/* }}} */
/* gradient_descent: implements the gradient descent algorithm. {{{ */
/*! This method implements the gradient descent algorithm.
*/
double gradient_descent(Options &options, double *ions) {
    int n=options.nTot;
    double *grad=new double[3*n];
    double diff=1e10;
    double epot=0;
    double eint=0;
    double pres=1.0;
    double Epres=1.0;
    double norm,tmp,max;
    int c=0;                                     //Loop counter
    int nc=0;
    int ns=0;
    ofstream file(options.monitor.c_str());
    bool monitor=true;
    if(!file.good()) {
        cerr << "[E] Can not open monitoring file (" << options.monitor.c_str()
            << ") !" << endl;
        monitor=false;
    }
    norm=max=tmp=0.0;
    nc=energy_grad(options,ions,grad,epot,eint); //Compute initial energy & grad
    if(monitor)
        file << "#step Energy[pot] Energy[int] displacement Nint Nsearch\n"
            << c << " " << epot << " " << eint << " " << max << " " << nc
            << " " << ns << endl;;
    display(options,ions,epot+eint);
    cerr << "[I] Initial energy : " << epot + eint 
        << "[" << epot << "|" << eint<<"]" << endl;
    while(true) {
        if(options.search==0) {
            /* Exact line search. {{{ */
            ns=linear_search(options,ions,grad,epot,eint,max);
            diff=epot+eint;
            shuffle(options,ions);
            nc=energy_grad(options,ions,grad,epot,eint);
            diff-=epot+eint;
            /* }}} */
        } else if(options.search==1) {
            /* Backtracking search. {{{ */
            double t=1.0;
            double beta=0.5;
            norm=0.0;
            max=0.0;
            for(int i=0;i<3*n;i++) {
                tmp=grad[i]*grad[i];
                norm+=tmp;
                if(tmp>max)
                    max=tmp;
            }
            diff=epot;
            double *np=new double[3*n];
            do {
                t*=beta;
                //stop=diff-alpha*t*norm;
                for(int i=0;i<3*n;i++) {
                    np[i]=ions[i]-t*grad[i];
                }
                shuffle(options,ions);
                nc=energy_grad(options,np,grad,epot,eint);
            } while(epot>diff);
            for(int i=0;i<3*n;i++)
                ions[i]=np[i];
            delete[] np;
            diff-=epot;
            max=sqrt(max)*t;
            norm=sqrt(norm)*t/n;
            /* }}} */
        }
        if(diff<0)
            cerr << "[I] Warning : Energy incresead at this step ["
                << c << "] !\n"
                << "[I] Was : " << epot+diff << ", is : " << epot << endl;
        else {
            if(c%10==0)
                display(options,ions,epot+eint);
            /*
            diff/=epot+eint;
            if(diff<Epres) {
                cerr << "[I] At step:" << c << ", Delta E/E was below " 
                    << Epres << " (" << diff << ") for the first time."
                    << endl;
                Epres/=10;
                display(options,ions,epot+eint);
            }
            if(diff==0.0)
                break;
                */
        }
        if(max<pres) {
            cerr << "[I] precision below: " << pres << " (" << max 
                << ") in : " << c << " iterations." << endl;
            pres/=10.0;
        }
        c++;
        if(monitor)
            file << c << " " << epot << " " << eint << " " << max << " " 
                << nc << " " << ns << endl;
        if(pres<options.dpres)
            break;
    }
    cerr << "[I] Gradient descent took " << c << " iterations to complete."
        << endl;
    display(options,ions,epot+eint);
    if(monitor) {
        file << endl;
        file.close();
    }
    delete[] grad;
    return epot;
}
/* }}} */
/*******************************************************************************
 * Miscellaneous functions.
 ******************************************************************************/
/* display: displays the ions coordinates on the standard output. {{{ */
void display(Options &options, double *ions, double frame) {
    int c=0;
    int max=0;
    cout << "#frame=" << frame << "\n";
    for(int i=0;i<options.nEsp;i++) {
        max+=options.n[i];
        cout << "#m=" << options.m[i] << "\n";
        while(c<max) {
            cout << ions[3*c] << " " << ions[3*c+1] << " " << ions[3*c+2] 
                << "\n";
            c++;
        }
        cout << "\n";
    }
    cout << endl;
    return;
}
/* }}} */
/* cloud_analysis: {{{ */
void cloud_analysis(Options &options, double *ions, double epot) {
    int n=options.nTot;
    double x,y,z,x2,y2,z2;                       //Ensemble average values
    x=y=z=x2=y2=z2=0.0;
    for(int i=0;i<n;i++) {
        int ii=3*i;
        x+=ions[ii];
        x2+=ions[ii]*ions[ii];
        ii++;
        y+=ions[ii];
        y2+=ions[ii]*ions[ii];
        ii++;
        z+=ions[ii];
        z2+=ions[ii]*ions[ii];
    }
    x/=n;
    y/=n;
    z/=n;
    x2/=n;
    y2/=n;
    z2/=n;
    x2-=x*x;
    y2-=y*y;
    z2-=z*z;
    x2=sqrt(x2);
    y2=sqrt(y2);
    z2=sqrt(z2);
    cerr << "************************************************************\n"
        << "Cloud Analysis\n";
    cerr << "Energy : " << epot << "\n";
    cerr << "First order :\n"
        << "\tMean positions : " << x << "," << y << "," << z << "\n"
        << "\tRms radius     : " << x2 << "," << y2 << "," << z2 << "\n";
    cerr << "End of analysis\n"
        << "************************************************************"
        << endl;
    return;
}
/* }}} */
/* getOptions: parse the options. {{{ */
void getOptions(Options &options, ConfigMap &config) {
    //Ions
    options.nEsp=getConfig(config,"ions::nEsp",1);
    options.size=getConfig(config,"ions::size",2e3);
    options.n=new int[options.nEsp];
    options.m=new double[options.nEsp];
    int *tmpDir=new int[options.nEsp];
    double *tmpInt=new double[options.nEsp];
    options.nTot=0;
    for(int i=0;i<options.nEsp;i++) {
        string tmp;
        tmp=string("ions::n")+toString(i+1);
        options.n[i]=getConfig(config,tmp.c_str(),1);
        options.nTot+=options.n[i];
        tmp=string("ions::m")+toString(i+1);
        options.m[i]=getConfig(config,tmp.c_str(),1.0);
        tmp=string("ions::c")+toString(i+1);
        tmpInt[i]=getConfig(config,tmp.c_str(),0.0);
        tmp=string("ions::d")+toString(i+1);
        tmpDir[i]=getConfig(config,tmp.c_str(),-1);
    }
    //Trap
    options.stdo[0]=getConfig(config,"trap::ox",1.0);
    options.stdo[1]=getConfig(config,"trap::oy",1.0);
    options.stdo[2]=getConfig(config,"trap::oz",1.0);
    options.eff=new double[3*options.nTot];
    options.cDir=new bool[3*options.nTot];
    options.cInt=new double[3*options.nTot];
    int c=0;
    int max=0;
    for(int i=0;i<options.nEsp;i++) {
        max+=options.n[i];
        while(c<max) {
            int ii=3*c;
            options.cDir[ii]=options.cDir[ii+1]=options.cDir[ii+2]=false;
            options.eff[ii]=coeff*pow(options.stdo[0],2.0)/options.m[i];
            options.eff[ii+1]=coeff*pow(options.stdo[1],2.0)/options.m[i];
            options.eff[ii+2]=coeff*options.stdo[2];
            if(tmpDir[i]>-1) {
                options.cDir[ii+tmpDir[i]]=true;
                options.cInt[ii+tmpDir[i]]=coeffCool*tmpInt[i];
            }
            c++;
        }
    }
    //Algortihm
    options.interactions=getConfig(config,"algorithm::interactions",0);
    options.search=getConfig(config,"algorithm::search",0);
    options.monitor=getConfig(config,"algorithm::monitor","monitor");
    options.dpres=getConfig(config,"algorithm::dpres",1e-3);
    options.gpres=getConfig(config,"algorithm::gpres",1e-3);
    options.seed=getConfig(config,"algorithm::seed",-1);
    if(options.seed==-1)
        options.seed=(int)time(NULL);
    //General
    options.file=false;
    options.init=getConfig(config,"general::init","");
    if(options.init.size()>0)
        options.file=true;
    return;
}
/* }}} */
/* initialization: initialize the system. {{{ */
bool initialization(Options &options, double *ions) {
    int n=options.nTot;
    if(options.file) {
        ifstream file(options.init.c_str(),ios::in);
        int l=0;
        while(file!=0) {
            string s;
            getline(file,s);
            if(s.size()>0 && s[0]!='#') {
                if(l==n) {
                    cerr << "[E] Error wrong number of particules !" << endl;
                    return false;
                }
                int i=s.find(' ');
                int ll=3*l;
                ions[ll]=atof(s.substr(0,i-1).c_str());
                s=s.substr(i+1);
                i=s.find(' ');
                ions[ll+1]=atof(s.substr(0,i-1).c_str());
                ions[ll+2]=atof(s.substr(i+1).c_str());
                l++;
            }
        }
        file.close();
    } else {
        srand(options.seed);                     //Initialize random generator
        for(int i=0;i<3*n;i+=3) {                //Position initialization
            //ions[i]=options.size*((double)rand()/RAND_MAX-0.5);
            double r=options.size*pow((double)rand()/RAND_MAX,1.0/3.0);
            double phi=2*pi*((double)rand()/RAND_MAX);
            double theta=acos(2*((double)rand()/RAND_MAX-0.5));
            ions[i]=r*sin(theta)*cos(phi);
            ions[i+1]=r*sin(theta)*sin(phi);
            ions[i+2]=r*cos(theta);
        }
    }
    return true;
}
/* }}} */
/* myFunction: main code. {{{ */
int myFunction(ConfigMap &config) {
    /* Insert here your code or call function and include the headers. */
    /* To retrieve the value of a given parameter, say an integer 'n', defined
     * in the [general] section of the configuration file under the name 'n' use
     * something as :
     * int n=getConfig(config,string("general::n"),1);
     * double n=getConfig(config,string("general::n"),1.0);
     * string n=getConfig(config,string("general::n"),string("foo"));
     */
    /* The return value of this function is the return value of the program. */
    /* Try and find equilibrium positions of a system of point like charges in a
     * potential well. */
    Options options;
    getOptions(options,config);
    int n=options.nTot;
    double *ions=new double[3*n];
    if(initialization(options,ions)) {
        if(config["bench"].size()>0)
            test(options,ions);
        else {
            double epot=gradient_descent(options,ions);
            cloud_analysis(options,ions,epot);   //Display synthetic summary
        }
    }
    delete[] ions;
    delete[] options.n;
    delete[] options.m;
    delete[] options.eff;
    return 0;
}
/* }}} */
/* test: {{{ */
/*! This method is here for testing purposes.
*/
void test(Options &options, double *ions) {
    int n=options.nTot;
    double *grad=new double[3*n];
    double epot,fp,fpp;
    epot=0;
    octTree tree(n,ions,grad);
    for(int i=0;i<3*n;i++)
        grad[i]=0;
    for(int i=0;i<100;i++) {
        tree.energy_grad(grad,epot);
        tree.guess(fp,fpp);
    }
    cerr << epot << endl;
    delete[] grad;
}
/* }}} */
/* mycode.cpp */
