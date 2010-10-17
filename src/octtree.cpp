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
#include <fstream>
#include <iostream>
#include <math.h>
#ifdef VECTORIZE
#include <xmmintrin.h>
#endif //VECTORIZE
#include "octtree.h"
using namespace std;
/* octTree: default constructor.  {{{ */
octTree::octTree(void) {
#ifdef VECTORIZE
    g2mean=delta=0.0f;
#else
    g2mean=delta=0.0;
#endif //VECTORIZE
    n_q=level=0;
    child=skip=next=0;
    indice=-1;
}
/* }}} */
/* octTree: constructs a tree from an array of charges. {{{ */
/*! This constructor takes an array of charges position and insert it into the
 * tree, using the "insert" method.
 */
octTree::octTree(int n, double *pos, double *grad) {
#ifdef VECTORIZE
    g2mean=delta=0.0f;
    /* Compute the size of the first cell. {{{ */
    size=0.0f;
    for(int i=0;i<n;i++) {
        float r=(float)sqrt(pos[3*i]*pos[3*i]+pos[3*i+1]*pos[3*i+1]
                +pos[3*i+2]*pos[3*i+2]);
        if(r>size)
            size=r;
    }
    size*=2;
    /* }}} */
    center[0]=center[1]=center[2]=0.0f;
#else
    g2mean=delta=0.0;
    /* Compute the size of the first cell. {{{ */
    size=0;
    for(int i=0;i<n;i++) {
        double r=sqrt(pos[3*i]*pos[3*i]+pos[3*i+1]*pos[3*i+1]
                +pos[3*i+2]*pos[3*i+2]);
        if(r>size)
            size=r;
    }
    size*=2;
    /* }}} */
    center[0]=center[1]=center[2]=0.0;
#endif //VECTORIZE
    n_q=level=0;
    child=skip=next=0;
    indice=-1;
    for(int i=0;i<n;i++)
        insert(&pos[3*i],&grad[3*i],i);
    flatten();
    finalize();
}
/* }}} */
/* ~octTree: destructor. {{{ */
octTree::~octTree(void) {
    if(child!=0)
        delete[] child;
}
/* }}} */
/* display: (flat). {{{ */
void octTree::display(void) {
    octTree *tmp=this;
    while (tmp!=0) {
        for(int i=0;i<tmp->level;i++)
            cerr << ".";
        cerr << tmp->level;
        cerr << "\t\t" << tmp << "\t" << tmp->indice << endl;
        tmp=tmp->next;
    }
    return;
}
/* }}} */
/* flatten: (flat). {{{ */
void octTree::flatten(void) {
    next=skip=0;
    if(n_q<2) {
        return;
    }
    octTree *tmp=this;
    do {
        if(tmp->n_q==1) {                        //One charge: next=skip.
            tmp->next=tmp->skip;
        } else {
            //Next is the first non empty child.
            int i=0;
            while(tmp->child[i].n_q==0)
                i++;
            tmp->next=&(tmp->child[i]);
            //Initialize children's skip pointer.
            for(int j=i+1;j<8;j++) {
                if(tmp->child[j].n_q!=0) {
                    tmp->child[i].skip=&(tmp->child[j]);
                    i=j;
                }
            }
            //The last non-empty child inherits its father skip pointer.
            tmp->child[i].skip=tmp->skip;
        }
        tmp=tmp->next;
    } while(tmp!=0);
    return;
}
/* }}} */
/* finalize: (flat). {{{ */
void octTree::finalize(void) {
    for(octTree *tmp=this;tmp!=0;tmp=tmp->next) {
        tmp->delta=tmp->g2mean=0.0f;
#ifdef VECTORIZE
        for(int i=0;i<3;i++) {
            ((float*)&(tmp->mean))[i]/=(float)tmp->n_q;
            ((float*)&(tmp->gmean))[i]/=(float)tmp->n_q;
            ((float*)&(tmp->gcmean))[i]/=(float)tmp->n_q;
            ((float*)&(tmp->gcmeanb))[i]/=(float)tmp->n_q;
        }
        for(int i=0;i<3;i++) {
            tmp->g2mean+=((float*)&(tmp->gcmean))[i];
            float tmpf=((float*)&(tmp->mean))[i]-tmp->center[i];
            tmp->delta+=tmpf*tmpf;
        }
#else
        for(int i=0;i<3;i++) {
            tmp->mean[i]/=(double)tmp->n_q;
            tmp->gmean[i]/=(double)tmp->n_q;
        }
        for(int i=0;i<6;i++) {
            tmp->gcmean[i]/=(double)tmp->n_q;
        }
        for(int i=0;i<3;i++) {
            tmp->g2mean+=tmp->gcmean[i];
            double tmpd=tmp->mean[i]-tmp->center[i];
            tmp->delta+=tmpd*tmpd;
        }
#endif //VECTORIZE
        tmp->delta=sqrt(tmp->delta);
        //Simpler open cell criterion
        tmp->delta+=tmp->size*TEST_CONST;
        tmp->delta*=tmp->delta;
        //Simpler gcmean[3:5]
#ifdef VECTORIZE
        for(int i=0;i<3;i++) {
            ((float*)&(tmp->gcmeanb))[i]*=2;
        }
#else
        tmp->gcmean[3]*=2;
        tmp->gcmean[4]*=2;
        tmp->gcmean[5]*=2;
#endif //VECTORIZE
    }
    return;
}
/* }}} */
#ifdef VECTORIZE
/* insert: {{{ */
/*! This method inserts a charge in the three, at the position stored in the
 * array "pos".
 * If the current cell is empty, the charge is inserted, and the node parameters
 * are updated accordingly.
 * If the current cell contains exactly one charge, the cell is splitted in
 * octants and the two charges are succesively inserted in the relevant
 * children.
 * If the current cell contains more than one charge, the charge is inserted in
 * the relevant octant.
 * In the last two cases, the node parameters are updated: charge increased by
 * one and mean position modified.
 *
 * Note that this function is recursive.
 *
 * The octree also stores information on the local gradient, provided by the
 * array "grad", see the optimization algorithm explaination for more details.
 */
void octTree::insert(double *pos, double *grad, int ind) {
    /* Empty leaf: inserting. {{{ */
    if(n_q==0) {
        n_q=1;
        indice=ind;
        for(int i=0;i<3;i++) {
            ((float*)&mean)[i]=(float)pos[i];
            ((float*)&gmean)[i]=(float)grad[i];
        }
        ((float*)&mean)[3]=0.0f;
        ((float*)&gmean)[3]=0.0f;
        gcmean=__builtin_ia32_mulps(gmean,gmean);
        ((float*)&gcmeanb)[0]=(float)(grad[0]*grad[1]);
        ((float*)&gcmeanb)[1]=(float)(grad[0]*grad[2]);
        ((float*)&gcmeanb)[2]=(float)(grad[1]*grad[2]);
        ((float*)&gcmeanb)[3]=0.0f;
        return;
    }
    /* }}} */
    /* Leaf with one charge. {{{ */
    if(n_q==1) {
        /* Creating the children. {{{ */
        child=new octTree[8];
        float halfsize=size/2;
        for(int i=0;i<8;i++) {
            child[i].child=0;
            child[i].size=halfsize;
            child[i].level=level+1;
        }
        halfsize/=2;
        child[0].center[0]=center[0]-halfsize;
        child[0].center[1]=center[1]-halfsize;
        child[0].center[2]=center[2]-halfsize;
        child[1].center[0]=center[0]-halfsize;
        child[1].center[1]=center[1]-halfsize;
        child[1].center[2]=center[2]+halfsize;
        child[2].center[0]=center[0]-halfsize;
        child[2].center[1]=center[1]+halfsize;
        child[2].center[2]=center[2]-halfsize;
        child[3].center[0]=center[0]-halfsize;
        child[3].center[1]=center[1]+halfsize;
        child[3].center[2]=center[2]+halfsize;
        child[4].center[0]=center[0]+halfsize;
        child[4].center[1]=center[1]-halfsize;
        child[4].center[2]=center[2]-halfsize;
        child[5].center[0]=center[0]+halfsize;
        child[5].center[1]=center[1]-halfsize;
        child[5].center[2]=center[2]+halfsize;
        child[6].center[0]=center[0]+halfsize;
        child[6].center[1]=center[1]+halfsize;
        child[6].center[2]=center[2]-halfsize;
        child[7].center[0]=center[0]+halfsize;
        child[7].center[1]=center[1]+halfsize;
        child[7].center[2]=center[2]+halfsize;
        /* }}} */
        /* Displace old charge. {{{ */
        double tmean[3];
        double tgmean[3];
        for(int i=0;i<3;i++) {
            tmean[i]=(double)((float*)&mean)[i];
            tgmean[i]=(double)((float*)&gmean)[i];
        }
        if(tmean[0]<center[0]) {       //Node 0-3
            if(tmean[1]<center[1]) {   //Node 0-1
                if(tmean[2]<center[2])
                    child[0].insert(&tmean[0],&tgmean[0],indice);
                else
                    child[1].insert(&tmean[0],&tgmean[0],indice);
            }
            else {                               //Node 2-3
                if(tmean[2]<center[2])
                    child[2].insert(&tmean[0],&tgmean[0],indice);
                else
                    child[3].insert(&tmean[0],&tgmean[0],indice);
            }
        }
        else {                                   //Node 4-7
            if(tmean[1]<center[1]) {   //Node 4-5
                if(tmean[2]<center[2])
                    child[4].insert(&tmean[0],&tgmean[0],indice);
                else
                    child[5].insert(&tmean[0],&tgmean[0],indice);
            }
            else {                               //Node 6-7
                if(tmean[2]<center[2])
                    child[6].insert(&tmean[0],&tgmean[0],indice);
                else
                    child[7].insert(&tmean[0],&tgmean[0],indice);
            }
        }
        /* }}} */
    }
    /* }}} */
    /* Update mean position. {{{ */
    for(int i=0;i<3;i++) {
        ((float*)&mean)[i]+=(float)pos[i];
        ((float*)&gmean)[i]+=(float)grad[i];
        ((float*)&gcmean)[i]+=(float)(grad[i]*grad[i]);
    }
    ((float*)&mean)[3]=0.0f;
    ((float*)&gmean)[3]=0.0f;
    ((float*)&gcmeanb)[0]+=(float)(grad[0]*grad[1]);
    ((float*)&gcmeanb)[1]+=(float)(grad[0]*grad[2]);
    ((float*)&gcmeanb)[2]+=(float)(grad[1]*grad[2]);
    ((float*)&gcmeanb)[3]=0.0f;
    n_q++;
    indice=-1;
    /* }}} */
    /* Insert new charge. {{{ */
    if(pos[0]<center[0]) {                      //Node 0-3
        if(pos[1]<center[1]) {                  //Node 0-1
            if(pos[2]<center[2])
                child[0].insert(pos,grad,ind);
            else
                child[1].insert(pos,grad,ind);
        }
        else {                                   //Node 2-3
            if(pos[2]<center[2])
                child[2].insert(pos,grad,ind);
            else
                child[3].insert(pos,grad,ind);
        }
    }
    else {                                       //Node 4-7
        if(pos[1]<center[1]) {                  //Node 4-5
            if(pos[2]<center[2])
                child[4].insert(pos,grad,ind);
            else
                child[5].insert(pos,grad,ind);
        }
        else {                                   //Node 6-7
            if(pos[2]<center[2])
                child[6].insert(pos,grad,ind);
            else
                child[7].insert(pos,grad,ind);
        }
    }
    /* }}} */
    return;
}
void octTree::insert(v4sf &pos, v4sf &grad, int ind) {
    n_q=1;
    indice=ind;
    mean=pos;
    gmean=grad;
    gcmean=__builtin_ia32_mulps(grad,grad);
    ((float*)&gcmeanb)[0]=((float*)&gmean)[0]*((float*)&gmean)[1];
    ((float*)&gcmeanb)[1]=((float*)&gmean)[0]*((float*)&gmean)[2];
    ((float*)&gcmeanb)[2]=((float*)&gmean)[1]*((float*)&gmean)[2];
    ((float*)&gcmeanb)[3]=0.0f;
    return;
}
/* }}} */
/* energy_grad: optimum (flat). {{{ */
/*!
 * More or less basic optimization performed:
 *  - prefetching tested but no improvement (why ?!?),
 *  - no left hand side assignement to array elements in the loop.
 *  TODO: read assembler output and try to find more optimizations...
 */
int octTree::energy_grad(double *grad, double &epot) {
    if(n_q<2) {                                  //Trivial case, nothing to do
        return 0;
    }
#ifdef DEBUG
    int res=0;
#endif
    int nint=0;
    float epot_t=0.0f;
    octTree *cur=this;
    do {                                         //Walk the tree
        if(cur->n_q==1) {
            int i=cur->indice;
            v4sf g=(v4sf){0.0f,0.0f,0.0f,0.0f};
            octTree *aux=this;
            do {
                int nq=aux->n_q;
                v4sf r=__builtin_ia32_subps(aux->mean,cur->mean);
                v4sf t=__builtin_ia32_mulps(r,r);//t.i=r.i*r.i
                float invr2=((float*)&t)[0]+((float*)&t)[1]+((float*)&t)[2];
                //Test cell
                if ((nq==1 && i!=aux->indice) || invr2>aux->delta) {
                    aux=aux->skip;
                    nint++;
                    invr2=1.0/invr2;
                    float invr=sqrt(invr2)*nq;
#ifdef DEBUG
                    res+=nq;
#endif
                    epot+=invr;
                    invr*=invr2;
                    t=(v4sf){invr,invr,invr,0};
                    r=__builtin_ia32_mulps(r,t);
                    g=__builtin_ia32_addps(g,r);
                } else {                         //Go deeper
                    aux=aux->next;
                }
            } while(aux!=0);
            int ii=i*3;
            grad[ii]+=((float*)&g)[0];
            grad[ii+1]+=((float*)&g)[1];
            grad[ii+2]+=((float*)&g)[2];
        }
        cur=cur->next;
    } while(cur!=0);
#ifdef DEBUG
    if(res!=n_q*(n_q-1)) {
        cerr << "[E] Error!" << endl;
        return -1;
    }
#endif
    epot+=(double)(epot_t*0.5f);
    return nint;
}
/* }}} */
/* guess: (flat). {{{ */
int octTree::guess(double &fp, double &fpp) {
    if(n_q<2) {                       //Trivial case, nothing to do
        return 0;
    }
    int res=0;
    octTree *cur=this;
    float fp_t,fpp_t;
    fp_t=fpp_t=0.0f;
    do {
        if(cur->n_q==1) {
            int i=cur->indice;
            octTree *aux=this;
            do {
                if(aux->n_q==1 && aux->indice!=i) {
                    v4sf r=__builtin_ia32_subps(cur->mean,aux->mean);
                    v4sf t=__builtin_ia32_mulps(r,r);//t.i=r.i*r.i
                    float invr2=((float*)&t)[0]+((float*)&t)[1]+((float*)&t)[2];
                    invr2=1.0/invr2;
                    float invr=sqrt(invr2);
                    float invr3=invr*invr2;
                    v4sf g=__builtin_ia32_subps(aux->gmean,cur->gmean);
                    t=__builtin_ia32_mulps(r,g);
                    float rg=((float*)&t)[0]+((float*)&t)[1]+((float*)&t)[2];
                    t=__builtin_ia32_mulps(g,g);
                    float gg=((float*)&t)[0]+((float*)&t)[1]+((float*)&t)[2];
                    fp+=invr3*rg;
                    fpp+=invr3*((rg*rg)*(3*invr2)-gg);
                    res++;
                    aux=aux->skip;
                } else {
                    v4sf r=__builtin_ia32_subps(cur->mean,aux->mean);
                    v4sf t=__builtin_ia32_mulps(r,r);//t.i=r.i*r.i
                    float invr2=((float*)&t)[0]+((float*)&t)[1]+((float*)&t)[2];
                    if(invr2>aux->delta) { //Test Cell
                        v4sf g=__builtin_ia32_subps(aux->gmean,cur->gmean);
                        t=__builtin_ia32_mulps(r,g);
                        float rg=((float*)&t)[0]+((float*)&t)[1]
                            +((float*)&t)[2];
                        t=__builtin_ia32_addps(g,aux->gmean);
                        t=__builtin_ia32_mulps(cur->gmean,t);
                        float gg=aux->g2mean;
                        gg-=((float*)&t)[0];
                        gg-=((float*)&t)[1];
                        gg-=((float*)&t)[2];
                        t=__builtin_ia32_mulps(r,cur->gmean);
                        float rg2=((float*)&t)[0]+((float*)&t)[1]
                            +((float*)&t)[2];
                        t=__builtin_ia32_mulps(r,aux->gmean);
                        rg2=rg2*(rg2-2*(((float*)&t)[0]+((float*)&t)[1]
                            +((float*)&t)[2]));
                        t=__builtin_ia32_mulps(r,r);
                        t=__builtin_ia32_mulps(t,aux->gcmean);
                        rg2+=((float*)&t)[0];
                        rg2+=((float*)&t)[1];
                        rg2+=((float*)&t)[2];
                        g=__builtin_ia32_shufps(r,r,_MM_SHUFFLE(0,1,2,3));//rx,rx,ry,0
                        t=__builtin_ia32_shufps(r,r,_MM_SHUFFLE(0,1,2,3));//ry,rz,rz,0
                        t=__builtin_ia32_mulps(t,aux->gcmeanb);
                        t=__builtin_ia32_mulps(t,g);
                        invr2=1.0/invr2;
                        float invr=sqrt(invr2);
                        float invr3=invr*invr2*aux->n_q;
                        fp+=invr3*rg;
                        fpp+=invr3*(3*invr2*rg2-gg);
                        res++;
                        aux=aux->skip;
                    } else {                         //Go deeper
                        aux=aux->next;
                    }
                }
            } while(aux!=0);
        }
        cur=cur->next;
    } while(cur!=0);
    fp/=2;
    fpp/=2;
    return res;
}
/* }}} */
#else
/* insert: {{{ */
/*! This method inserts a charge in the three, at the position stored in the
 * array "pos".
 * If the current cell is empty, the charge is inserted, and the node parameters
 * are updated accordingly.
 * If the current cell contains exactly one charge, the cell is splitted in
 * octants and the two charges are succesively inserted in the relevant
 * children.
 * If the current cell contains more than one charge, the charge is inserted in
 * the relevant octant.
 * In the last two cases, the node parameters are updated: charge increased by
 * one and mean position modified.
 *
 * Note that this function is recursive.
 *
 * The octree also stores information on the local gradient, provided by the
 * array "grad", see the optimization algorithm explaination for more details.
 */
void octTree::insert(double *pos, double *grad, int ind) {
    /* Empty leaf: inserting. {{{ */
    if(n_q==0) {
        n_q=1;
        mean[0]=pos[0];
        mean[1]=pos[1];
        mean[2]=pos[2];
        gmean[0]=grad[0];
        gmean[1]=grad[1];
        gmean[2]=grad[2];
        gcmean[0]=grad[0]*grad[0];
        gcmean[1]=grad[1]*grad[1];
        gcmean[2]=grad[2]*grad[2];
        gcmean[3]=grad[0]*grad[1];
        gcmean[4]=grad[0]*grad[2];
        gcmean[5]=grad[1]*grad[2];
        indice=ind;
        return;
    }
    /* }}} */
    /* Leaf with one charge. {{{ */
    if(n_q==1) {
        /* Creating the children. {{{ */
        child=new octTree[8];
        double halfsize=size/2;
        for(int i=0;i<8;i++) {
            child[i].size=halfsize;
            child[i].level=level+1;
        }
        halfsize/=2;
        child[0].center[0]=center[0]-halfsize;
        child[0].center[1]=center[1]-halfsize;
        child[0].center[2]=center[2]-halfsize;
        child[1].center[0]=center[0]-halfsize;
        child[1].center[1]=center[1]-halfsize;
        child[1].center[2]=center[2]+halfsize;
        child[2].center[0]=center[0]-halfsize;
        child[2].center[1]=center[1]+halfsize;
        child[2].center[2]=center[2]-halfsize;
        child[3].center[0]=center[0]-halfsize;
        child[3].center[1]=center[1]+halfsize;
        child[3].center[2]=center[2]+halfsize;
        child[4].center[0]=center[0]+halfsize;
        child[4].center[1]=center[1]-halfsize;
        child[4].center[2]=center[2]-halfsize;
        child[5].center[0]=center[0]+halfsize;
        child[5].center[1]=center[1]-halfsize;
        child[5].center[2]=center[2]+halfsize;
        child[6].center[0]=center[0]+halfsize;
        child[6].center[1]=center[1]+halfsize;
        child[6].center[2]=center[2]-halfsize;
        child[7].center[0]=center[0]+halfsize;
        child[7].center[1]=center[1]+halfsize;
        child[7].center[2]=center[2]+halfsize;
        /* }}} */
        /* Displace old charge. {{{ */
        if(mean[0]<center[0]) {                  //Node 0-3
            if(mean[1]<center[1]) {              //Node 0-1
                if(mean[2]<center[2])
                    child[0].insert(mean,gmean,indice);
                else
                    child[1].insert(mean,gmean,indice);
            }
            else {                               //Node 2-3
                if(mean[2]<center[2])
                    child[2].insert(mean,gmean,indice);
                else
                    child[3].insert(mean,gmean,indice);
            }
        }
        else {                                   //Node 4-7
            if(mean[1]<center[1]) {              //Node 4-5
                if(mean[2]<center[2])
                    child[4].insert(mean,gmean,indice);
                else
                    child[5].insert(mean,gmean,indice);
            }
            else {                               //Node 6-7
                if(mean[2]<center[2])
                    child[6].insert(mean,gmean,indice);
                else
                    child[7].insert(mean,gmean,indice);
            }
        }
        /* }}} */
    }
    /* }}} */
    /* Update mean position. {{{ */
    for(int i=0;i<3;i++) {
        mean[i]+=pos[i];
        gmean[i]+=grad[i];
        gcmean[i]+=grad[i]*grad[i];
    }
    gcmean[3]+=grad[0]*grad[1];
    gcmean[4]+=grad[0]*grad[2];
    gcmean[5]+=grad[1]*grad[2];
    n_q++;
    indice=-1;
    /* }}} */
    /* Insert new charge. {{{ */
    if(pos[0]<center[0]) {                      //Node 0-3
        if(pos[1]<center[1]) {                  //Node 0-1
            if(pos[2]<center[2])
                child[0].insert(pos,grad,ind);
            else
                child[1].insert(pos,grad,ind);
        }
        else {                                   //Node 2-3
            if(pos[2]<center[2])
                child[2].insert(pos,grad,ind);
            else
                child[3].insert(pos,grad,ind);
        }
    }
    else {                                       //Node 4-7
        if(pos[1]<center[1]) {                  //Node 4-5
            if(pos[2]<center[2])
                child[4].insert(pos,grad,ind);
            else
                child[5].insert(pos,grad,ind);
        }
        else {                                   //Node 6-7
            if(pos[2]<center[2])
                child[6].insert(pos,grad,ind);
            else
                child[7].insert(pos,grad,ind);
        }
    }
    /* }}} */
    return;
}
/* }}} */
/* energy_grad: optimum (flat). {{{ */
/*!
 * More or less basic optimization performed:
 *  - prefetching tested but no improvement (why ?!?),
 *  - no left hand side assignement to array elements in the loop.
 *  TODO: read assembler output and try to find more optimizations...
 */
int octTree::energy_grad(double *grad, double &epot) {
    if(n_q<2) {                                  //Trivial case, nothing to do
        return 0;
    }
#ifdef DEBUG
    int res=0;
#endif
    int nint=0;
    epot*=2;                                     //Avoid divisions inside loop
    octTree *cur=this;
    do {                                         //Walk the tree
        if(cur->n_q==1) {
            int i=cur->indice;
            double gx,gy,gz;
            gx=gy=gz=0.0;
            octTree *aux=this;
            do {
                int nq=aux->n_q;
                double rx,ry,rz;
                rx=aux->mean[0]-cur->mean[0];
                double invr2=rx*rx;
                ry=aux->mean[1]-cur->mean[1];
                invr2+=ry*ry;
                rz=aux->mean[2]-cur->mean[2];
                invr2+=rz*rz;
                //Test cell
                if ((nq==1 && i!=aux->indice) || invr2>aux->delta) {
                    aux=aux->skip;
                    nint++;
                    invr2=1.0/invr2;
                    double invr=sqrt(invr2)*nq;
                    invr2*=invr;
#ifdef DEBUG
                    res+=nq;
#endif
                    epot+=invr;
                    rx*=invr2; gx+=rx;
                    ry*=invr2; gy+=ry;
                    rz*=invr2; gz+=rz;
                } else {                         //Go deeper
                    aux=aux->next;
                }
            } while(aux!=0);
            int ii=i*3;
            grad[ii]+=gx;
            grad[ii+1]+=gy;
            grad[ii+2]+=gz;
        }
        cur=cur->next;
    } while(cur!=0);
#ifdef DEBUG
    if(res!=n_q*(n_q-1)) {
        cerr << "[E] Error!" << endl;
        return -1;
    }
#endif
    epot/=2;
    return nint;
}
/* }}} */
/* guess: (flat). {{{ */
int octTree::guess(double &fp, double &fpp) {
    if(n_q<2) {                       //Trivial case, nothing to do
        return 0;
    }
    int res=0;
    //Avoid lots of divisions inside the loop
    fp*=2;
    fpp*=2;
    octTree *cur=this;
    do {
        if(cur->n_q==1) {
            int i=cur->indice;
            octTree *aux=this;
            do {
                if(aux->n_q==1 && aux->indice!=i) {
                    double rx=cur->mean[0]-aux->mean[0];
                    double invr2=rx*rx;
                    double ry=cur->mean[1]-aux->mean[1];
                    invr2+=ry*ry;
                    double rz=cur->mean[2]-aux->mean[2];
                    invr2+=rz*rz;
                    invr2=1.0/invr2;
                    double invr=sqrt(invr2);
                    double invr3=invr*invr2;
                    double gx=aux->gmean[0]-cur->gmean[0];
                    double rg=rx*gx;
                    double gg=gx*gx;
                    double gy=aux->gmean[1]-cur->gmean[1];
                    rg+=ry*gy;
                    gg+=gy*gy;
                    double gz=aux->gmean[2]-cur->gmean[2];
                    rg+=rz*gz;
                    gg+=gz*gz;
                    fp+=invr3*rg;
                    fpp+=invr3*((rg*rg)*(3*invr2)-gg);
                    res++;
                    aux=aux->skip;
                } else {
                    double rx=cur->mean[0]-aux->mean[0];
                    double invr2=rx*rx;
                    double ry=cur->mean[1]-aux->mean[1];
                    invr2+=ry*ry;
                    double rz=cur->mean[2]-aux->mean[2];
                    invr2+=rz*rz;
                    if (invr2>aux->delta) { //Test Cell
                        double gg=aux->g2mean;
                        double gx=aux->gmean[0]-cur->gmean[0];
                        double rg=rx*gx;
                        gg-=cur->gmean[0]*(gx+aux->gmean[0]);
                        double gy=aux->gmean[1]-cur->gmean[1];
                        rg+=ry*gy;
                        gg-=cur->gmean[1]*(gy+aux->gmean[1]);
                        double gz=aux->gmean[2]-cur->gmean[2];
                        rg+=rz*gz;
                        gg-=cur->gmean[2]*(gz+aux->gmean[2]);
                        double rg2=rx*cur->gmean[0]+ry*cur->gmean[1]
                            +rz*cur->gmean[2];
                        rg2=rg2*(rg2-2*(rx*aux->gmean[0]+ry*aux->gmean[1]
                                    +rz*aux->gmean[2]));
                        rg2+=rx*((rx*aux->gcmean[0])+(ry*aux->gcmean[3]));
                        rg2+=ry*((ry*aux->gcmean[1])+(rz*aux->gcmean[5]));
                        rg2+=rz*((rz*aux->gcmean[2])+(rx*aux->gcmean[4]));
                        invr2=1.0/invr2;
                        double invr=sqrt(invr2);
                        double invr3=invr*invr2*aux->n_q;
                        fp+=invr3*rg;
                        fpp+=invr3*(3*invr2*rg2-gg);
                        res++;
                        aux=aux->skip;
                    } else {                         //Go deeper
                        aux=aux->next;
                    }
                }
            } while(aux!=0);
        }
        cur=cur->next;
    } while(cur!=0);
    fp/=2;
    fpp/=2;
    return res;
}
/* }}} */
#endif //VECTORIZE
/* octtree.cpp */
