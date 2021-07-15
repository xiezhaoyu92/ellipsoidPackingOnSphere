//
//  particles.c
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#include "particles.h"
#include <math.h>
#include "operation.h"
#include "surface.h"


//define contact function of two ellispoids, calculate its derivative at lambda ranging from 0 to 1. This is a setup for finding lambda where the contact function derivative is 0.
//This is not used because it's not accurate, sometimes judging two particles overlap even if they are not.
/*double contactFunctionDerivative(particle *p1, particle *p2, double lambda) {
    double *u1 = p1->director;
    double *x1 = p1->position;
    double *u2 = p2->director;
    double *x2 = p2->position;
    double a1 = p1->aa;
    double b1 = p1->bb;
    double a2 = p2->aa;
    double b2 = p2->bb;
    double sx, sy, sz, sxy, syz, sxz, dsx, dsy, dsz, dsxy, dsyz,dsxz;
    double numerator, denominator, dnumerator, ddenominator, dpotential;
    
    sx = (1-lambda)*(b1*b1+(a1*a1-b1*b1)*u1[0]*u1[0])+lambda*(b2*b2+(a2*a2-b2*b2)*u2[0]*u2[0]);
    sy = (1-lambda)*(b1*b1+(a1*a1-b1*b1)*u1[1]*u1[1])+lambda*(b2*b2+(a2*a2-b2*b2)*u2[1]*u2[1]);
    sz = (1-lambda)*(b1*b1+(a1*a1-b1*b1)*u1[2]*u1[2])+lambda*(b2*b2+(a2*a2-b2*b2)*u2[2]*u2[2]);
    sxy = (1-lambda)*(a1*a1-b1*b1)*u1[0]*u1[1]+lambda*(a2*a2-b2*b2)*u2[0]*u2[1];
    syz = (1-lambda)*(a1*a1-b1*b1)*u1[1]*u1[2]+lambda*(a2*a2-b2*b2)*u2[1]*u2[2];
    sxz = (1-lambda)*(a1*a1-b1*b1)*u1[0]*u1[2]+lambda*(a2*a2-b2*b2)*u2[0]*u2[2];
    dsx = (a2*a2-b2*b2)*u2[0]*u2[0]-(a1*a1-b1*b1)*u1[0]*u1[0];
    dsy = (a2*a2-b2*b2)*u2[1]*u2[1]-(a1*a1-b1*b1)*u1[1]*u1[1];
    dsz = (a2*a2-b2*b2)*u2[2]*u2[2]-(a1*a1-b1*b1)*u1[2]*u1[2];
    dsxy = (a2*a2-b2*b2)*u2[0]*u2[1]-(a1*a1-b1*b1)*u1[0]*u1[1];
    dsyz = (a2*a2-b2*b2)*u2[1]*u2[2]-(a1*a1-b1*b1)*u1[1]*u1[2];
    dsxz = (a2*a2-b2*b2)*u2[0]*u2[2]-(a1*a1-b1*b1)*u1[0]*u1[2];
    
    numerator = (x1[0]-x2[0])*(x1[0]-x2[0])*(sy*sz-syz*syz)+(x1[1]-x2[1])*(x1[1]-x2[1])*(sx*sz-sxz*sxz)+(x1[2]-x2[2])*(x1[2]-x2[2])*(sx*sy-sxy*sxy)+2*(x1[0]-x2[0])*(x1[1]-x2[1])*(sxz*syz-sz*sxy)+2*(x1[1]-x2[1])*(x1[2]-x2[2])*(sxy*sxz-sx*syz)+2*(x1[0]-x2[0])*(x1[2]-x2[2])*(sxy*syz-sy*sxz);
    dnumerator = (x1[0]-x2[0])*(x1[0]-x2[0])*(sy*dsz+dsy*sz-2*syz*dsyz)+(x1[1]-x2[1])*(x1[1]-x2[1])*(sx*dsz+dsx*sz-2*sxz*dsxz)+(x1[2]-x2[2])*(x1[2]-x2[2])*(sx*dsy+dsx*sy-2*sxy*dsxy)+2*(x1[0]-x2[0])*(x1[1]-x2[1])*(sxz*dsyz+dsxz*syz-sz*dsxy-dsxy*sz)+2*(x1[1]-x2[1])*(x1[2]-x2[2])*(sxy*dsxz+dsxy*sxz-sx*dsyz-dsx*syz)+2*(x1[0]-x2[0])*(x1[2]-x2[2])*(sxy*dsyz+dsxy*syz-sy*dsxz-dsy*sxz);
    denominator = 2*sxy*syz*sxz+sx*sy*sz-sx*syz*syz-sy*sxz*sxz-sz*sxy*sxy;
    ddenominator = 2*(sxy*syz*dsxz+sxy*dsyz*sxz+dsxy*syz*sxz)+(sx*sy*dsz+sx*dsy*sz+dsx*sy*sz)-dsx*syz*syz-2*sx*syz*dsyz-dsy*sxz*sxz-2*sy*sxz*dsxz-dsz*sxy*sxy-2*sz*sxy*dsxy;
    
    dpotential=(((1-2*lambda)*numerator+lambda*(1-lambda)*dnumerator)*denominator-lambda*(1-lambda)*numerator*ddenominator)/(denominator*denominator);
    return dpotential;
}*/

//calculate contact function at lambda ranging from 0 to 1.
double contactFunction(particle *p1, particle *p2, double lambda) {
    double *u1 = p1->director;
    double *x1 = p1->position;
    double *u2 = p2->director;
    double *x2 = p2->position;
    double a1 = p1->aa;
    double b1 = p1->bb;
    double a2 = p2->aa;
    double b2 = p2->bb;
    double sx, sy, sz, sxy, syz, sxz, numerator, denominator, potential;
    
    sx = (1-lambda)*(b1*b1+(a1*a1-b1*b1)*u1[0]*u1[0])+lambda*(b2*b2+(a2*a2-b2*b2)*u2[0]*u2[0]);
    sy = (1-lambda)*(b1*b1+(a1*a1-b1*b1)*u1[1]*u1[1])+lambda*(b2*b2+(a2*a2-b2*b2)*u2[1]*u2[1]);
    sz = (1-lambda)*(b1*b1+(a1*a1-b1*b1)*u1[2]*u1[2])+lambda*(b2*b2+(a2*a2-b2*b2)*u2[2]*u2[2]);
    sxy = (1-lambda)*(a1*a1-b1*b1)*u1[0]*u1[1]+lambda*(a2*a2-b2*b2)*u2[0]*u2[1];
    syz = (1-lambda)*(a1*a1-b1*b1)*u1[1]*u1[2]+lambda*(a2*a2-b2*b2)*u2[1]*u2[2];
    sxz = (1-lambda)*(a1*a1-b1*b1)*u1[0]*u1[2]+lambda*(a2*a2-b2*b2)*u2[0]*u2[2];
    
    numerator = (x1[0]-x2[0])*(x1[0]-x2[0])*(sy*sz-syz*syz)+(x1[1]-x2[1])*(x1[1]-x2[1])*(sx*sz-sxz*sxz)+(x1[2]-x2[2])*(x1[2]-x2[2])*(sx*sy-sxy*sxy)+2*(x1[0]-x2[0])*(x1[1]-x2[1])*(sxz*syz-sz*sxy)+2*(x1[1]-x2[1])*(x1[2]-x2[2])*(sxy*sxz-sx*syz)+2*(x1[0]-x2[0])*(x1[2]-x2[2])*(sxy*syz-sy*sxz);
    denominator = 2*sxy*syz*sxz+sx*sy*sz-sx*syz*syz-sy*sxz*sxz-sz*sxy*sxy;
    
    potential=lambda*(1-lambda)*numerator/denominator;
    return potential;
}
//calculate the derivative of contact funtion by nemeric method
double contactFunctionDerivative(particle *p1, particle *p2, double lambda) {
    double delta = 0.000001;
    return (contactFunction(p1, p2, lambda+delta)-contactFunction(p1, p2, lambda-delta))/2/delta;
}

//Check if two ellipsoids overlap, from the max of contactFunction, if max >= 1, the two ellipsoids don't overlap. Using Brent method to find lambda where the derivative is 0, and calculate the related contact function as max.
/*int overlapQ(particle *p1, particle *p2) {
    double c, d, e, fa, fb, fc, m, p, q, r, s, sa, sb, tol;
    double machep = MACHINE_EPSILON;
    double t = MACHINE_EPSILON;
    sa = 0;
    sb = 1;
    fa = contactFunctionDerivative(p1, p2, sa);
    fb = contactFunctionDerivative(p1, p2, sb);
    c = sa;
    fc = fa;
    e = sb - sa;
    d = e;
    for ( ; ; ) {
        if (fabs(fc)<fabs(fb)) {
            sa = sb;
            sb = c;
            c = sa;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        tol = 2.0*machep*fabs(sb)+t;
        m = 0.5*(c-sb);
        if (fabs(m)<=tol||fabs(fb)<=MACHINE_EPSILON)
            break;
        if (fabs(e)<tol||fabs(fa)<=fabs(fb)) {
            e = m;
            d = e;
        }
        else {
            s = fb/fa;
            if (sa==c){
                p = 2.0*m*s;
                q = 1.0-s;
            }
            else {
                q = fa/fc;
                r = fb/fc;
                p = s*(2.0*m*q*(q-r)-(sb-sa)*(r-1.0));
                q = (q-1.0)*(r-1.0)*(s-1.0);
            }
            if (0.0<p)
                q = -q;
            else
                p = -p;
            s = e;
            e = d;
            if (2.0*p<3.0*m*q-fabs(tol*q) && p<fabs(0.5*s*q))
                d = p/q;
            else {
                e = m;
                d = e;
            }
        }
        sa = sb;
        fa = fb;
        if (tol<fabs(d))
            sb = sb+d;
        else if (0.0<m)
            sb = sb+tol;
        else
            sb = sb-tol;
        fb = contactFunctionDerivative(p1, p2, sb);
        if ((0.0<fb&&0.0<fc)||(fb<=0.0&&fc<=0.0)) {
            c = sa;
            fc = fa;
            e = sb-sa;
            d = e;
        }
    }
    
    double max;
    max = contactFunction(p1, p2, sb);
    if (max>=1)
        return FALSE;
    else
        return TRUE;
}*/
//compare contact distance and center-to-center distance to judge overlap, faster than previous method.
int overlapQ(particle *p1, particle *p2) {
    double *position1 = p1->position;
    double *position2 = p2->position;
    double *u1 = p1->director;
    double *u2 = p2->director;
    double a = p1->aa;
    double b = p1->bb;
    
    double r[3];
    
    vectorSubtract(position2, position1, r);
    double R = vectorNorm(r);
    vectorNormalize(r, r);
    
    double chi = (a*a-b*b)/(a*a+b*b);
    double n1 = vectorDotProduct(r, u1);
    double n2 = vectorDotProduct(r, u2);
    double nu = vectorDotProduct(u1, u2);
    double contactDistance = 2*b/sqrt(1-chi/2*((n1+n2)*(n1+n2)/(1+chi*nu)+(n1-n2)*(n1-n2)/(1-chi*nu)));
    
    if(R < contactDistance)
        return TRUE;
    else
        return FALSE;
}

/*//for a specific particle i in all particles p, check if it overlaps with its surrounding particles.
int anyOverlapQ(int i, particle p[], int np) {
    int overlaps = FALSE;
    for (int j=0; j<np; j++) {
        if(j==i)
            continue;
        if(overlapQ(&(p[i]), &(p[j]))) {
            overlaps = TRUE;
            break;
        }
    }
    if(overlaps)
        return TRUE;
    else
        return FALSE;
}*/

//for a specific particle i in all particles p, check if it overlaps with its surrounding particles. We define surrounding particles as those in the neighboring partitions
int anyOverlapQ(int i, particle p[], int np) {
    int overlaps = FALSE;
    for(int j=0; j<np; j++) {
        if(j==i)
            continue;
        if((p[j].coord[0]-p[i].coord[0]==-1||p[j].coord[0]-p[i].coord[0]==0||p[j].coord[0]-p[i].coord[0]==1)&&(p[j].coord[1]-p[i].coord[1]==-1||p[j].coord[1]-p[i].coord[1]==0||p[j].coord[1]-p[i].coord[1]==1)&&(p[j].coord[2]-p[i].coord[2]==-1||p[j].coord[2]-p[i].coord[2]==0||p[j].coord[2]-p[i].coord[2]==1)) {
            if(overlapQ(&(p[i]), &(p[j]))) {
                overlaps = TRUE;
                break;
            }
        }
    }
    if(overlaps)
        return TRUE;
    else
        return FALSE;
}
//undoOverlaps
//define the interaction of two overlapped ellipsoids, based on the Gaussian model potential of Berne and Pechukas.
double gaussianModelPotential(double x1, double y1, double z1, double theta1, double x2, double y2, double z2, double theta2, double a, double b, double Rad) {
    double epsilon0 = 10; //Control the magnitude of potential
    double position1[4] = {x1, y1, z1, theta1};
    double position2[4] = {x2, y2, z2, theta2};
    double u1[3];
    double u2[3];
    directorCalc(Rad, position1, u1);
    directorCalc(Rad, position2, u2);
    
    double r[3];
    
    vectorSubtract(position2, position1, r);
    double R = vectorNorm(r);
    vectorNormalize(r, r);
    
    double chi = (a*a-b*b)/(a*a+b*b);
    double n1 = vectorDotProduct(r, u1);
    double n2 = vectorDotProduct(r, u2);
    double nu = vectorDotProduct(u1, u2);
    double sigma = 2*b*b/(1-chi/2*((n1+n2)*(n1+n2)/(1+chi*nu)+(n1-n2)*(n1-n2)/(1-chi*nu)));
    double epsilon = epsilon0/sqrt(1-chi*chi*nu*nu);
    return epsilon*exp(-R*R/sigma);
}
//use the gradient of potential for different coordinates to get the force
void addOverlapForce(particle *p1, particle *p2, double a, double b, double Rad) {
    double x1 = p1->position[0];
    double y1 = p1->position[1];
    double z1 = p1->position[2];
    double theta1 = p1->position[3];
    double x2 = p2->position[0];
    double y2 = p2->position[1];
    double z2 = p2->position[2];
    double theta2 = p2->position[3];
    
    double delta = 0.0001;
    p1->force[0] = p1->force[0]+(gaussianModelPotential(x1-delta,y1,z1,theta1,x2,y2,z2,theta2,a,b,Rad)-gaussianModelPotential(x1+delta,y1,z1,theta1,x2,y2,z2,theta2,a,b,Rad))/2/delta;
    p1->force[1] = p1->force[1]+(gaussianModelPotential(x1,y1-delta,z1,theta1,x2,y2,z2,theta2,a,b,Rad)-gaussianModelPotential(x1,y1+delta,z1,theta1,x2,y2,z2,theta2,a,b,Rad))/2/delta;
    p1->force[2] = p1->force[2]+(gaussianModelPotential(x1,y1,z1-delta,theta1,x2,y2,z2,theta2,a,b,Rad)-gaussianModelPotential(x1,y1,z1+delta,theta1,x2,y2,z2,theta2,a,b,Rad))/2/delta;
    p1->force[3] = p1->force[3]+10*(gaussianModelPotential(x1,y1,z1,theta1-delta,x2,y2,z2,theta2,a,b,Rad)-gaussianModelPotential(x1,y1,z1,theta1+delta,x2,y2,z2,theta2,a,b,Rad))/2/delta;
    p2->force[0] = p2->force[0]+(gaussianModelPotential(x1,y1,z1,theta1,x2-delta,y2,z2,theta2,a,b,Rad)-gaussianModelPotential(x1,y1,z1,theta1,x2+delta,y2,z2,theta2,a,b,Rad))/2/delta;
    p2->force[1] = p2->force[1]+(gaussianModelPotential(x1,y1,z1,theta1,x2,y2-delta,z2,theta2,a,b,Rad)-gaussianModelPotential(x1,y1,z1,theta1,x2,y2+delta,z2,theta2,a,b,Rad))/2/delta;
    p2->force[2] = p2->force[2]+(gaussianModelPotential(x1,y1,z1,theta1,x2,y2,z2-delta,theta2,a,b,Rad)-gaussianModelPotential(x1,y1,z1,theta1,x2,y2,z2+delta,theta2,a,b,Rad))/2/delta;
    p2->force[3] = p2->force[3]+10*(gaussianModelPotential(x1,y1,z1,theta1,x2,y2,z2,theta2-delta,a,b,Rad)-gaussianModelPotential(x1,y1,z1,theta1,x2,y2,z2,theta2+delta,a,b,Rad))/2/delta;
}

//project particle onto current sphere surface, and update the director u
void projectOntoSphere(particle *p, double Rad) {
    double r = vectorNorm(p->position);
    for(int i=0; i<3; i++)
        p->position[i] = p->position[i]/r*Rad;
    directorCalc(Rad, p->position, p->director);
}

//move particle, and project onto the surface. Update the information of particle
void gradDescentConstrainedIntegrate(particle *p, double dt, double Rad) {
    //move particle
    for(int i=0; i<4; i++)
        p->position[i] = p->position[i] + (p->force[i])*dt + (p->forceStoc[i])*sqrt(dt);
    if(p->position[3]>=PI)
        p->position[3] = p->position[3]-PI;
    else if(p->position[3]<0)
        p->position[3] = p->position[3]+PI;
    
    //project onto sphere
    projectOntoSphere(p, Rad);
}

/*//undooverlaps
void undoOverlaps(particle p[], int maxUndoSteps, double dtOverlap, int *nearJammed, double a, double b, double Rad, int np) {
    for(int i=0; i<np; i++)
        for(int j=0; j<4; j++)
            p[i].forceStoc[j] = 0;
    
    int undoSteps = 0;
    *nearJammed = FALSE;
    int totalOverlapQ;
    do {
        totalOverlapQ = FALSE;
        undoSteps++;
        if(undoSteps>maxUndoSteps) {
            *nearJammed = TRUE;
            break;
        }
        for(int i=0; i<np; i++)
            for(int j=0; j<4; j++)
                p[i].force[j] = 0;
        for(int i=0; i<np; i++)
            for(int j=i+1; j<np; j++)
                if(overlapQ(&(p[i]), &(p[j]))) {
                    addOverlapForce(&(p[i]), &(p[j]), a, b, Rad);
                    totalOverlapQ = TRUE;
                }
        for(int i=0; i<np; i++)
            gradDescentConstrainedIntegrate(&(p[i]), dtOverlap, Rad);
    
    } while(totalOverlapQ);
    //printf("undoOverlapSteps are %d\n", undoSteps);
}*/
//undooverlaps by checking particle with those in surrounding partitions
void undoOverlaps(particle p[], int maxUndoSteps, double dtOverlap, int *nearJammed, double a, double b, double Rad, int np, double partBoundX[], double partBoundY[], double partBoundZ[], int nPart[]) {
    for(int i=0; i<np; i++)
        for(int j=0; j<4; j++)
            p[i].forceStoc[j] = 0;
    
    int undoSteps = 0;
    *nearJammed = FALSE;
    int totalOverlapQ;
    do {
        totalOverlapQ = FALSE;
        undoSteps++;
        if(undoSteps>maxUndoSteps) {
            *nearJammed = TRUE;
            break;
        }
        for(int i=0; i<np; i++)
            for(int j=0; j<4; j++)
                p[i].force[j] = 0;
        for(int i=0; i<np; i++)
            for(int j=i+1; j<np; j++)
                if((p[j].coord[0]-p[i].coord[0]==-1||p[j].coord[0]-p[i].coord[0]==0||p[j].coord[0]-p[i].coord[0]==1)&&(p[j].coord[1]-p[i].coord[1]==-1||p[j].coord[1]-p[i].coord[1]==0||p[j].coord[1]-p[i].coord[1]==1)&&(p[j].coord[2]-p[i].coord[2]==-1||p[j].coord[2]-p[i].coord[2]==0||p[j].coord[2]-p[i].coord[2]==1)&&overlapQ(&(p[i]), &(p[j]))) {
                    addOverlapForce(&(p[i]), &(p[j]), a, b, Rad);
                    totalOverlapQ = TRUE;
                }
        for(int i=0; i<np; i++) {
            vectorCopy(p[i].position, p[i].oldPosition, 4);
            vectorCopy(p[i].director, p[i].oldDirector, 3);
            gradDescentConstrainedIntegrate(&(p[i]), dtOverlap, Rad);
            if(isnan(p[i].position[0])){
                vectorCopy(p[i].oldPosition, p[i].position, 4);
                vectorCopy(p[i].oldDirector, p[i].director, 3);
            }
        }
        for(int i=0; i<np; i++)
            partitionOneParticle(&(p[i]), partBoundX, partBoundY, partBoundZ, nPart);
        
    } while(totalOverlapQ);
    //printf("undoOverlapSteps are %d\n", undoSteps);
}

//Add stochastic force
void addStochasticForce(particle *p, double Da, double Db, double Dth, double Rad) {
    double theta = p->position[3];
    double t1[3];
    double t2[3];
    tangent1(Rad, p->position, t1);
    tangent2(Rad, p->position, t2);
    double ETAa, ETAb, ETAth;
    ETAa = randNormal();
    ETAb = randNormal();
    ETAth = randNormal();
    
    for(int i=0; i<3; i++)
        p->forceStoc[i] = (ETAa*cos(theta)*sqrt(2*Da)-ETAb*sin(theta)*sqrt(2*Db))*t1[i] + (ETAa*sin(theta)*sqrt(2*Da)+ETAb*cos(theta)*sqrt(2*Db))*t2[i];
    p->forceStoc[3] = ETAth*sqrt(2*Dth);
}

//Put one particle into its related partition, setting the value of coord[3]
void partitionOneParticle(particle *p, double partBoundX[], double partBoundY[], double partBoundZ[], int nPart[]) {
    int i=0, j=0, k=0;
    while(i<nPart[0]&&p->position[0]>=partBoundX[i])
        i++;
    p->coord[0] = i-1;
    
    while(j<nPart[1]&&p->position[1]>=partBoundY[j])
        j++;
    p->coord[1] = j-1;
    
    while(k<nPart[2]&&p->position[2]>=partBoundZ[k])
        k++;
    p->coord[2] = k-1;
}

