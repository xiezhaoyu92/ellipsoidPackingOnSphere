//
//  surface.c
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#include <math.h>
#include "surface.h"

double conformalInvert(const double u2) {
    return 2*u2-1;
}

void tangent1(double Rad, double *position, double *t1) {
    double u = position[2]/Rad;
    double rho = sqrt(position[0]*position[0]+position[1]*position[1]);
    if (fabs(rho) < MACHINE_EPSILON) {
        t1[0] = 1;
        t1[1] = 0;
        t1[2] = 0;
    }
    else {
        t1[0] = -position[1]/rho;
        t1[1] = position[0]/rho;
        t1[2] = 0;
    }
}

void tangent2(double Rad, double *position, double *t2) {
    double u = position[2]/Rad;
    double rho = sqrt(position[0]*position[0]+position[1]*position[1]);
    if (fabs(rho) < MACHINE_EPSILON){
        t2[0] = 0;
        t2[1] = 1;
        t2[2] = 0;
    }
    else {
        t2[0] = -u*position[0]/rho;
        t2[1] = -u*position[1]/rho;
        t2[2] = sqrt(1 - u*u);
    }

}
