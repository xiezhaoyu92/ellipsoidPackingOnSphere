//
//  surface.h
//  ellipsoidalPacking
//
//  Created by Zhaoyu Xie on 10/26/17.
//  Copyright © 2017 Zhaoyu Xie. All rights reserved.
//

#ifndef surface_h
#define surface_h

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON 1e-15
#endif

#include <stdio.h>

double conformalInvert(const double u2);
void tangent1(double Rad, double *position, double *t1);
void tangent2(double Rad, double *position, double *t2);
#endif /* surface_h */
