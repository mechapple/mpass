/*
 * num_int_2D_test.cxx
 * 
 * Copyright 2017 Anirban <anirban@ZeroPointEnergy>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#include <iostream>
#include <cmath>
#include <ctime>
#include "cubature.h"


double func(double x1,double x2) { 
	double y = sqrt(1 + x1*x1 + x2*x2);
	return 1.0/y; 
}

int f1(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
    // compute the output value: note that fdim should == 1 from below
    fval[0] = 1.0/sqrt(1 + x[0]*x[0] + x[1]*x[1]);
    return 0; // success
}

int main(int argc, char **argv)
{
	int start_s=clock();
	int Nsim = 1000;
	double h = 1/(2.0*Nsim), sum = 0.0, lsum[3], mid1, mid2;
	double lam1 = 3.0*(35.0+sqrt(385.0))/140.0, lam2 = 3.0*(35.0-sqrt(385.0))/140.0;
	double b1 = (77.0-3.0*sqrt(385.0))/891.0, b2 = (77.0+3.0*sqrt(385.0))/891.0;
	double mu = sqrt(3.0/5.0), c = 25.0/324.0;
	double m[2][2] = {{4,8},{8,16}};
	
	for(int i=0;i<Nsim;i++) 
		for(int j=0;j<Nsim;j++) {
			mid1 = i*2*h + h; mid2 = j*2*h + h;
			
			lsum[0] = func(mid1 + lam1*h,mid2) + func(mid1 - lam1*h,mid2) + func(mid1,mid2 + lam1*h) + func(mid1,mid2 - lam1*h);
			lsum[1] = func(mid1 + lam2*h,mid2) + func(mid1 - lam2*h,mid2) + func(mid1,mid2 + lam2*h) + func(mid1,mid2 - lam2*h);
			lsum[2] = func(mid1 - mu*h,mid2 - mu*h) + func(mid1 - mu*h,mid2 + mu*h) + func(mid1 + mu*h,mid2 - mu*h) + func(mid1 + mu*h,mid2 + mu*h);
			sum += b1*lsum[0] + b2*lsum[1] + c*lsum[2];
	}
	
	sum *= 4.0*h*h;
	double error = fabs(sum - 0.793359121);
	//std::cout << "Integral = " << sum*h*h << "\n";
	printf("Integral %.6f %.15f %.15g\n",h,sum,error);
	
	int stop_s=clock();
	std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;
	
	double xmin[2] = {0,0}, xmax[2] = {1,1}, sigma = 0.5, val[1], err[1];
    hcubature(1, f1, &sigma, 2, xmin, xmax, 0, 0, 1e-8, ERROR_INDIVIDUAL, val, err);
    printf("Computed integral = %0.15f +/- %.15g\n", val[0], err[0]);
	
	int stop_s1=clock();
	std::cout << "time: " << (stop_s1-stop_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

	return 0;
}

