#pragma once

/*
// DEFINE ALL THE PARAMETERS HERE 
//#define iter		1
//#define VER			32		 //  width of lattice matrix. for ex 6*6 lattice
//#define DIM         2		// dimension of lattice // here 2; square lattice
//#define latt_pc     1		 // percentage of lattice points where atom exist
//#define upspin_pc   0.5		 // percentage of upspin in lattice
//#define w           1        // omega(w) for bimmodal distribution of Bi in eq 6.1    
//#define N			VER * VER   // width of edge matrix. for ex 36*36 edges
//#define V			N+2       // source + width of edge matrix + sink
//#define J			1		
//#define INFINITE 100000
//#define HIGH INFINITE

//delta is multipling factor with Bmat
//#define del 0.7
//#define del_beg     .7
//#define del_end     1.5
//#define del_inc     .1

 
// text substitute
//#define tab			"\t"
//#define d			" : "

// used in push relabel
//#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
*/

#define fixed_float(x) std::fixed <<std::setprecision(3)<<(x)

#define VER 40
#define N VER*VER
#define tab	"\t"
#define d ":"
#define J 1
#define INFINITE 100000
#define V N+2

#define ITERS 1
#define delta_beg 1.8
#define delta_end 1.81
#define delta_inc 0.1
#define si_beg 1.6
#define si_end 2.01
#define si_inc 0.1

#define NTEMP 15
#define NMEAS 1000
#define NS VER*VER
