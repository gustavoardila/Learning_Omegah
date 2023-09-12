#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"
#include <math.h>
#include <cmath>
#include <stdlib.h> 
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <bits/stdc++.h> 
using namespace std;



//--------------------------------------------------------------------------------------------------------
//------------This part of the code allows the product between a complex and an int-----------------------
//--------------------------------------------------------------------------------------------------------

// Trick to allow type promotion below
template <typename T>
struct identity_t { typedef T type; };

/// Make working with std::complex<> nubmers suck less... allow promotion.
#define COMPLEX_OPS(OP)                                                 \
  template <typename _Tp>                                               \
  std::complex<_Tp>                                                     \
  operator OP(std::complex<_Tp> lhs, const typename identity_t<_Tp>::type & rhs) \
  {                                                                     \
    return lhs OP rhs;                                                  \
  }                                                                     \
  template <typename _Tp>                                               \
  std::complex<_Tp>                                                     \
  operator OP(const typename identity_t<_Tp>::type & lhs, const std::complex<_Tp> & rhs) \
  {                                                                     \
    return lhs OP rhs;                                                  \
  }
COMPLEX_OPS(+)
COMPLEX_OPS(-)
COMPLEX_OPS(*)
COMPLEX_OPS(/)
#undef COMPLEX_OPS
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------




//-----------------------Main----------------------------------



int main(int argc,char** argv)
{
   srand(time(NULL)); 	//initialise la sequence de la fonction rand afin que la serie ne commence pas tjs au meme pt (va prendre l'heure)
   int err;
   char cdmName[10];
   int spin2, charge3, cdim;

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */
  VZdecay=0; VWdecay=0;  
  
   int taille = 100;
   double Q;
   char chaine[taille] ;
   int errSph;

	int nblines = 324;
	ifstream inFile;

//--------------------------------------------------------
//                Parameter definitions
//--------------------------------------------------------

	//COUPLINGS FROM THE SCALAR SECTOR
	
  double l1 ; //Higgs coupling, it has to be changed to obtain Mh = 125Gev
  double lSEt;    //Singlet Scalar coupling with H
  double l4Si,l4Eta;  
  double l31EtH,l32EtH,l33EtH ;    
  double lHSi, lEtSi ;   


  //coupling of : Yn L Et N
  double Yn11, Yn11c ; //Real and Imaginary parts
  double Yn12, Yn12c ; 
  double Yn13, Yn13c ;  
  
  double Yn21, Yn21c ; 
  double Yn22, Yn22c ; 
  double Yn23, Yn23c ;  

  double Yn31, Yn31c ; 
  double Yn32, Yn32c ; 
  double Yn33, Yn33c ;  


  //Kappa couplings
  double Kap11, Kap11c ; 
  double Kap12, Kap12c ; 
  double Kap13, Kap13c ;  

  double Kap21, Kap21c ; 
  double Kap22, Kap22c ; 
  double Kap23, Kap23c ;

  double Kap31, Kap31c ; 
  double Kap32, Kap32c ; 
  double Kap33, Kap33c ;
  

  //Masses: tree level definitons and parameters from the Lagrangien
  double MSi2, MEt2;
  double mh01_tree, mh02_tree, met0R_tree, met0I_tree, metc_tree;

  //Physical Masses definition: Will be the masses computed at 1-loop by SPheno
  double MEtac,MetaR,MetaI;//charged scalar
  double MX01,MX02,MX03;

  double v = 246.22; //vev [GeV]
  double Mh1, Mh2;
  double vSi;   

  //Neutrino masses
  double mnu1 ;   
  double mnu2 ;
  double mnu3 ;
  
	//MATRICES OF THE MODEL
  
  //Element of the Scalar Mixing Matrix
  double ZH_11,ZH_12,ZH_21,ZH_22 ;   
  double ImZH_11,ImZH_12,ImZH_21,ImZH_22 ; 
   
  double ZA_11,ZA_12,ZA_21,ZA_22 ;   
  double ImZA_11,ImZA_12,ImZA_21,ImZA_22 ; 

  //Element of the Fermion Mixing Matrix
  double ZX_11,ZX_12,ZX_13,ZX_21,ZX_22,ZX_23,ZX_31,ZX_32,ZX_33 ;   
  double ImZX_11,ImZX_12,ImZX_13,ImZX_21,ImZX_22,ImZX_23 ,ImZX_31,ImZX_32,ImZX_33;

  double ZDL11,ZDL12,ZDL13,ZDL21,ZDL22,ZDL23,ZDL31,ZDL32,ZDL33, IZDL11,IZDL12,IZDL13,IZDL21,IZDL22,IZDL23,IZDL31,IZDL32,IZDL33;
  double ZDR11,ZDR12,ZDR13,ZDR21,ZDR22,ZDR23,ZDR31,ZDR32,ZDR33, IZDR11,IZDR12,IZDR13,IZDR21,IZDR22,IZDR23,IZDR31,IZDR32,IZDR33;
  double ZUL11,ZUL12,ZUL13,ZUL21,ZUL22,ZUL23,ZUL31,ZUL32,ZUL33, IZUL11,IZUL12,IZUL13,IZUL21,IZUL22,IZUL23,IZUL31,IZUL32,IZUL33;
  double ZUR11,ZUR12,ZUR13,ZUR21,ZUR22,ZUR23,ZUR31,ZUR32,ZUR33, IZUR11,IZUR12,IZUR13,IZUR21,IZUR22,IZUR23,IZUR31,IZUR32,IZUR33;
  double ZEL11,ZEL12,ZEL13,ZEL21,ZEL22,ZEL23,ZEL31,ZEL32,ZEL33, IZEL11,IZEL12,IZEL13,IZEL21,IZEL22,IZEL23,IZEL31,IZEL32,IZEL33;
  double ZER11,ZER12,ZER13,ZER21,ZER22,ZER23,ZER31,ZER32,ZER33, IZER11,IZER12,IZER13,IZER21,IZER22,IZER23,IZER31,IZER32,IZER33;     		
  double UV11,UV12,UV13, IUV11,IUV12,IUV13;
  double UV21,UV22,UV23, IUV21,IUV22,IUV23;
  double UV31,UV32,UV33, IUV31,IUV32,IUV33;
  double HPP1,HGG1,APP2;
  double HPP2,HGG2,AGG2; 

  double BrMEG ;  //Branching ratio : mu -> e gamma
  double BrTEG ;  //Branching ratio : tau -> e gamma
  double BrTMG ;  //Branching ratio : tau -> mu gamma
  double BrM3e ;  //Branching ratio : mu -> 3e
  double BrT3e ;  //Branching ratio : tau -> 3e 
  double BrT3m ;  //Branching ratio : tau -> 3mu


  double g_2e;  //(g-2)electron 
  double g_2m;  //muon
  double g_2t;  //tau

 double MAh2; 
 double MZ; 




//--------------------------------------------------------------------
// 			SPheno commands
//	     readings of the values in SPheno Output
//--------------------------------------------------------------------


    string name = ""; 
	if(argc == 2){ //Check if there is smthg given in the command
	name = argv[1];}
 
    int n = name.length(); 
  
    // declaring character array 
    char filename[n + 1]; 
  
    // copying the contents of the 
    // string to char array 
    strcpy(filename, name.c_str()); 

    FILE *output = NULL; 
    output = fopen(filename,"r");


   if(output != NULL) {   

   errSph = slhaRead(filename,0);

   Mh1 = slhaVal("MASS",Q,1,25);
   Mh2 = slhaVal("MASS",Q,1,1111);

   //Masses
   MAh2    = slhaVal("MASS",Q,1,200002);
   MetaR  = slhaVal("MASS",Q,1,1001);
   MetaI  = slhaVal("MASS",Q,1,1002);
   MEtac  = slhaVal("MASS",Q,1,1003);   

   MX01 = slhaVal("MASS",Q,1,3001); 
   MX02 = slhaVal("MASS",Q,1,3002); 
   MX03 = slhaVal("MASS",Q,1,3003); 
   //MX04 = slhaVal("MASS",Q,1,3004); 


   mnu1 = slhaVal("MASS",Q,1,12);
   mnu2 = slhaVal("MASS",Q,1,14);
   mnu3 = slhaVal("MASS",Q,1,16);	

   //Flavor violation
   BrMEG = slhaVal("FlavorKitLFV",Q,1,701);
   BrTEG = slhaVal("FlavorKitLFV",Q,1,702); // tau e gamma
   BrTMG = slhaVal("FlavorKitLFV",Q,1,703);
   BrM3e = slhaVal("FlavorKitLFV",Q,1,901); // mu -> 3e
   BrT3e = slhaVal("FlavorKitLFV",Q,1,902); // tau -> 3e 
   BrT3m = slhaVal("FlavorKitLFV",Q,1,903); // tau -> 3mu


   //Scalar parameters
  l1 = slhaVal("MINPAR",Q,1,1);
  lSEt = slhaVal("MINPAR",Q,1,8);
  l4Si = slhaVal("MINPAR",Q,1,6);
  l4Eta = slhaVal("MINPAR",Q,1,2);
  l31EtH = slhaVal("MINPAR",Q,1,3);
  l32EtH = slhaVal("MINPAR",Q,1,4);
  l33EtH = slhaVal("MINPAR",Q,1,5);
  lHSi = slhaVal("MINPAR",Q,1,7);
  lEtSi = slhaVal("MINPAR",Q,1,8);
  vSi = slhaVal("MINPAR",Q,1,11);

  //New Fermion Couplings Real Part
  Yn11 = slhaVal("COUPLINGSYN",Q,2,1,1);
  Yn12 = slhaVal("COUPLINGSYN",Q,2,1,2);
  Yn13 = slhaVal("COUPLINGSYN",Q,2,1,3);

  Yn21 = slhaVal("COUPLINGSYN",Q,2,2,1);
  Yn22 = slhaVal("COUPLINGSYN",Q,2,2,2);
  Yn23 = slhaVal("COUPLINGSYN",Q,2,2,3);
  
  Yn31 = slhaVal("COUPLINGSYN",Q,2,3,1);
  Yn32 = slhaVal("COUPLINGSYN",Q,2,3,2);
  Yn33 = slhaVal("COUPLINGSYN",Q,2,3,3);
  
  //Imaginary part
/*  Yn11c = slhaVal("IMCOUPLINGSYN",Q,2,1,1);
  Yn12c = slhaVal("IMCOUPLINGSYN",Q,2,1,2);
  Yn13c = slhaVal("IMCOUPLINGSYN",Q,2,1,3);

  Yn21c = slhaVal("IMCOUPLINGSYN",Q,2,2,1);
  Yn22c = slhaVal("IMCOUPLINGSYN",Q,2,2,2);
  Yn23c = slhaVal("IMCOUPLINGSYN",Q,2,2,3);
  
  Yn31c = slhaVal("IMCOUPLINGSYN",Q,2,3,1);
  Yn32c = slhaVal("IMCOUPLINGSYN",Q,2,3,2);
  Yn33c = slhaVal("IMCOUPLINGSYN",Q,2,3,3);*/

  
  //New Fermion Couplings Imaginary Part
  Kap11 = slhaVal("KAPPA",Q,2,1,1);
  Kap12 = slhaVal("KAPPA",Q,2,1,2);
  Kap13 = slhaVal("KAPPA",Q,2,1,3);

  Kap21 = slhaVal("KAPPA",Q,2,2,1);
  Kap22 = slhaVal("KAPPA",Q,2,2,2);
  Kap23 = slhaVal("KAPPA",Q,2,2,3);

  Kap31 = slhaVal("KAPPA",Q,2,3,1);
  Kap32 = slhaVal("KAPPA",Q,2,3,2);
  Kap33 = slhaVal("KAPPA",Q,2,3,3);

   
   //Element matrice de melange des scalaires
   ZH_11  = slhaVal("ZH0",Q,2,1,1);
   ZH_12  = slhaVal("ZH0",Q,2,1,2);
   ZH_21  = slhaVal("ZH0",Q,2,2,1);
   ZH_22  = slhaVal("ZH0",Q,2,2,2);
   
/*   ImZH_11 = slhaVal("IMZH0",Q,2,1,1);
   ImZH_12 = slhaVal("IMZH0",Q,2,1,2);
   ImZH_21 = slhaVal("IMZH0",Q,2,2,1);
   ImZH_22 = slhaVal("IMZH0",Q,2,2,2);*/
   
   ZA_11 = slhaVal("ZA0",Q,2,1,1);
   ZA_12 = slhaVal("ZA0",Q,2,1,2);
   ZA_21 = slhaVal("ZA0",Q,2,2,1);
   ZA_22 = slhaVal("ZA0",Q,2,2,2);
 

  //Element matrice de melange des fermions
   ZX_11  = slhaVal("ZX",Q,2,1,1);
   ZX_12  = slhaVal("ZX",Q,2,1,2);
   ZX_13  = slhaVal("ZX",Q,2,1,3);
   ZX_21  = slhaVal("ZX",Q,2,2,1);
   ZX_22  = slhaVal("ZX",Q,2,2,2);
   ZX_23  = slhaVal("ZX",Q,2,2,3);
   ZX_31  = slhaVal("ZX",Q,2,3,1);
   ZX_32  = slhaVal("ZX",Q,2,3,2);
   ZX_33  = slhaVal("ZX",Q,2,3,3);

   ImZX_11  = slhaVal("ZX",Q,2,1,1);
   ImZX_12  = slhaVal("ZX",Q,2,1,2);
   ImZX_13  = slhaVal("ZX",Q,2,1,3);
   ImZX_21  = slhaVal("ZX",Q,2,2,1);
   ImZX_22  = slhaVal("ZX",Q,2,2,2);
   ImZX_23  = slhaVal("ZX",Q,2,2,3);
   ImZX_31  = slhaVal("ZX",Q,2,3,1);
   ImZX_32  = slhaVal("ZX",Q,2,3,2);
   ImZX_33  = slhaVal("ZX",Q,2,3,3);

   //Anomalous magnetic moment
   g_2e = slhaVal("SPhenoLowEnergy",Q,1,20);
   g_2m = slhaVal("SPhenoLowEnergy",Q,1,21);
   g_2t = slhaVal("SPhenoLowEnergy",Q,1,22);

   ZDL11 = slhaVal("UDLMIX",Q,2,1,1) ;
   ZDL12 = slhaVal("UDLMIX",Q,2,1,2) ;
   ZDL13 = slhaVal("UDLMIX",Q,2,1,3); 
   ZDL21 = slhaVal("UDLMIX",Q,2,2,1);
   ZDL22 = slhaVal("UDLMIX",Q,2,2,2); 
   ZDL23 = slhaVal("UDLMIX",Q,2,2,3); 
   ZDL31 = slhaVal("UDLMIX",Q,2,3,1); 
   ZDL32 = slhaVal("UDLMIX",Q,2,3,2); 
   ZDL33 = slhaVal("UDLMIX",Q,2,3,3); 
   ZDR11 = slhaVal("UDRMIX",Q,2,1,1); 
   ZDR12 = slhaVal("UDRMIX",Q,2,1,2); 
   ZDR13 = slhaVal("UDRMIX",Q,2,1,3); 
   ZDR21 = slhaVal("UDRMIX",Q,2,2,1); 
   ZDR22 = slhaVal("UDRMIX",Q,2,2,2); 
   ZDR23 = slhaVal("UDRMIX",Q,2,2,3); 
   ZDR31 = slhaVal("UDRMIX",Q,2,3,1); 
   ZDR32 = slhaVal("UDRMIX",Q,2,3,2); 
   ZDR33 = slhaVal("UDRMIX",Q,2,3,3); 
   ZUL11 = slhaVal("UULMIX",Q,2,1,1); 
   ZUL12 = slhaVal("UULMIX",Q,2,1,2); 
   ZUL13 = slhaVal("UULMIX",Q,2,1,3); 
   ZUL21 = slhaVal("UULMIX",Q,2,2,1); 
   ZUL22 = slhaVal("UULMIX",Q,2,2,2); 
   ZUL23 = slhaVal("UULMIX",Q,2,2,3); 
   ZUL31 = slhaVal("UULMIX",Q,2,3,1); 
   ZUL32 = slhaVal("UULMIX",Q,2,3,2); 
   ZUL33 = slhaVal("UULMIX",Q,2,3,3); 
   ZUR11 = slhaVal("UURMIX",Q,2,1,1); 
   ZUR12 = slhaVal("UURMIX",Q,2,1,2); 
   ZUR13 = slhaVal("UURMIX",Q,2,1,3); 
   ZUR21 = slhaVal("UURMIX",Q,2,2,1); 
   ZUR22 = slhaVal("UURMIX",Q,2,2,2); 
   ZUR23 = slhaVal("UURMIX",Q,2,2,3); 
   ZUR31 = slhaVal("UURMIX",Q,2,3,1); 
   ZUR32 = slhaVal("UURMIX",Q,2,3,2); 
   ZUR33 = slhaVal("UURMIX",Q,2,3,3); 
   ZEL11 = slhaVal("UELMIX",Q,2,1,1); 
   ZEL12 = slhaVal("UELMIX",Q,2,1,2); 
   ZEL13 = slhaVal("UELMIX",Q,2,1,3); 
   ZEL21 = slhaVal("UELMIX",Q,2,2,1); 
   ZEL22 = slhaVal("UELMIX",Q,2,2,2); 
   ZEL23 = slhaVal("UELMIX",Q,2,2,3); 
   ZEL31 = slhaVal("UELMIX",Q,2,3,1); 
   ZEL32 = slhaVal("UELMIX",Q,2,3,2); 
   ZEL33 = slhaVal("UELMIX",Q,2,3,3); 
   ZER11 = slhaVal("UERMIX",Q,2,1,1); 
   ZER12 = slhaVal("UERMIX",Q,2,1,2); 
   ZER13 = slhaVal("UERMIX",Q,2,1,3); 
   ZER21 = slhaVal("UERMIX",Q,2,2,1); 
   ZER22 = slhaVal("UERMIX",Q,2,2,2); 
   ZER23 = slhaVal("UERMIX",Q,2,2,3); 
   ZER31 = slhaVal("UERMIX",Q,2,3,1); 
   ZER32 = slhaVal("UERMIX",Q,2,3,2); 
   ZER33 = slhaVal("UERMIX",Q,2,3,3); 
   UV11  = slhaVal("UVMIX",Q,2,1,1) ;
   UV12  = slhaVal("UVMIX",Q,2,1,2); 
   UV13  = slhaVal("UVMIX",Q,2,1,3); 
   UV21  = slhaVal("UVMIX",Q,2,2,1); 
   UV22  = slhaVal("UVMIX",Q,2,2,2); 
   UV23  = slhaVal("UVMIX",Q,2,2,3);
   UV31  = slhaVal("UVMIX",Q,2,3,1); 
   UV32  = slhaVal("UVMIX",Q,2,3,2);  
   UV33  = slhaVal("UVMIX",Q,2,3,3); 

/*   IZDL11 = slhaVal("IMUDLMIX",Q,2,1,1) ;
   IZDL12 = slhaVal("IMUDLMIX",Q,2,1,2) ;
   IZDL13 = slhaVal("IMUDLMIX",Q,2,1,3); 
   IZDL21 = slhaVal("IMUDLMIX",Q,2,2,1);
   IZDL22 = slhaVal("IMUDLMIX",Q,2,2,2); 
   IZDL23 = slhaVal("IMUDLMIX",Q,2,2,3); 
   IZDL31 = slhaVal("IMUDLMIX",Q,2,3,1); 
   IZDL32 = slhaVal("IMUDLMIX",Q,2,3,2); 
   IZDL33 = slhaVal("IMUDLMIX",Q,2,3,3); 
   IZDR11 = slhaVal("IMUDRMIX",Q,2,1,1); 
   IZDR12 = slhaVal("IMUDRMIX",Q,2,1,2); 
   IZDR13 = slhaVal("IMUDRMIX",Q,2,1,3); 
   IZDR21 = slhaVal("IMUDRMIX",Q,2,2,1); 
   IZDR22 = slhaVal("IMUDRMIX",Q,2,2,2); 
   IZDR23 = slhaVal("IMUDRMIX",Q,2,2,3); 
   IZDR31 = slhaVal("IMUDRMIX",Q,2,3,1); 
   IZDR32 = slhaVal("IMUDRMIX",Q,2,3,2); 
   IZDR33 = slhaVal("IMUDRMIX",Q,2,3,3); 
   IZUL11 = slhaVal("IMUULMIX",Q,2,1,1); 
   IZUL12 = slhaVal("IMUULMIX",Q,2,1,2); 
   IZUL13 = slhaVal("IMUULMIX",Q,2,1,3); 
   IZUL21 = slhaVal("IMUULMIX",Q,2,2,1); 
   IZUL22 = slhaVal("IMUULMIX",Q,2,2,2); 
   IZUL23 = slhaVal("IMUULMIX",Q,2,2,3); 
   IZUL31 = slhaVal("IMUULMIX",Q,2,3,1); 
   IZUL32 = slhaVal("IMUULMIX",Q,2,3,2); 
   IZUL33 = slhaVal("IMUULMIX",Q,2,3,3); 
   IZUR11 = slhaVal("IMUURMIX",Q,2,1,1); 
   IZUR12 = slhaVal("IMUURMIX",Q,2,1,2); 
   IZUR13 = slhaVal("IMUURMIX",Q,2,1,3); 
   IZUR21 = slhaVal("IMUURMIX",Q,2,2,1); 
   IZUR22 = slhaVal("IMUURMIX",Q,2,2,2); 
   IZUR23 = slhaVal("IMUURMIX",Q,2,2,3); 
   IZUR31 = slhaVal("IMUURMIX",Q,2,3,1); 
   IZUR32 = slhaVal("IMUURMIX",Q,2,3,2); 
   IZUR33 = slhaVal("IMUURMIX",Q,2,3,3); 
   IZEL11 = slhaVal("IMUELMIX",Q,2,1,1); 
   IZEL12 = slhaVal("IMUELMIX",Q,2,1,2); 
   IZEL13 = slhaVal("IMUELMIX",Q,2,1,3); 
   IZEL21 = slhaVal("IMUELMIX",Q,2,2,1); 
   IZEL22 = slhaVal("IMUELMIX",Q,2,2,2); 
   IZEL23 = slhaVal("IMUELMIX",Q,2,2,3); 
   IZEL31 = slhaVal("IMUELMIX",Q,2,3,1); 
   IZEL32 = slhaVal("IMUELMIX",Q,2,3,2); 
   IZEL33 = slhaVal("IMUELMIX",Q,2,3,3); 
   IZER11 = slhaVal("IMUERMIX",Q,2,1,1); 
   IZER12 = slhaVal("IMUERMIX",Q,2,1,2); 
   IZER13 = slhaVal("IMUERMIX",Q,2,1,3); 
   IZER21 = slhaVal("IMUERMIX",Q,2,2,1); 
   IZER22 = slhaVal("IMUERMIX",Q,2,2,2); 
   IZER23 = slhaVal("IMUERMIX",Q,2,2,3); 
   IZER31 = slhaVal("IMUERMIX",Q,2,3,1); 
   IZER32 = slhaVal("IMUERMIX",Q,2,3,2); 
   IZER33 = slhaVal("IMUERMIX",Q,2,3,3); 
   IUV11  = slhaVal("IMUVMIX",Q,2,1,1) ;
   IUV12  = slhaVal("IMUVMIX",Q,2,1,2); 
   IUV13  = slhaVal("IMUVMIX",Q,2,1,3); 
   IUV21  = slhaVal("IMUVMIX",Q,2,2,1); 
   IUV22  = slhaVal("IMUVMIX",Q,2,2,2); 
   IUV23  = slhaVal("IMUVMIX",Q,2,2,3);
   IUV31  = slhaVal("IMUVMIX",Q,2,3,1); 
   IUV32  = slhaVal("IMUVMIX",Q,2,3,2);  
   IUV33  = slhaVal("IMUVMIX",Q,2,3,3); */

   HPP1  = slhaVal("EFFHIGGSCOUPLINGS",Q,3,25,22,22); 
   HGG1  = slhaVal("EFFHIGGSCOUPLINGS",Q,3,25,21,21); 
   HPP2  = slhaVal("EFFHIGGSCOUPLINGS",Q,3,1111,22,22);        
   HGG2  = slhaVal("EFFHIGGSCOUPLINGS",Q,3,1111,21,21);
   APP2  = slhaVal("EFFHIGGSCOUPLINGS",Q,3,200002,22,22); 
   AGG2  = slhaVal("EFFHIGGSCOUPLINGS",Q,3,200002,21,21);

	double Md1 = slhaVal("MASS",Q,1,1) ;
	double Md2 = slhaVal("MASS",Q,1,3) ;
	double Md3 = slhaVal("MASS",Q,1,5) ;
	double Mu1 = slhaVal("MASS",Q,1,2) ;
	double Mu2 = slhaVal("MASS",Q,1,4) ;
	double Mu3 = slhaVal("MASS",Q,1,6) ;
	double Me1 = slhaVal("MASS",Q,1,11) ;
	double Me2 = slhaVal("MASS",Q,1,13) ;
	double Me3 = slhaVal("MASS",Q,1,15) ;

//-------Extract PMNS angles-------
//PMNS*Majorana phases = ( A B C
//                         D E F
//                         G H I )
// PMNS = Ue^dagger * Unu

	complex<double> Ue11(ZEL11,IZEL11);
	complex<double> Ue12(ZEL12,IZEL12);
	complex<double> Ue13(ZEL13,IZEL13);
	complex<double> Ue21(ZEL21,IZEL21);
	complex<double> Ue22(ZEL22,IZEL22);
	complex<double> Ue23(ZEL23,IZEL23);
	complex<double> Ue31(ZEL31,IZEL31);
	complex<double> Ue32(ZEL32,IZEL32);
	complex<double> Ue33(ZEL33,IZEL33);

	complex<double> Unu11(UV11,IUV11);
	complex<double> Unu12(UV12,IUV12);
	complex<double> Unu13(UV13,IUV13);
	complex<double> Unu21(UV21,IUV21);
	complex<double> Unu22(UV22,IUV22);
	complex<double> Unu23(UV23,IUV23);
	complex<double> Unu31(UV31,IUV31);
	complex<double> Unu32(UV32,IUV32);
	complex<double> Unu33(UV33,IUV33);

	complex<double> A = conj(Ue11)*Unu11 + conj(Ue21)*Unu21 + conj(Ue31)*Unu31;
	complex<double> D = conj(Ue12)*Unu11 + conj(Ue22)*Unu21 + conj(Ue32)*Unu31;
	complex<double> G = conj(Ue13)*Unu11 + conj(Ue23)*Unu21 + conj(Ue33)*Unu31;
	complex<double> B = conj(Ue11)*Unu12 + conj(Ue21)*Unu22 + conj(Ue31)*Unu32;
	complex<double> E = conj(Ue12)*Unu12 + conj(Ue22)*Unu22 + conj(Ue32)*Unu32;
	complex<double> H = conj(Ue13)*Unu12 + conj(Ue23)*Unu22 + conj(Ue33)*Unu32;
	complex<double> C = conj(Ue11)*Unu13 + conj(Ue21)*Unu23 + conj(Ue31)*Unu33;
	complex<double> F = conj(Ue12)*Unu13 + conj(Ue22)*Unu23 + conj(Ue32)*Unu33;
	complex<double> I = conj(Ue13)*Unu13 + conj(Ue23)*Unu23 + conj(Ue33)*Unu33;

	complex<double> PMNS[3][3] = { {A, B, C},
				       {D, E, F},
				       {G, H, I} };




	//complex<double> nb1(2,3);
	//complex<double> nb2(5,1);
	//cout <<" nb cplx 1 : " << nb1 << " et nb cplex 2 : " << nb2 << endl;
	//cout << "nb1 * nb2 = " << nb1*nb2 << endl;


	complex<double> Theta_13_cpp = asin(abs(PMNS[0][2]))	;
	complex<double> Theta_12_cpp = atan(abs(PMNS[0][1]/PMNS[0][0]));
	complex<double> Theta_23_cpp = atan(abs(PMNS[1][2]/PMNS[2][2]));



   fclose(output); 



//--------------------------------------------------
//                Micromegas Code
//--------------------------------------------------

  // SM particles
  assignValW("Mh1",Mh1);
  assignValW("Mh2",Mh2);
  assignValW("MZ",91.2);
  assignValW("Mu1",Mu1);
  assignValW("Mu2",1.3);
  assignValW("Mu3",173.2);
  assignValW("Md1",Md1);
  assignValW("Md2",Md2);
  assignValW("Md3",4.2);
  assignValW("Me1",Me1);
  assignValW("Me2",Me2);
  assignValW("Me3",1.7);
  assignValW("Mnu1",mnu1);
  assignValW("Mnu2",mnu2);
  assignValW("Mnu3",mnu3);


  // New particles
  assignValW("MetaR",MetaR);
  assignValW("MAh2",MAh2);
  assignValW("MX01",MX01);
  assignValW("MX02",MX02);
  assignValW("MX03",MX03);
  assignValW("Metp",MEtac);
  assignValW("MetaI",MetaI);

  
  // New coupling Real Part
  assignValW("Yn11",Yn11);
  assignValW("Yn12",Yn12);
  assignValW("Yn13",Yn13);

  assignValW("Yn21",Yn21);
  assignValW("Yn22",Yn22);
  assignValW("Yn23",Yn23);
  
  assignValW("Yn31",Yn31);
  assignValW("Yn32",Yn32);
  assignValW("Yn33",Yn33);


  assignValW("Kap11",Kap11);
  assignValW("Kap12",Kap12);  
  assignValW("Kap13",Kap13);  
 
  assignValW("Kap21",Kap21);
  assignValW("Kap22",Kap22);  
  assignValW("Kap23",Kap23);  
  
  assignValW("Kap31",Kap31);
  assignValW("Kap32",Kap32);  
  assignValW("Kap33",Kap33);  
 
/*  
  // New coupling Imaginary Part
  assignValW("IYn11",Yn11c);
  assignValW("IYn12",Yn12c);
  assignValW("IYn13",Yn13c);

  assignValW("IYn21",Yn21c);
  assignValW("IYn22",Yn22c);
  assignValW("IYn23",Yn23c);
  
  assignValW("IYn31",Yn31c);
  assignValW("IYn32",Yn32c);
  assignValW("IYn33",Yn33c);

  assignValW("IKap11",Kap11c);
  assignValW("IKap12",Kap12c);  
  assignValW("IKap13",Kap13c);  
 
  assignValW("IKap21",Kap21c);
  assignValW("IKap22",Kap22c);  
  assignValW("IKap23",Kap23c);  
  
  assignValW("IKap31",Kap31c);
  assignValW("IKap32",Kap32c);  
  assignValW("IKap33",Kap33c);  */

  //Scalar parameters
  assignValW("L4Et",l4Eta);
  //assignValW("IL4Et",0.0);
  assignValW("L4S",l4Si);
  assignValW("LHE32",l32EtH);
  //assignValW("ILHE32",0.0);
  assignValW("LHE31",l31EtH);
  //assignValW("ILHE31",0.0);
  assignValW("LHE33",l33EtH);
  assignValW("LHS",lHSi);
  assignValW("LES",lEtSi);  
  assignValW("vSi",vSi);  
 

  //Matrices
  assignValW("ZX11",ZX_11);
  assignValW("ZX12",ZX_12);
  assignValW("ZX13",ZX_13);
  assignValW("ZX21",ZX_21);
  assignValW("ZX22",ZX_22);
  assignValW("ZX23",ZX_23);
  assignValW("ZX31",ZX_31);
  assignValW("ZX32",ZX_32);
  assignValW("ZX33",ZX_33);

/*  assignValW("IZX11",ImZX_11);
  assignValW("IZX12",ImZX_12);
  assignValW("IZX13",ImZX_13);
  assignValW("IZX21",ImZX_21);
  assignValW("IZX22",ImZX_22);
  assignValW("IZX23",ImZX_23);
  assignValW("IZX31",ImZX_31);
  assignValW("IZX32",ImZX_32);
  assignValW("IZX33",ImZX_33);*/


  assignValW("ZA11",ZA_11);
  assignValW("ZA12",ZA_12);
  assignValW("ZA21",ZA_21);
  assignValW("ZA22",ZA_22);
  assignValW("ZH11",ZH_11);
  assignValW("ZH12",ZH_12);
  assignValW("ZH21",ZH_21);
  assignValW("ZH22",ZH_22);

  assignValW("ZDL11",ZDL11);
  assignValW("ZDL12",ZDL12);
  assignValW("ZDL13",ZDL13);
  assignValW("ZDL21",ZDL21);
  assignValW("ZDL22",ZDL22);
  assignValW("ZDL23",ZDL23);
  assignValW("ZDL31",ZDL31);
  assignValW("ZDL32",ZDL32);
  assignValW("ZDL33",ZDL33);
  assignValW("ZDR11",ZDR11);
  assignValW("ZDR12",ZDR12);
  assignValW("ZDR13",ZDR13);
  assignValW("ZDR21",ZDR21);
  assignValW("ZDR22",ZDR22);
  assignValW("ZDR23",ZDR23);
  assignValW("ZDR31",ZDR31);
  assignValW("ZDR32",ZDR32);
  assignValW("ZDR33",ZDR33);
  assignValW("ZUL11",ZUL11);
  assignValW("ZUL12",ZUL12);
  assignValW("ZUL13",ZUL13);
  assignValW("ZUL21",ZUL21);
  assignValW("ZUL22",ZUL22);
  assignValW("ZUL23",ZUL23);
  assignValW("ZUL31",ZUL31);
  assignValW("ZUL32",ZUL32);
  assignValW("ZUL33",ZUL33);
  assignValW("ZUR11",ZUR11);
  assignValW("ZUR12",ZUR12);
  assignValW("ZUR13",ZUR13);
  assignValW("ZUR21",ZUR21);
  assignValW("ZUR22",ZUR22);
  assignValW("ZUR23",ZUR23);
  assignValW("ZUR31",ZUR31);
  assignValW("ZUR32",ZUR32);
  assignValW("ZUR33",ZUR33);
  assignValW("ZEL11",ZEL11);
  assignValW("ZEL12",ZEL12);
  assignValW("ZEL13",ZEL13);
  assignValW("ZEL21",ZEL21);
  assignValW("ZEL22",ZEL22);
  assignValW("ZEL23",ZEL23);
  assignValW("ZEL31",ZEL31);
  assignValW("ZEL32",ZEL32);
  assignValW("ZEL33",ZEL33);
  assignValW("ZER11",ZER11);
  assignValW("ZER12",ZER12);
  assignValW("ZER13",ZER13);
  assignValW("ZER21",ZER21);
  assignValW("ZER22",ZER22);
  assignValW("ZER23",ZER23);
  assignValW("ZER31",ZER31);
  assignValW("ZER32",ZER32);
  assignValW("ZER33",ZER33);
  assignValW("UV11",UV11);
  assignValW("UV12",UV12);
  assignValW("UV13",UV13);
  assignValW("UV21",UV21);
  assignValW("UV22",UV22);
  assignValW("UV23",UV23);
  assignValW("UV31",UV31);
  assignValW("UV32",UV32);
  assignValW("UV33",UV33);

/*  assignValW("IZDL11",IZDL11);
  assignValW("IZDL12",IZDL12);
  assignValW("IZDL13",IZDL13);
  assignValW("IZDL21",IZDL21);
  assignValW("IZDL22",IZDL22);
  assignValW("IZDL23",IZDL23);
  assignValW("IZDL31",IZDL31);
  assignValW("IZDL32",IZDL32);
  assignValW("IZDL33",IZDL33);
  assignValW("IZDR11",IZDR11);
  assignValW("IZDR12",IZDR12);
  assignValW("IZDR13",IZDR13);
  assignValW("IZDR21",IZDR21);
  assignValW("IZDR22",IZDR22);
  assignValW("IZDR23",IZDR23);
  assignValW("IZDR31",IZDR31);
  assignValW("IZDR32",IZDR32);
  assignValW("IZDR33",IZDR33);
  assignValW("IZUL11",IZUL11);
  assignValW("IZUL12",IZUL12);
  assignValW("IZUL13",IZUL13);
  assignValW("IZUL21",IZUL21);
  assignValW("IZUL22",IZUL22);
  assignValW("IZUL23",IZUL23);
  assignValW("IZUL31",IZUL31);
  assignValW("IZUL32",IZUL32);
  assignValW("IZUL33",IZUL33);
  assignValW("IZUR11",IZUR11);
  assignValW("IZUR12",IZUR12);
  assignValW("IZUR13",IZUR13);
  assignValW("IZUR21",IZUR21);
  assignValW("IZUR22",IZUR22);
  assignValW("IZUR23",IZUR23);
  assignValW("IZUR31",IZUR31);
  assignValW("IZUR32",IZUR32);
  assignValW("IZUR33",IZUR33);
  assignValW("IZEL11",IZEL11);
  assignValW("IZEL12",IZEL12);
  assignValW("IZEL13",IZEL13);
  assignValW("IZEL21",IZEL21);
  assignValW("IZEL22",IZEL22);
  assignValW("IZEL23",IZEL23);
  assignValW("IZEL31",IZEL31);
  assignValW("IZEL32",IZEL32);
  assignValW("IZEL33",IZEL33);
  assignValW("IZER11",IZER11);
  assignValW("IZER12",IZER12);
  assignValW("IZER13",IZER13);
  assignValW("IZER21",IZER21);
  assignValW("IZER22",IZER22);
  assignValW("IZER23",IZER23);
  assignValW("IZER31",IZER31);
  assignValW("IZER32",IZER32);
  assignValW("IZER33",IZER33);
  assignValW("IUV11",IUV11);
  assignValW("IUV12",IUV12);
  assignValW("IUV13",IUV13);
  assignValW("IUV21",IUV21);
  assignValW("IUV22",IUV22);
  assignValW("IUV23",IUV23);
  assignValW("IUV31",IUV31);
  assignValW("IUV32",IUV32);*/
 
  assignValW("HPP1",HPP1);
  assignValW("HGG1",HGG1);
  assignValW("HPP2",HPP2);
  assignValW("HGG2",HGG2);
  assignValW("APP2",APP2);
  assignValW("AGG2",AGG2);

//---------------------------------------------------------------
//             Prepare relic density calculation
//---------------------------------------------------------------
                 err=sortOddParticles(cdmName);
      
     printf("MAh2 = %e",MAh2 );	
     //printf("Mass of X01 = %e",pMass("~X01"));
     qNumbers(CDM1, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E and charge %d \n",CDM1,  spin2,Mcdm1, charge3); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n"); 
  
  	int chargeDM; //DM should be a neutral particle
	if(Mcdm1 == pMass("~etaI") || Mcdm1 == pMass("~etaR") || Mcdm1 == pMass("~X01")){
	chargeDM = 1;}
	else{
	chargeDM = 100;};

    //Compute DM relic density
	printf("\n==== Calculation of relic density =====\n");   
	double Omegah2;
 	double Beps= 1.E-4;
 	int fast=1; 
	double Xf;
    double cut=0.01;
	qNumbers(CDM1, &spin2, &charge3, &cdim);
    Omegah2=darkOmega(&Xf,fast,Beps,&err);
    printf("Xf = %lf       Omegah2 = %lf       \n", Xf, Omegah2 );
   if(Omegah2>0)printChannels(Xf,cut,Beps,1,stdout);




    printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
    printHiggs(stdout);
    printMasses(stdout,1);


/*

//----------------
//    X01 DM
//----------------

double x01x01_nimp = oneChannel(Xf, Beps, "~FX01", "~FX01", "*", "*"); //somme toutes les possibilites de X01 X01 -> Tout possible
double xpxm_nimp   = oneChannel(Xf, Beps, "~FXm", "~fXm", "*", "*");
double x01x02_nimp = oneChannel(Xf, Beps, "~FX01", "~FX02", "*", "*");
double xpx01_nimp = oneChannel(Xf, Beps, "~fXm", "~FX01", "*", "*");
double xmx01_nimp = oneChannel(Xf, Beps, "~FXm", "~FX01", "*", "*");
double xmxm_nimp   = oneChannel(Xf, Beps, "~FXm", "~FXm", "*", "*");
double xpxp_nimp   = oneChannel(Xf, Beps, "~fXm", "~fXm", "*", "*");

//printf("Test X01 X01 -> Tout  : %e  \n",x01x01_nimp);

double xpxm_e1E1 = oneChannel(Xf, Beps, "~FXm", "~fXm", "e1", "E1");
double xpxm_e2E2 = oneChannel(Xf, Beps, "~FXm", "~fXm", "e2", "E2");
double xpxm_e3E3 = oneChannel(Xf, Beps, "~FXm", "~fXm", "e3", "E3");

double xpx01_d1U1 = oneChannel(Xf, Beps, "~FXm", "~FX01", "d1", "U1");
double xpx01_d2U2 = oneChannel(Xf, Beps, "~FXm", "~FX01", "d2", "U2");
double xpx01_d3U3 = oneChannel(Xf, Beps, "~FXm", "~FX01", "d3", "U3");

double x01x02_e1E1 = oneChannel(Xf, Beps, "~FX01", "~FX02", "e1", "E1");
double x01x02_e2E2 = oneChannel(Xf, Beps, "~FX01", "~FX02", "e2", "E2");
double x01x02_e3E3 = oneChannel(Xf, Beps, "~FX01", "~FX02", "e3", "E3");

//LFV
double x01x02_e1E2 = oneChannel(Xf, Beps, "~FX01", "~FX02", "e1", "E2");
double x01x02_e1E3 = oneChannel(Xf, Beps, "~FX01", "~FX02", "e1", "E3");

double x01x02_e2E1 = oneChannel(Xf, Beps, "~FX01", "~FX02", "e2", "E1");
double x01x02_e2E3 = oneChannel(Xf, Beps, "~FX01", "~FX02", "e2", "E3");

double x01x02_e3E1 = oneChannel(Xf, Beps, "~FX01", "~FX02", "e3", "E1");
double x01x02_e3E2 = oneChannel(Xf, Beps, "~FX01", "~FX02", "e3", "E2");

double xpxp_wmwm = oneChannel(Xf, Beps, "~FXm", "~FXm", "Wm", "Wm");
double xpxp_wpwp = oneChannel(Xf, Beps, "~FXm", "~FXm", "Wp", "Wp");

double x01x01_hz = oneChannel(Xf, Beps, "~FX01", "~FX01", "h", "Z");
double x01x01_hh = oneChannel(Xf, Beps, "~FX01", "~FX01", "h", "h");
double x01x01_wpwm = oneChannel(Xf, Beps, "~FX01", "~FX01", "Wp", "Wm");
double x01x01_e1E1 = oneChannel(Xf, Beps, "~FX01", "~FX01", "e1", "E1");
double x01x01_e2E2 = oneChannel(Xf, Beps, "~FX01", "~FX01", "e2", "E2");
double x01x01_e3E3 = oneChannel(Xf, Beps, "~FX01", "~FX01", "e3", "E3");
double x01x01_n1n1 = oneChannel(Xf, Beps, "~FX01", "~FX01", "nu1", "nu1");
double x01x01_n2n2 = oneChannel(Xf, Beps, "~FX01", "~FX01", "nu2", "nu2");
double x01x01_n3n3 = oneChannel(Xf, Beps, "~FX01", "~FX01", "nu3", "nu3");*/

/*
double totx = x01x02_nimp + x01x01_nimp + xpxm_nimp + xpx01_nimp + xmx01_nimp + xpxp_nimp + xmxm_nimp  ;
//----------------
//    P01 DM
//----------------

double p01p01_nimp = oneChannel(Xf, Beps, "~P01", "~P01", "*", "*"); //somme toutes les possibilites de X01 X01 -> Tout possible
double pppm_nimp   = oneChannel(Xf, Beps, "~Phip", "~Phim", "*", "*");
double p01a0_nimp  = oneChannel(Xf, Beps, "~P01", "~A0", "*", "*");
double a0a0_nimp   = oneChannel(Xf, Beps, "~A0", "~A0", "*", "*");

double pppp_nimp   = oneChannel(Xf, Beps, "~Phip", "~Phip", "*", "*");
double pmpm_nimp   = oneChannel(Xf, Beps, "~Phim", "~Phim", "*", "*");

double ppa0_nimp   = oneChannel(Xf, Beps, "~Phip", "~A0", "*", "*");
double pma0_nimp   = oneChannel(Xf, Beps, "~Phim", "~A0", "*", "*");
double ppp01_nimp   = oneChannel(Xf, Beps, "~Phip", "~P01", "*", "*");
double pmp01_nimp   = oneChannel(Xf, Beps, "~Phim", "~P01", "*", "*");



double p01p01_wpwm = oneChannel(Xf, Beps, "~P01", "~P01", "Wp", "Wm");
double p01p01_zz = oneChannel(Xf, Beps, "~P01", "~P01", "Z", "Z");
double p01p01_hwp = oneChannel(Xf, Beps, "~P01", "~P01", "h", "Wp");
double p01p01_hh = oneChannel(Xf, Beps, "~P01", "~P01", "h", "h");

double p01a0_hZ = oneChannel(Xf, Beps, "~P01", "~A0", "h", "Z");
double ppa0_wpZ = oneChannel(Xf, Beps, "~Phip", "~A0", "Wp", "Z");
double a0a0_wpwm = oneChannel(Xf, Beps, "~A0", "~A0", "Wp", "Wm");

double pppp_wpwp = oneChannel(Xf, Beps, "~Phip", "~Phip", "Wp", "Wp");
double pppm_zz = oneChannel(Xf, Beps, "~Phip", "~Phim", "Z", "Z");
double pppm_hh = oneChannel(Xf, Beps, "~Phip", "~Phim", "h", "h");

double totp = a0a0_nimp + p01a0_nimp + pppm_nimp + pppp_nimp + p01p01_nimp + ppa0_nimp + pma0_nimp + pmpm_nimp + pmp01_nimp + ppp01_nimp;

double TheTot = totp + totx ;
//----------------
//     A0 DM
//----------------
*/


//------------Prepare calculation of direct detection cross-section and limit---------------
//------------------------------------------------------------------------------------------
 
  double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

  double dNdE[300];
  double nEvents;

 	 //printf("\n======== Direct Detection ========\n");    

 	nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,dNdE);
	nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
       // printf("\nCDM[antiCDM]-nucleon micrOMEGAs amplitudes for %s \n",CDM1);
       // printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
       // printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] );

        SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
       // printf("\nCDM[antiCDM]-nucleon cross sections[pb]:\n");
       // printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       //    SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
       // printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       // SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);

        double sigSIP, sigSIN, sigSDP, sigSDN;
        sigSIP = SCcoeff*pA0[0]*pA0[0]; //cross-section proton, Spin Independent
        sigSIN = SCcoeff*nA0[0]*nA0[0]; //cross-section neutron, Spin Independent
        sigSDP = 3*SCcoeff*pA5[0]*pA5[0]; //cross-section proton, Spin Dependent
        sigSDN = 3*SCcoeff*nA5[0]*nA5[0]; //cross-section neutron, Spin Dependent

       	//Conversion pb -> cm2 : 1pb = 1.0*10^-36 cm2
	double sigSIP_cm2, sigSIN_cm2, sigSDP_cm2, sigSDN_cm2; //cm2
	sigSIP_cm2 = sigSIP * 1.0*pow(10,-36);
	sigSIN_cm2 = sigSIN * 1.0*pow(10,-36);
	sigSDP_cm2 = sigSDP * 1.0*pow(10,-36);
	sigSDN_cm2 = sigSDN * 1.0*pow(10,-36);




//---------------------------------------------------------------------------------------------------------
//-----------------------------------Print the values into files-------------------------------------------
//---------------------------------------------------------------------------------------------------------
//---------------------------------Test of NAN---------------------------------------
//---------------------------------0 no nan, -1 nan----------------------------------

	//mA0_tree = sqrt(MPhi2 + (lP+lPp-lPpp)*pow(v,2)/2);
   	//mPc_tree = sqrt(2*MPhi2 + (lP)*pow(v,2));
   	//mP0_tree = sqrt(MPhi2 + (lP+lPp+lPpp)*pow(v,2)/2)   ;
   	//mS_tree  = sqrt(MS2 + (lS)*pow(v,2)/2) ;
	
	

	FILE *SPhenoend = NULL;
	SPhenoend = fopen(filename,"a");
	fprintf(SPhenoend ," \nBlock RELIC DENSITY    #from micromegas \n");
	fprintf(SPhenoend ,"   1    %lf    #Omega_h2 \n",Omegah2);
	fprintf(SPhenoend ,"   2    %d     #Test for the charge of the DM 1=neutral, 100=not neutral \n",chargeDM);
	fprintf(SPhenoend ,"   3    %lf    #Mass of DM \n",Mcdm1);
	fprintf(SPhenoend ," \nBlock DIRECT DETECTION  #from micromegas in [cm2] \n");
	fprintf(SPhenoend ,"   1      %e   #Proton Spin Dependent \n",sigSDP_cm2);
	fprintf(SPhenoend ,"   2      %e   #Proton Spin Independent \n",sigSIP_cm2);
	fprintf(SPhenoend ,"   3      %e   #Neutron Spin Dependent \n",sigSDN_cm2);
	fprintf(SPhenoend ,"   4      %e   #Neutron Spin Independent \n",sigSIN_cm2);
	/*fprintf(SPhenoend ," \nBlock CANALS FERMION  #from micromegas  \n");
	fprintf(SPhenoend ,"   1    %e     #X01 X01 -> ~ ~ \n",x01x01_nimp);
	fprintf(SPhenoend ,"   2    %e     #X01 X02 -> ~ ~ \n",x01x02_nimp);
	fprintf(SPhenoend ,"   3    %e     #X^+ X^- -> ~ ~ \n",xpxm_nimp);
	fprintf(SPhenoend ,"   4    %e     #X^+ X01- -> ~ ~ \n",xpx01_nimp);
	fprintf(SPhenoend ,"   5    %e     #X^- X01 -> ~ ~ \n",xmx01_nimp);
	fprintf(SPhenoend ,"   6    %e     #X^- X^- -> ~ ~ \n",xmxm_nimp);
	fprintf(SPhenoend ,"   7    %e     #X^+ X^+ -> ~ ~ \n",xpxp_nimp);
	fprintf(SPhenoend ,"   8    %e    #Total \n",totx);
	fprintf(SPhenoend ," \nBlock CANALS SCALAR  #from micromegas  \n");
	fprintf(SPhenoend ,"   1    %e      #A0 A0   -> ~ ~ \n",a0a0_nimp);
	fprintf(SPhenoend ,"   2    %e      #P01 A0  -> ~ ~  \n",p01a0_nimp);
	fprintf(SPhenoend ,"   3    %e      #P^+ P^- -> ~ ~  \n",pppm_nimp);
	fprintf(SPhenoend ,"   4    %e      #P^+ P^+ -> ~ ~  \n",pppp_nimp);
	fprintf(SPhenoend ,"   5    %e      #P^- P^- -> ~ ~  \n",pmpm_nimp);
	fprintf(SPhenoend ,"   6    %e      #P01 P01 -> ~ ~  \n",p01p01_nimp);
	fprintf(SPhenoend ,"   7    %e      #P^- P01 -> ~ ~  \n",pmp01_nimp);
	fprintf(SPhenoend ,"   8    %e      #P^+ P01 -> ~ ~  \n",ppp01_nimp);
	fprintf(SPhenoend ,"   9    %e      #P^+  A0 -> ~ ~  \n",ppa0_nimp);
	fprintf(SPhenoend ,"   10    %e      #P^-  A0 -> ~ ~  \n",pma0_nimp);
	fprintf(SPhenoend ,"   11   %e      #Tot   \n",totp);*/
	fclose(SPhenoend);  




 };  //enf if output non null

	  


  
  // End of program
  return 0;
}

