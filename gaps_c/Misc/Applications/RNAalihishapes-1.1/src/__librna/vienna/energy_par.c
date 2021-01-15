/* 
    Current free energy parameters are summarized in:

    D.H.Mathews, J. Sabina, M. ZUker, D.H. Turner
    "Expanded sequence dependence of thermodynamic parameters improves
    prediction of RNA secondary structure"
    JMB, 288, pp 911-940, 1999

    Enthalpies taken from:
    
    A. Walter, D Turner, J Kim, M Lyttle, P M"uller, D Mathews, M Zuker
    "Coaxial stckaing of helices enhances binding of oligoribonucleotides.."
    PNAS, 91, pp 9218-9222, 1994
    
    D.H. Turner, N. Sugimoto, and S.M. Freier.
    "RNA Structure Prediction",
    Ann. Rev. Biophys. Biophys. Chem. 17, 167-192, 1988.

    John A.Jaeger, Douglas H.Turner, and Michael Zuker.
    "Improved predictions of secondary structures for RNA",
    PNAS, 86, 7706-7710, October 1989.
    
    L. He, R. Kierzek, J. SantaLucia, A.E. Walter, D.H. Turner
    "Nearest-Neughbor Parameters for GU Mismatches...."
    Biochemistry 1991, 30 11124-11132

    A.E. Peritz, R. Kierzek, N, Sugimoto, D.H. Turner
    "Thermodynamic Study of Internal Loops in Oligoribonucleotides..."
    Biochemistry 1991, 30, 6428--6435

    
*/

#include "energy_const.h"
/*@unused@*/
static char rcsid[] = "$Id: energy_par.c,v 1.6 2004/08/12 12:11:57 ivo Exp $";

#define NST 0     /* Energy for nonstandard stacked pairs */
#define DEF -50   /* Default terminal mismatch, used for I */
                  /* and any non_pairing bases */
#define NSM 0     /* terminal mismatch for non standard pairs */

#define PUBLIC

PUBLIC double Tmeasure = 37+K0;  /* temperature of param measurements */




//@@ stack[block name in file] --> S --> stack37[data_structure_name]
PUBLIC int stack37[NBPAIRS+1][NBPAIRS+1] =
/*          CG     GC     GU     UG     AU     UA  */
{
  {  INF,   INF,   INF,   INF,   INF,   INF,   INF, INF},
    { INF,-240,-330,-210,-140,-210,-210, 0},
    { INF,-330,-340,-250,-150,-220,-240, 0},
    { INF,-210,-250, 130, DEF,-140,-130, 0},
    { INF,-140,-150, DEF,30, -60,-100, 0},
    { INF,-210,-220,-140, -60,-110, -90, 0},
    { INF,-210,-240,-130,-100, -90,-130, 0},
    { INF, 0, 0, 0, 0, 0, 0, 0}
};
// { {  INF,   INF,   INF,   INF,   INF,   INF,   INF, INF},
//   {  INF,  -240,  -330,  -210,  -140,  -210,  -210, NST},
//   {  INF,  -330,  -340,  -250,  -150,  -220,  -240, NST},
//   {  INF,  -210,  -250,   130,   -50,  -140,  -130, NST},
//   {  INF,  -140,  -150,   -50,    30,   -60,  -100, NST},
//   {  INF,  -210,  -220,  -140,   -60,  -110,   -90, NST},
//   {  INF,  -210,  -240,  -130,  -100,   -90,  -130, NST},
//   {  INF,   NST,   NST,   NST,   NST,   NST,   NST, NST}};

//@@ stack_enthalpies[block name in file] --> S_H --> stackdH[data_structure_name]
/* enthalpies (0.01*kcal/mol at 37 C) for stacked pairs */
/* different from mfold-2.3, which uses values from mfold-2.2 */
PUBLIC int stackdH[NBPAIRS+1][NBPAIRS+1] = 
/*          CG     GC     GU     UG     AU     UA  */
{
  {  INF,   INF,   INF,   INF,   INF,   INF,   INF, INF},
    { INF,-1060 -1340 -1210,-560 -1050 -1040, 0},
    { INF,-1340 -1490 -1260,-830 -1140 -1240, 0},
    { INF,-1210 -1260 -1460 -1350,-880 -1280, 0},
    { INF,-560,-830 -1350,-930,-320,-700, 0},
    { INF,-1050 -1140,-880,-320,-940,-680, 0},
    { INF,-1040 -1240 -1280,-700,-680,-770, 0},
    { INF, 0, 0, 0, 0, 0, 0, 0}
};
// { {  INF,   INF,   INF,   INF,   INF,   INF,   INF, INF}, 
//   {  INF, -1060, -1340, -1210,  -560, -1050, -1040, NST},
//   {  INF, -1340, -1490, -1260,  -830, -1140, -1240, NST},
//   {  INF, -1210, -1260, -1460, -1350,  -880, -1280, NST},
//   {  INF,  -560,  -830, -1350,  -930,  -320,  -700, NST},
//   {  INF, -1050, -1140,  -880,  -320,  -940,  -680, NST},
//   {  INF, -1040, -1240, -1280,  -700,  -680,  -770, NST},
//   {  INF,   NST,   NST,   NST,   NST,   NST,   NST, NST}};

//@@ hairpin[block name in file] --> HP --> hairpin37[data_structure_name] 
PUBLIC int hairpin37[31] = {
  INF, INF, INF, 570, 560, 560, 540, 590, 560, 640, 650,
       660, 670, 678, 686, 694, 701, 707, 713, 719, 725,
       730, 735, 740, 744, 749, 753, 757, 761, 765, 769};

//@@ hairpin_enthalpies[block name in file] --> HP_H --> hairpindH[data_structure_name]
PUBLIC int hairpindH[31] = {
  INF, INF, INF, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

//@@ bulge [in file] --> B --> bulge37 [in program]
PUBLIC int bulge37[31] = {
  INF, 380, 280, 320, 360, 400, 440, 459, 470, 480, 490,
       500, 510, 519, 527, 534, 541, 548, 554, 560, 565,
  571, 576, 580, 585, 589, 594, 598, 602, 605, 609};

//@@ bulge_enthalpies [in file] --> B_H --> bulgedH [in program]
PUBLIC int bulgedH[31] = {
  INF, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

//@@ interior [in file] --> IL --> internal_loop37 [in program]
PUBLIC int internal_loop37[31] = {
  INF, INF, 410, 510, 170, 180, 200, 220, 230, 240, 250,
       260, 270, 278, 286, 294, 301, 307, 313, 319, 325,
       330, 335, 340, 345, 349, 353, 357, 361, 365, 369};
       
       
//@@ interior_enthalpies [in file] --> IL_H --> internal_loopdH [in program]
PUBLIC int internal_loopdH[31] = {
  INF, INF, INF, INF, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0};



       
       
       
//@@ mismatch_exterior [in file] --> MME --> mismatchExt37 [in program]       
PUBLIC int mismatchExt37[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0,-110, -40,-130, -60},
    { DEF,-160, -90,-180,-110},
    { -30,-140, -70,-160, -90},
    { -20,-130, -60,-150, -80},
    { -10,-120, DEF,-140, -70}},
  { /* GC */
    { 0,-170, -80,-170,-120},
    { -20,-190,-100,-190,-140},
    { -30,-200,-110,-200,-150},
    { 0,-170, -80,-170,-120},
    { 0,-170, -80,-170,-120}},
  { /* GU */
    { 0, -70, -10, -70, -10},
    { -30,-100, -40,-100, -40},
    { -30,-100, -40,-100, -40},
    { -40,-110, DEF,-110, DEF},
    { -20, -90, -30, -90, -30}},
  { /* UG */
    { 0, -80, DEF, -80, -60},
    { -30,-110, -80,-110, -90},
    { -10, -90, -60, -90, -70},
    { -20,-100, -70,-100, -80},
    { -20,-100, -70,-100, -80}},
  { /* AU */
    { 0, -70, -10, -70, -10},
    { -30,-100, -40,-100, -40},
    { -30,-100, -40,-100, -40},
    { -40,-110, DEF,-110, DEF},
    { -20, -90, -30, -90, -30}},
  { /* UA */
    { 0, -80, DEF, -80, -60},
    { -30,-110, -80,-110, -90},
    { -10, -90, -60, -90, -70},
    { -20,-100, -70,-100, -80},
    { -20,-100, -70,-100, -80}},
  { /* @@ */
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0}}
};
       
       
//@@ mismatch_exterior_enthalpies [in file] --> MME_H --> mismatchExtdH [in program]       
PUBLIC int mismatchExtdH[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0,-740,-280,-640,-360},
    {-240,-980,-520,-880,-600},
    { 330,-410,50,-310, -30},
    {80,-660,-200,-560,-280},
    {-140,-880,-420,-780,-500}},
  { /* GC */
    { 0,-900,-410,-860,-750},
    {-160 -1060,-570 -1020,-910},
    {70,-830,-340,-790,-680},
    {-460 -1360,-870 -1320 -1210},
    { -40,-940,-450,-900,-790}},
  { /* GU */
    { 0,-740,-240,-720,-490},
    { 160,-580, -80,-560,-330},
    { 220,-520, -20,-500,-270},
    {70,-670,-170,-650,-420},
    { 310,-430,70,-410,-180}},
  { /* UG */
    { 0,-490, -90,-550,-230},
    {-150,-640,-240,-700,-380},
    { 510,20, 420, -40, 280},
    {10,-480, -80,-540,-220},
    { 100,-390,10,-450,-130}},
  { /* AU */
    { 0,-570, -70,-580,-220},
    { 160,-410,90,-420, -60},
    { 220,-350, 150,-360, 0},
    {70,-500, 0,-510,-150},
    { 310,-260, 240,-270,90}},
  { /* UA */
    { 0,-490, -90,-550,-230},
    { DEF,-540,-140,-600,-280},
    { 690, 200, 600, 140, 460},
    { -60,-550,-150,-610,-290},
    { -60,-550,-150,-610,-290}},
  { /* @@ */
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0}}
};



// mismatch_hairpin ==> MMH  ==> mismatchH37
/* mismatch free energies for hairpins at 37C */
PUBLIC int mismatchH37[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0, 0, 0, 0, 0},
    { -90,-150,-150,-140,-180},
    { -90,-100, -90,-290, -80},
    { -90,-220,-200,-160,-110},
    { -90,-170,-140,-180,-200}},
  { /* GC */
    { 0, 0, 0, 0, 0},
    { -70,-110,-150,-130,-210},
    { -70,-110, -70,-240, DEF},
    { -70,-240,-290,-140,-120},
    { -70,-190,-100,-220,-150}},
  { /* GU */
    { 0, 0, 0, 0, 0},
    { 0,20, DEF, -30, -30},
    { 0, -10, -20,-150, -20},
    { 0, -90,-110, -30, 0},
    { 0, -30, -30, -40,-110}},
  { /* UG */
    { 0, 0, 0, 0, 0},
    { 0, DEF, -30, -60, DEF},
    { 0, -20, -10,-170, 0},
    { 0, -80,-120, -30, -70},
    { 0, -60, -10, -60, -80}},
  { /* AU */
    { 0, 0, 0, 0, 0},
    { 0, -30, DEF, -30, -30},
    { 0, -10, -20,-150, -20},
    { 0,-110,-120, -20,20},
    { 0, -30, -30, -60,-110}},
  { /* UA */
    { 0, 0, 0, 0, 0},
    { 0, DEF, -30, -60, DEF},
    { 0, -20, -10,-120, 0},
    { 0,-140,-120, -70, -20},
    { 0, -30, -10, DEF, -80}},
  { /* @@ */
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0}}
};


// mismatch_hairpin_enthalpies ==> MMH_H  ==> mismatchHdH
/* mismatch free energies for hairpins at 37C */
PUBLIC int mismatchHdH[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0, 0, 0, 0, 0},
    { DEF -1030,-950 -1030 -1030},
    { DEF,-520,-450,-520,-670},
    { DEF,-940,-940,-940,-940},
    { DEF,-810,-740,-810,-860}},
  { /* GC */
    { 0, 0, 0, 0, 0},
    { DEF,-520,-880,-560,-880},
    { DEF,-720,-310,-310,-390},
    { DEF,-710,-740,-620,-740},
    { DEF,-500,-500,-500,-570}},
  { /* GU */
    { 0, 0, 0, 0, 0},
    { DEF,-430,-600,-600,-600},
    { DEF,-260,-240,-240,-240},
    { DEF,-340,-690,-690,-690},
    { DEF,-330,-330,-330,-330}},
  { /* UG */
    { 0, 0, 0, 0, 0},
    { DEF,-720,-790,-960,-810},
    { DEF,-480,-480,-360,-480},
    { DEF,-660,-810,-920,-810},
    { DEF,-550,-440,-550,-360}},
  { /* AU */
    { 0, 0, 0, 0, 0},
    { DEF,-430,-600,-600,-600},
    { DEF,-260,-240,-240,-240},
    { DEF,-340,-690,-690,-690},
    { DEF,-330,-330,-330,-330}},
  { /* UA */
    { 0, 0, 0, 0, 0},
    { DEF,-400,-630,-890,-590},
    { DEF,-430,-510,-200,-180},
    { DEF,-380,-680,-890,-680},
    { DEF,-280,-140,-280,-140}},
  { /* @@ */
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF}}
};


//@@ mismatch_interior [in file] --> MMI --> mismatchI37 [in program]       
/* terminal mismatches */
/* mismatch free energies for interior loops at 37C */
PUBLIC int mismatchI37[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,    0,    0, -110,    0}, /* A@  AA  AC  AG  AU */
   {   0,    0,    0,    0,    0}, /* C@  CA  CC  CG  CU */
   {   0, -110,    0,    0,    0}, /* G@  GA  GC  GG  GU */
   {   0,    0,    0,    0,  -70}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,    0,    0, -110,    0}, /* A@  AA  AC  AG  AU */
   {   0,    0,    0,    0,    0}, /* C@  CA  CC  CG  CU */
   {   0, -110,    0,    0,    0}, /* G@  GA  GC  GG  GU */
   {   0,    0,    0,    0,  -70}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* UG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* AU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* UA */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   70,   70,  -40,   70}, /* A@  AA  AC  AG  AU */
   {   0,   70,   70,   70,   70}, /* C@  CA  CC  CG  CU */
   {   0,  -40,   70,   70,   70}, /* G@  GA  GC  GG  GU */
   {   0,   70,   70,   70,    0}},/* U@  UA  UC  UG  UU */
  { /* @@ */
   { 90, 90, 90, 90, 90},{ 90, 90, 90, 90,-20},{ 90, 90, 90, 90, 90},
   { 90,-20, 90, 90, 90},{ 90, 90, 90, 90, 20}}
};


//@@ mismatch_interior_enthalpies [in file] --> MMI_H --> mismatchIdH [in program]       
PUBLIC int mismatchIdH[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0, 0, 0, 0, 0},
    { DEF -1030,-950 -1030 -1030},
    { DEF,-520,-450,-520,-670},
    { DEF,-940,-940,-940,-940},
    { DEF,-810,-740,-810,-860}},
  { /* GC */
    { 0, 0, 0, 0, 0},
    { DEF,-520,-880,-560,-880},
    { DEF,-720,-310,-310,-390},
    { DEF,-710,-740,-620,-740},
    { DEF,-500,-500,-500,-570}},
  { /* GU */
    { 0, 0, 0, 0, 0},
    { DEF,-430,-600,-600,-600},
    { DEF,-260,-240,-240,-240},
    { DEF,-340,-690,-690,-690},
    { DEF,-330,-330,-330,-330}},
  { /* UG */
    { 0, 0, 0, 0, 0},
    { DEF,-720,-790,-960,-810},
    { DEF,-480,-480,-360,-480},
    { DEF,-660,-810,-920,-810},
    { DEF,-550,-440,-550,-360}},
  { /* AU */
    { 0, 0, 0, 0, 0},
    { DEF,-430,-600,-600,-600},
    { DEF,-260,-240,-240,-240},
    { DEF,-340,-690,-690,-690},
    { DEF,-330,-330,-330,-330}},
  { /* UA */
    { 0, 0, 0, 0, 0},
    { DEF,-400,-630,-890,-590},
    { DEF,-430,-510,-200,-180},
    { DEF,-380,-680,-890,-680},
    { DEF,-280,-140,-280,-140}},
  { /* @@ */
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF}}
};




//@@ mismatch_interior_1n [in file] --> MMI1N --> mismatch1nI37 [in program]       
PUBLIC int mismatch1nI37[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0, 0, 0, 0, 0},
    { 0, 0, 0,-110, 0},
    { 0, 0, 0, 0, 0},
    { 0,-110, 0, 0, 0},
    { 0, 0, 0, 0, -70}},
  { /* GC */
    { 0, 0, 0, 0, 0},
    { 0, 0, 0,-110, 0},
    { 0, 0, 0, 0, 0},
    { 0,-110, 0, 0, 0},
    { 0, 0, 0, 0, -70}},
  { /* GU */
    { 0, 0, 0, 0, 0},
    { 0,70,70, -40,70},
    { 0,70,70,70,70},
    { 0, -40,70,70,70},
    { 0,70,70,70, 0}},
  { /* UG */
    { 0, 0, 0, 0, 0},
    { 0,70,70, -40,70},
    { 0,70,70,70,70},
    { 0, -40,70,70,70},
    { 0,70,70,70, 0}},
  { /* AU */
    { 0, 0, 0, 0, 0},
    { 0,70,70, -40,70},
    { 0,70,70,70,70},
    { 0, -40,70,70,70},
    { 0,70,70,70, 0}},
  { /* UA */
    { 0, 0, 0, 0, 0},
    { 0,70,70, -40,70},
    { 0,70,70,70,70},
    { 0, -40,70,70,70},
    { 0,70,70,70, 0}},
  { /* @@ */
    {90,90,90,90,90},
    {90,90,90,90, -20},
    {90,90,90,90,90},
    {90, -20,90,90,90},
    {90,90,90,90,20}}
};

//@@ mismatch_interior_1n_enthalpies [in file] --> MMI1N_H --> mismatch1nIdH [in program]       
PUBLIC int mismatch1nIdH[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0, 0, 0, 0, 0},
    { DEF -1030,-950 -1030 -1030},
    { DEF,-520,-450,-520,-670},
    { DEF,-940,-940,-940,-940},
    { DEF,-810,-740,-810,-860}},
  { /* GC */
    { 0, 0, 0, 0, 0},
    { DEF,-520,-880,-560,-880},
    { DEF,-720,-310,-310,-390},
    { DEF,-710,-740,-620,-740},
    { DEF,-500,-500,-500,-570}},
  { /* GU */
    { 0, 0, 0, 0, 0},
    { DEF,-430,-600,-600,-600},
    { DEF,-260,-240,-240,-240},
    { DEF,-340,-690,-690,-690},
    { DEF,-330,-330,-330,-330}},
  { /* UG */
    { 0, 0, 0, 0, 0},
    { DEF,-720,-790,-960,-810},
    { DEF,-480,-480,-360,-480},
    { DEF,-660,-810,-920,-810},
    { DEF,-550,-440,-550,-360}},
  { /* AU */
    { 0, 0, 0, 0, 0},
    { DEF,-430,-600,-600,-600},
    { DEF,-260,-240,-240,-240},
    { DEF,-340,-690,-690,-690},
    { DEF,-330,-330,-330,-330}},
  { /* UA */
    { 0, 0, 0, 0, 0},
    { DEF,-400,-630,-890,-590},
    { DEF,-430,-510,-200,-180},
    { DEF,-380,-680,-890,-680},
    { DEF,-280,-140,-280,-140}},
  { /* @@ */
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF}}
};

//@@ mismatch_interior_23 [in file] --> MMI23 --> mismatch23I37 [in program]       
PUBLIC int mismatch23I37[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0, 0, 0, 0, 0},
    { 0, 0, 0,-110, 0},
    { 0, 0, 0, 0, 0},
    { 0,-110, 0, 0, 0},
    { 0, 0, 0, 0, -70}},
  { /* GC */
    { 0, 0, 0, 0, 0},
    { 0, 0, 0,-110, 0},
    { 0, 0, 0, 0, 0},
    { 0,-110, 0, 0, 0},
    { 0, 0, 0, 0, -70}},
  { /* GU */
    { 0, 0, 0, 0, 0},
    { 0,70,70, -40,70},
    { 0,70,70,70,70},
    { 0, -40,70,70,70},
    { 0,70,70,70, 0}},
  { /* UG */
    { 0, 0, 0, 0, 0},
    { 0,70,70, -40,70},
    { 0,70,70,70,70},
    { 0, -40,70,70,70},
    { 0,70,70,70, 0}},
  { /* AU */
    { 0, 0, 0, 0, 0},
    { 0,70,70, -40,70},
    { 0,70,70,70,70},
    { 0, -40,70,70,70},
    { 0,70,70,70, 0}},
  { /* UA */
    { 0, 0, 0, 0, 0},
    { 0,70,70, -40,70},
    { 0,70,70,70,70},
    { 0, -40,70,70,70},
    { 0,70,70,70, 0}},
  { /* @@ */
    {90,90,90,90,90},
    {90,90,90,90, -20},
    {90,90,90,90,90},
    {90, -20,90,90,90},
    {90,90,90,90,20}}
};

//@@ mismatch_interior_23_enthalpies [in file] --> MMI23_H --> mismatch23IdH [in program]       
PUBLIC int mismatch23IdH[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0, 0, 0, 0, 0},
    { DEF -1030,-950 -1030 -1030},
    { DEF,-520,-450,-520,-670},
    { DEF,-940,-940,-940,-940},
    { DEF,-810,-740,-810,-860}},
  { /* GC */
    { 0, 0, 0, 0, 0},
    { DEF,-520,-880,-560,-880},
    { DEF,-720,-310,-310,-390},
    { DEF,-710,-740,-620,-740},
    { DEF,-500,-500,-500,-570}},
  { /* GU */
    { 0, 0, 0, 0, 0},
    { DEF,-430,-600,-600,-600},
    { DEF,-260,-240,-240,-240},
    { DEF,-340,-690,-690,-690},
    { DEF,-330,-330,-330,-330}},
  { /* UG */
    { 0, 0, 0, 0, 0},
    { DEF,-720,-790,-960,-810},
    { DEF,-480,-480,-360,-480},
    { DEF,-660,-810,-920,-810},
    { DEF,-550,-440,-550,-360}},
  { /* AU */
    { 0, 0, 0, 0, 0},
    { DEF,-430,-600,-600,-600},
    { DEF,-260,-240,-240,-240},
    { DEF,-340,-690,-690,-690},
    { DEF,-330,-330,-330,-330}},
  { /* UA */
    { 0, 0, 0, 0, 0},
    { DEF,-400,-630,-890,-590},
    { DEF,-430,-510,-200,-180},
    { DEF,-380,-680,-890,-680},
    { DEF,-280,-140,-280,-140}},
  { /* @@ */
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF},
    { DEF, DEF, DEF, DEF, DEF}}
};

//@@ mismatch_multi [in file] --> MMM --> mismatchM37 [in program]       
PUBLIC int mismatchM37[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0,-110, -40,-130, -60},
    { DEF,-160, -90,-180,-110},
    { -30,-140, -70,-160, -90},
    { -20,-130, -60,-150, -80},
    { -10,-120, DEF,-140, -70}},
  { /* GC */
    { 0,-170, -80,-170,-120},
    { -20,-190,-100,-190,-140},
    { -30,-200,-110,-200,-150},
    { 0,-170, -80,-170,-120},
    { 0,-170, -80,-170,-120}},
  { /* GU */
    { 0, -70, -10, -70, -10},
    { -30,-100, -40,-100, -40},
    { -30,-100, -40,-100, -40},
    { -40,-110, DEF,-110, DEF},
    { -20, -90, -30, -90, -30}},
  { /* UG */
    { 0, -80, DEF, -80, -60},
    { -30,-110, -80,-110, -90},
    { -10, -90, -60, -90, -70},
    { -20,-100, -70,-100, -80},
    { -20,-100, -70,-100, -80}},
  { /* AU */
    { 0, -70, -10, -70, -10},
    { -30,-100, -40,-100, -40},
    { -30,-100, -40,-100, -40},
    { -40,-110, DEF,-110, DEF},
    { -20, -90, -30, -90, -30}},
  { /* UA */
    { 0, -80, DEF, -80, -60},
    { -30,-110, -80,-110, -90},
    { -10, -90, -60, -90, -70},
    { -20,-100, -70,-100, -80},
    { -20,-100, -70,-100, -80}},
  { /* @@ */
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0}}
};

//@@ mismatch_multi_enthalpies [in file] --> MMM_H --> mismatchMdH [in program]       
PUBLIC int mismatchMdH[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
    { 0,-740,-280,-640,-360},
    {-240,-980,-520,-880,-600},
    { 330,-410,50,-310, -30},
    {80,-660,-200,-560,-280},
    {-140,-880,-420,-780,-500}},
  { /* GC */
    { 0,-900,-410,-860,-750},
    {-160 -1060,-570 -1020,-910},
    {70,-830,-340,-790,-680},
    {-460 -1360,-870 -1320 -1210},
    { -40,-940,-450,-900,-790}},
  { /* GU */
    { 0,-740,-240,-720,-490},
    { 160,-580, -80,-560,-330},
    { 220,-520, -20,-500,-270},
    {70,-670,-170,-650,-420},
    { 310,-430,70,-410,-180}},
  { /* UG */
    { 0,-490, -90,-550,-230},
    {-150,-640,-240,-700,-380},
    { 510,20, 420, -40, 280},
    {10,-480, -80,-540,-220},
    { 100,-390,10,-450,-130}},
  { /* AU */
    { 0,-570, -70,-580,-220},
    { 160,-410,90,-420, -60},
    { 220,-350, 150,-360, 0},
    {70,-500, 0,-510,-150},
    { 310,-260, 240,-270,90}},
  { /* UA */
    { 0,-490, -90,-550,-230},
    { DEF,-540,-140,-600,-280},
    { 690, 200, 600, 140, 460},
    { -60,-550,-150,-610,-290},
    { -60,-550,-150,-610,-290}},
  { /* @@ */
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0}}
};




//@@ dangle5 ==> D5 ==> dangle5_37
/* 5' dangling ends (unpaird base stacks on first paired base) */
PUBLIC int dangle5_37[NBPAIRS+1][5]=
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF}, /* no pair */
    { INF, DEF, -30, -20, -10},
    { INF, -20, -30, 0, 0},
    { INF, -30, -30, -40, -20},
    { INF, -30, -10, -20, -20},
    { INF, -30, -30, -40, -20},
    { INF, -30, -10, -20, -20},
    { 0, 0, 0, 0, 0}
};

//@@ dangle5_enthalpies ==> D5_H ==> dangle5_dH
PUBLIC int dangle5_dH[NBPAIRS+1][5]=
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF}, /* no pair */
    { 0,-240, 330,80,-140},
    { 0,-160,70,-460, -40},
    { 0, 160, 220,70, 310},
    { 0,-150, 510,10, 100},
    { 0, 160, 220,70, 310},
    { 0, DEF, 690, -60, -60},
    { 0, 0, 0, 0, 0}
};
// PUBLIC int dangle5_H[NBPAIRS+1][5] =
// {/*   @     A     C     G     U   */
//    { INF,  INF,  INF,  INF,  INF},  /* no pair */
//    {   0, -240,  330,   80, -140},
//    {   0, -160,   70, -460,  -40},
//    {   0,  160,  220,   70,  310},
//    {   0, -150,  510,   10,  100},
//    {   0,  160,  220,   70,  310},
//    {   0,  -50,  690,  -60,  -60},
//    {   0,    0,    0,    0,   0}
// };


//@@ dangle3 ==> D3 ==> dangle3_37
/* 3' dangling ends (unpaired base stacks on second paired base */
PUBLIC int dangle3_37[NBPAIRS+1][5]=
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF},  /* no pair */
    { INF,-110, -40,-130, -60},
    { INF,-170, -80,-170,-120},
    { INF, -70, -10, -70, -10},
    { INF, -80, DEF, -80, -60},
    { INF, -70, -10, -70, -10},
    { INF, -80, DEF, -80, -60},
    { 0, 0, 0, 0, 0}
};

//@@ dangle3_enthalpies ==> D3_H ==> dangle3_dH
/* enthalpies for temperature scaling */
PUBLIC int dangle3_dH[NBPAIRS+1][5] =
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF},  /* no pair */
    { 0,-740,-280,-640,-360},
    { 0,-900,-410,-860,-750},
    { 0,-740,-240,-720,-490},
    { 0,-490, -90,-550,-230},
    { 0,-570, -70,-580,-220},
    { 0,-490, -90,-550,-230},
    { 0, 0, 0, 0, 0}
};


//@@ NINIO ==> NIN ==> ninio37, niniodH, MAX_NINIO
PUBLIC int ninio37=50;
PUBLIC int niniodH=0;
PUBLIC int MAX_NINIO=300;
//@@ "ML_params" ==> ML ==> ML_BASE37, ML_BASEdH, ML_closing37, ML_closingdH, ML_intern37, ML_interndH
PUBLIC int ML_BASE37=0;
PUBLIC int ML_BASEdH=0;
PUBLIC int ML_closing37=340;
PUBLIC int ML_closingdH=0;
PUBLIC int ML_intern37=40; 
PUBLIC int ML_interndH=0;
//@@ "Misc") ==> MISC ==> DuplexInit37, DuplexInitdH, TerminalAU37, TerminalAUdH, lxc37, 0
PUBLIC int DuplexInit37=410;
PUBLIC int DuplexInitdH=0;
PUBLIC int TerminalAU37=50;
PUBLIC int TerminalAUdH=0;
PUBLIC double lxc37=107.856;
//0

//Following 6 records do exist in both Turner1999 and Turner2004
PUBLIC int TripleC37=100;
PUBLIC int TripleCdH=1860;
PUBLIC int MultipleCA37=30;
PUBLIC int MultipleCAdH=340;
PUBLIC int MultipleCB37=160;
PUBLIC int MultipleCBdH=760;




//@@ "Triloops" ==> TRI ==> Triloops+c*6, Triloop37[c], TriloopdH[c]
//Triloops are empty in Turner1999
PUBLIC char Triloops[241] = "";
PUBLIC int Triloop37[40] = {   };
PUBLIC int TriloopdH[40] = {   };

//@@ "Tetraloops" ==> TL ==> Tetraloops+c*7, Tetraloop37[c], TetraloopdH[c]
// Tetraloops
PUBLIC char Tetraloops[281] =
    "GGGGAC "
    "GGUGAC "
    "CGAAAG "
    "GGAGAC "
    "CGCAAG "
    "GGAAAC "
    "CGGAAG "
    "CUUCGG "
    "CGUGAG "
    "CGAAGG "
    "CUACGG "
    "GGCAAC "
    "CGCGAG "
    "UGAGAG "
    "CGAGAG "
    "AGAAAU "
    "CGUAAG "
    "CUAACG "
    "UGAAAG "
    "GGAAGC "
    "GGGAAC "
    "UGAAAA "
    "AGCAAU "
    "AGUAAU "
    "CGGGAG "
    "AGUGAU "
    "GGCGAC "
    "GGGAGC "
    "GUGAAC "
    "UGGAAA "
;
PUBLIC int Tetraloop37[40] = {   
  	20,  	20,  	40,  	20,  	40,
	20,  	40,  	80,  	40,  	150,
	130,  	70,  	90,  	230,  	140,
	250,  	140,  	220,  	280,  	270,
	170,  	270,  	300,  	300,  	190,
	300,  	170,  	270,  	220,  	270
};
PUBLIC int TetraloopdH[40] = {   
  	-1110,	-1110,	-1340,	-1110,	-1340,	
	-1110,	-1340,	-1210,	-1340,	-1340,
	-1210,	-1110,	-1340,	-1060,	-1340,	 
	-740,	-1340,  -1140,	-1060,	-1020,	
	-1110,	-780,   -740,	-740,	-1340,	
	-740,	-1110,	-1020,  -900,	-780
};


//@@ "Hexaloops") ==> HEX ==> Hexaloops+c*9, Hexaloop37[c], HexaloopdH[c]
//Hexaloops are empty in Turner1999
PUBLIC char Hexaloops[361] = "";
PUBLIC int Hexaloop37[40] = {  };
PUBLIC int HexaloopdH[40] = {  };


//@@ "END") ==> QUIT
//@@  else return UNKNOWN;



#include "intloops.h"