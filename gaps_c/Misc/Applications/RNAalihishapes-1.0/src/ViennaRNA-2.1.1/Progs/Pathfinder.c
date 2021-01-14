#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "findpath.h"
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include <time.h>
#include <float.h>
#include <math.h>
//KEY
#include <unistd.h>
#include "pair_mat.h"
//#include "part_func.h"
//#include "PS_dot.h"
//#include "H/pair_mat.h"
#include "part_func.h"
#include "PS_dot.h"
#include "MEA.h"

#include "aln_util.h"
#include "alifold.h"


#define GRAPHSIZE 2048
#define INFINITY GRAPHSIZE*GRAPHSIZE
#define MAX(a, b) ((a > b) ? (a) : (b))
#define HISHREPLENGTH 1024
#define HISHREPSIZE 2048

//int e; /* The number of nonzero edges in the graph */
int n_node; /* The number of nodes in the graph */
float distances[GRAPHSIZE][GRAPHSIZE]; /* distances[i][j] is the distance between node i and j; or 0 if there is no direct connection */
float d[GRAPHSIZE]; /* d[i] is the length of the shortest path between the source (s) and node i */
int prev[GRAPHSIZE]; /* prev[i] is the node that comes right before i in the shortest path from the source to i*/
char seq[HISHREPLENGTH], struc1[HISHREPLENGTH], struc2[HISHREPLENGTH];
int anchors[100];


typedef char** stringArray;
stringArray MallocStringArray(size_t SizeOfOneString, size_t StringCount)
{
  char** t=malloc(StringCount*sizeof(char*));
  size_t i;
  for(i=0;i<StringCount;++i)
    t[i]=malloc(SizeOfOneString);
  return t;
}
stringArray strings;

void FreeStringArray(stringArray StringArray, size_t StringCount)
{
  size_t i;
  for(i=0;i<StringCount;++i)
    free(StringArray[i]);
  free(StringArray);
} 

// typedef struct move {
//     int i;  
//     int j;
//     int E;
//     int when;  
// } move_t;

// typedef struct move {
//   int i;  /* i,j>0 insert; i,j<0 delete */
//   int j;
//   int when;  /* 0 if still available, else resulting distance from start */
//   int E;
// } move_t;


int getBasePairDistance(char *ss1, char *ss2) {
    short *pairTable1, *pairTable2;
    move_t *move_tList;
    int i, length, distance=0;
    pairTable1 = make_pair_table(ss1);
    pairTable2 = make_pair_table(ss2);
    length = (int) strlen(ss1);
    move_tList = (move_t *) space(sizeof(move_t)*length); 

    for ( i=1; i<=length; i++ ) {
        if (pairTable1[i] != pairTable2[i]) {
            if ( i<pairTable1[i] ) {
	        move_tList[distance].i = -i;
	        move_tList[distance].j = -pairTable1[i];
	        move_tList[distance++].when = 0;
            }
            if ( i<pairTable2[i] ) {
	        move_tList[distance].i = i;
	        move_tList[distance].j = pairTable2[i];
	        move_tList[distance++].when = 0;
            }
       }
  }
  free(pairTable1);
  free(pairTable2);
  free(move_tList);
  return distance;
}

void printD() {
	int i;

	printf("Distances:\n");
	for (i = 0; i < n_node; i++)
		printf("%10d", i);
	printf("\n");
	for (i = 0; i < n_node; i++) {
		printf("%10f", d[i]);
	}
	printf("\n");
}

/*
 * Prints the shortest path from the source to dest.
 *
 * dijkstra(int) MUST be run at least once BEFORE
 * this is called
 */
void printPath(int dest) {
	if (prev[dest] != -1)
		printPath(prev[dest]);
	printf("%d ", dest);
}

float calculateDistance(int i, int j, int maxkeep)
{

      int k_max=0;
      float en_max=-INFINITY;
      char struc_max[HISHREPLENGTH];
      int dist;
      path_t *route;
      int k;

      route=get_path(seq, strings[i], strings[j], maxkeep);
      dist=getBasePairDistance(strings[i], strings[j]);

      for (k=0;k<dist+1;k++){
	  if (route[k].en > en_max)
	  {
	      k_max = k;
	      en_max = route[k].en;
	      strcpy(struc_max, route[k].s);
	  }
      }

      distances[i][j] = en_max;
      return distances[i][j];
}

int dijkstra(int s, int t, int maxkeep) {
	int i, k, mini;
	int visited[GRAPHSIZE];
	int dist;
	path_t *route;

	// initialize d[], prev[], visited[]
	// empty the three arrays
	for (i = 0; i < n_node; ++i) {
		d[i] = -INFINITY;
		prev[i] = -1; /* no path has yet been found to i */
        	visited[i] = 0; /* the i-th element has not yet been visited */
	}
	//d[s] = 0;



	d[s] = -INFINITY;
	prev[s] = -1;
	for (i = 0; i < n_node; i++) {
	    if (i==s)
	    {
		continue;
	    }
	    d[i] = calculateDistance(0, i, maxkeep);
            //d[i] = distances[0][i];
	    prev[i] = s;
	    //printf("Path to %f[0][%d]: ", distances[0][i], i);
	}
	visited[s] = 1;
	
	
	for (k = 0; k < n_node; ++k) {
	        if (k==s)
		{
		  continue;
		}
		mini = -1;
		for (i = 0; i < n_node; ++i)
			if (!visited[i] && ((mini == -1) || (d[i] < d[mini])))
				mini = i;
  
		visited[mini] = 1;
		//printf("%d\t",mini);
		
		
		for (i = 0; i < n_node; ++i)
		{
			if (i==s)
			{
			    continue;
			}
			
                        //float distance_mini_i = distances[mini][i];
			//if (distance_mini_i == -INFINITY)
			//{
			      float distance_mini_i = calculateDistance(mini, i, maxkeep);
			//}
			if (MAX(d[mini],distance_mini_i) < d[i])
			{
				d[i] = MAX(d[mini],distance_mini_i);
				prev[i] = mini;
			}
		}
		
		if (mini==t)
		{
			int anchor_index = 0;
			
			int x=t;
			while (prev[x] != -1)
			{
			    //printf("%d ", x);
			    //printf("%s ", strings[x]);
			    anchors[anchor_index] = x;
			    anchor_index++;
			    x = prev[x];
			}
			
			// add the start item
			anchors[anchor_index] = s;
			anchor_index++;
			return anchor_index;
		}

	}
	
}







//KEY
extern void  read_parameter_file(const char fname[]);
extern float energy_of_circ_struct(const char *seq, const char *structure);
#define MAX_NUM_NAMES    500 
static const char scale[] = "....,....1....,....2....,....3....,....4"
			    "....,....5....,....6....,....7....,....8";

static void /*@exits@*/ usage(void);
static char **annote(const char *structure, const char *AS[]);
static void print_pi(const pair_info pi, FILE *file);
static void print_aliout(char **AS, plist *pl, int n_seq, char * mfe, FILE *aliout);
static void mark_endgaps(char *seq, char egap);
static cpair *make_color_pinfo(char **sequences, plist *pl, int n_seq, bondT *mfe);

  
static void /*@exits@*/ usage(void);
//static char *seq;
static short *S, *S1;
//static int BP_dist;
static move_t *path;
static int path_fwd; /* 1: struc1->struc2, else struc2 -> struc1 */

int main(int argc, char *argv[]) {


  char *string;    // consensus sequence
  char *structure=NULL, *cstruc=NULL;
  char *structure2=NULL, *cstruc2=NULL;
  char  ffname[20], gfname[20], fname[13]="";
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   n_seq, length, sym;//i, r;
  int   endgaps=0, mis=0;
  double min_en, min_en2, real_en, real_en2, sfact=1.07;
  int   pf=0, noPS=0, istty;    // istty means if it uses typing input
  char     *AS[MAX_NUM_NAMES];          /* aligned sequences */
  char     *names[MAX_NUM_NAMES];       /* sequence names */
  FILE     *clust_file = stdin;
  int circ=0;
  int doAlnPS=0;
  int doColor=0;
  int n_back=0;
  int eval_energy = 0;
  int doMEA=0;
  double MEAgamma = 1.;
  
  do_backtrack = 1;
  string=NULL;
  dangles=2;
  oldAliEn=0;
  
  // NOTE: r is not integer ==> conflicts TODO
  //KEY char *seq; 
  char *s1, *s2;
  int E, maxkeep=10;
  int verbose=0, i;
  path_t *route, *r;

  
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-') {
      switch ( argv[i][1] )
	{
	//case 'm': if (strcmp(argv[i],"-m")==0)
	//    sscanf(argv[++i], "%d", &maxkeep);
	//  break;
	//case 'v':  verbose = !strcmp(argv[i],"-v");
	//  break;
	//case 'd': if (strcmp(argv[i],"-d")==0)
	//    sscanf(argv[++i], "%d", &dangles);
	//  break;
	  
	case 'T':  if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &temperature);
	  if (!r) usage();
	  break;
	//NOTE ## the p option is not implemented in Pathfinder
	//case 'p':  pf=1;
	//  if (argv[i][2]!='\0')
	//    (void) sscanf(argv[i]+2, "%d", &do_backtrack);
	//  break;
	case 'n':
	  if ( strcmp(argv[i], "-noGU")==0) noGU=1;
	  if ( strcmp(argv[i], "-noCloseGU")==0) no_closingGU=1;
	  if ( strcmp(argv[i], "-noLP")==0) noLonelyPairs=1;
	  if ( strcmp(argv[i], "-nsp") ==0) {
	    if (i==argc-1) usage();
	    ns_bases = argv[++i];
	  }
	  if ( strcmp(argv[i], "-nc")==0) {
	    r=sscanf(argv[++i], "%lf", &nc_fact);
	    if (!r) usage();
	  }
	  break;
	case 'o':
	  if ( strcmp(argv[i], "-old")==0) oldAliEn=1;
	  break;
	case 'm':
	  if ( strcmp(argv[i], "-mis")==0) mis=1;
	  else usage();
	  break;
	case '4':
	  tetra_loop=0;
	  break;
	case 'e':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%d", &energy_set);
	  if (!r) usage();
	  break;
	case 'E':
	  endgaps=1;
	  break;
	case 'C':
	  fold_constrained=1;
	  break;
	case 'S':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%lf", &sfact);
	  if (!r) usage();
	  break;
	case 'd': dangles=0;
	  if (argv[i][2]!='\0') {
	    r=sscanf(argv[i]+2, "%d", &dangles);
	    if (r!=1) usage();
	  }
	  break;
	case 'P':
	  if (i==argc-1) usage();
	  ParamFile = argv[++i];
	  break;
	case 'M':
	  if (strcmp(argv[i], "-MEA")==0) pf=doMEA=1;
	  if (i<argc-1)
	    if (isdigit(argv[i+1][0]) || argv[i+1][0]=='.') {
	      r=sscanf(argv[++i], "%lf", &MEAgamma);
	      if (r!=1) usage();
	    }
	  break;
	case 'c':
	  if ( strcmp(argv[i], "-cv")==0) {
	    r=sscanf(argv[++i], "%lf", &cv_fact);
	    if (!r) usage();
	  } else {
	    if (strcmp(argv[i], "-circ")==0)
	      circ=1;
	    else
	      if (strcmp(argv[i], "-color")==0)
		doColor=1;
	  }
	  break;
	case 'a':
	  if ( strcmp(argv[i], "-aln")==0) {
	    doAlnPS=1;
	  }
	  break;
	case 'R':
	  if (i==argc-1) usage();
	  RibosumFile = argv[++i];
	  ribo=1;
	  break;
	 case 'r':
	  RibosumFile = NULL;
	  ribo=1;
	  break;
	case 's':
	  if (argv[i][2]=='e') eval_energy = 1;
	  else if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r= sscanf(argv[++i], "%d", &n_back);
	  if (!r) usage();
	  do_backtrack=0;
	  pf=1;
	  init_rand();
	  break;
	default: usage();
	}
    }
    else { /* doesn't start with '-' should be filename */
      if (i!=argc-1) usage();
      clust_file = fopen(argv[i], "r");
      if (clust_file == NULL) {
	fprintf(stderr, "can't open %s\n", argv[i]);
	usage();
      }

    }
  }
  

  
  
  //seq = get_line(stdin);
  //s1 = get_line(stdin);
  //s2 = get_line(stdin);

  //E = find_saddle(seq, s1, s2, maxkeep);
  //printf("saddle_energy = %6.2f\n", E/100.);
  
  
  //KEY
  make_pair_matrix();

  if (circ && noLonelyPairs)
    fprintf(stderr,
	    "warning, depending on the origin of the circular sequence, "
	    "some structures may be missed when using -noLP\n"
	    "Try rotating your sequence a few times\n");
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL) {
    nonstandards = space(33);
    c=ns_bases;
    i=sym=0;
    if (*c=='-') {
      sym=1; c++;
    }
    while (*c!='\0') {
      if (*c!=',') {
	nonstandards[i++]=*c++;
	nonstandards[i++]=*c;
	if ((sym)&&(*c!=*(c-1))) {
	  nonstandards[i++]=*c;
	  nonstandards[i++]=*(c-1);
	}
      }
      c++;
    }
  }
  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
  if ((fold_constrained)&&(istty)) {
    printf("Input constraints using the following notation:\n");
    printf("| : paired with another base\n");
    printf(". : no constraint at all\n");
    printf("x : base must not pair\n");
    printf("< : base i is paired with a base j<i\n");
    printf("> : base i is paired with a base j>i\n");
    printf("matching brackets ( ): base i pairs base j\n");
  }
  if (fold_constrained) {
    if (istty) printf("%s\n", scale);
    cstruc = get_line(stdin);
    cstruc2 = get_line(stdin);
  }

  if (istty && (clust_file == stdin)) {
    printf("\nInput aligned sequences in clustalw format\n");
    if (!fold_constrained) printf("%s\n", scale);
  }

  n_seq = read_clustal(clust_file, AS, names);

  if (endgaps)
    for (i=0; i<n_seq; i++) mark_endgaps(AS[i], '~');
  if (clust_file != stdin) fclose(clust_file);
  if (n_seq==0)
    nrerror("no sequences found");

  length = (int) strlen(AS[0]);
  structure = (char *) space((unsigned) length+1);
  structure2 = (char *) space((unsigned) length+1);
  if (fold_constrained) {
    if (cstruc!=NULL)
      strncpy(structure, cstruc, length);
    else
      fprintf(stderr, "constraints for start sturcture missing\n");
    if (cstruc2!=NULL)
      strncpy(structure2, cstruc2, length);
    else
      fprintf(stderr, "constraints for target structure missing\n");
  }

  if (circ) {
    min_en = circalifold((const char **)AS, structure);
    min_en2 = circalifold((const char **)AS, structure2);
  }
  else
  {
    min_en = alifold(AS, structure);
    min_en2 = alifold(AS, structure2);
  }

  if (circ) {
    int i; double s=0, s2=0;
    extern int eos_debug;
    eos_debug=-1; /* shut off warnings about nonstandard pairs */
    for (i=0; AS[i]!=NULL; i++) {
      s +=energy_of_circ_struct(AS[i], structure);
      s2 +=energy_of_circ_struct(AS[i], structure2);
    }
    real_en = s/i;
    real_en2 = s2/i;
  } else {
    float CVen, CVen2;
    real_en = energy_of_alistruct(AS, structure, n_seq, &CVen);
    real_en2 = energy_of_alistruct(AS, structure2, n_seq, &CVen2);
    printf(" (%6.2f) \n", real_en );
    printf(" (%6.2f) \n", real_en2 );
  }

  string = (mis) ?
    consens_mis((const char **) AS) : consensus((const char **) AS);
  printf("%s\n%s", string, structure);    // string is consensus sequence
  
  if (istty) {
    printf("\n minimum free energy = %6.2f kcal/mol (%6.2f + %6.2f)\n",
	   min_en, real_en, min_en - real_en);
  } else {
    printf(" (%6.2f = %6.2f + %6.2f) \n", min_en, real_en, min_en-real_en );
  }
  
  printf("%s", structure2);
  if (istty) {
    printf("\n minimum free energy = %6.2f kcal/mol (%6.2f + %6.2f)\n",
	   min_en2, real_en2, min_en2 - real_en2);
  } else {
    printf(" (%6.2f = %6.2f + %6.2f) \n", min_en2, real_en2, min_en2-real_en2 );
  }
  
  // only draw alirna.ps  aln.ps for start structure
  if (fname[0]!='\0') {
    strcpy(ffname, fname);
    strcat(ffname, "_ss.ps");
    strcpy(gfname, fname);
    strcat(gfname, "_ss.g");
  } else {
    strcpy(ffname, "alirna.ps");
    strcpy(gfname, "alirna.g");
  }
  if (!noPS) {
    char **A;
    A = annote(structure, (const char**) AS);
    if (doColor)
      (void) PS_rna_plot_a(string, structure, ffname, A[0], A[1]);
    else
      (void) PS_rna_plot_a(string, structure, ffname, NULL, A[1]);
    free(A[0]); free(A[1]); free(A);
  }
  if (doAlnPS)
    PS_color_aln(structure, "aln.ps", (const char const **) AS, (const char const **) names);

  { /* free mfe arrays but preserve base_pair for PS_dot_plot */
    struct bond  *bp;
    bp = base_pair; base_pair = space(16);
    free_alifold_arrays();  /* free's base_pair */
    base_pair = bp;
  }

  
// //#### initializing the parameter required by BFS algorithm ####  
//   char filename[100], struc1[HISHREPLENGTH], struc2[HISHREPLENGTH];
//   int maxkeep;
//   double myTime;
//   int dist;
//   path_t *route;
// 
//   int n = 0;
// 
//   if (argc!=3){
//     printf("Usage: %s inputfile maxkeep\n", argv[0]);
//     printf("\n\tProgram will then find hipath using findPath\n");
//     printf("\talgorithm.\n");
//     exit(1);
//   }
//   
//   myTime=clock();
//   maxkeep=10;
//   strcpy(filename, argv[1]);
//   maxkeep=atoi(argv[2]);
//   
//   int i, j, k, s, t;
// 
//   
//   strings=MallocStringArray(HISHREPLENGTH,HISHREPSIZE);
//   char acBuffer[HISHREPLENGTH];
//   FILE *fin2 = fopen(filename, "r");
//   fgets returns NULL at the end of the file
// 
//   while (fgets( acBuffer, HISHREPLENGTH, fin2 ) != NULL) 
//   {
//       printf("Read from file: %s\n", acBuffer)
//       int len = strlen(acBuffer) -1;
//       if( acBuffer[len] == '\n')
// 	  acBuffer[len] = '\0';
//       strcpy(strings[n], acBuffer);
//       n++;
//   }
//   
//   strcpy(seq, strings[n-1]);
//   strcpy(struc1, strings[0]);
//   strcpy(struc2, strings[1]);
//   n_node = n-1;
//   
//   /*
//   for (i = 0; i < n_node; ++i)
//   {
// 	  for (j = 0; j < n_node; ++j)
// 	  {
// 		distances[i][j] = -INFINITY;
// 	  }
//   }*/
// 
//   s = 0;
//   t = 1;
//   printf("s=%d, t=%d\n", s, t);
//   
// 	
//   int anchor_no = dijkstra(s,t,maxkeep);
//   
//   /*
//   printf("Hipath: ");
//   for (i = anchor_no-1; i >= 0; i--) {
//       printf("%d\t", anchors[i]);
//   }
//   printf("\n");*/
//   
//   for (i = anchor_no-1; i >= 1; i--) {
//       int k_max=0;
//       float en_max=-INFINITY;
//       char struc_max[HISHREPLENGTH];
//       route=get_path(seq, strings[anchors[i]], strings[anchors[i-1]], maxkeep);
//       dist=getBasePairDistance(strings[anchors[i]], strings[anchors[i-1]]);
//       for (k=0;k<dist+1;k++){
// 	  printf("%d: %g %s\n", k, route[k].en, route[k].s);
//       }
//   }
// 
// 
//   myTime=clock()-myTime;
//   FreeStringArray(strings,HISHREPSIZE);



//#### free resources ####
//   if (verbose) {
//     if (path_fwd)
//       print_path(seq,s1);
//     else
//       print_path(seq,s2);
//     free(path);
//     route = get_path(seq, s1, s2, maxkeep);
//     for (r=route; r->s; r++) {
//       printf("%s %6.2f - %6.2f\n", r->s, energy_of_struct(seq,r->s), r->en);
//       free(r->s);
//     }
//   }

  free(route);
  //free(seq);  no seq exists any more
  free(s1); free(s2);

  //KEY
   if (cstruc!=NULL) free(cstruc);
   if (cstruc2!=NULL) free(cstruc2);
   free(base_pair);
  (void) fflush(stdout);
  free(string);
  free(structure);
  free(structure2);
  for (i=0; AS[i]; i++) {
    free(AS[i]); free(names[i]);
  }
  
  return(EXIT_SUCCESS);
}













// int main(int argc, char *argv[]) {
//         char filename[100], struc1[HISHREPLENGTH], struc2[HISHREPLENGTH];
// 	int maxkeep;
// 	double myTime;
// 	int dist;
// 	path_t *route;
// 
// 	int n = 0;
// 
// 	if (argc!=3){
// 	  printf("Usage: %s inputfile maxkeep\n", argv[0]);
// 	  printf("\n\tProgram will then find hipath using findPath\n");
// 	  printf("\talgorithm.\n");
// 	  exit(1);
// 	}
// 	
// 	myTime=clock();
// 	maxkeep=10;
// 	strcpy(filename, argv[1]);
//         maxkeep=atoi(argv[2]);
// 	
// 	int i, j, k, s, t;
// 
// 	
// 	strings=MallocStringArray(HISHREPLENGTH,HISHREPSIZE);
// 	char acBuffer[HISHREPLENGTH];
// 	FILE *fin2 = fopen(filename, "r");
// 	// fgets returns NULL at the end of the file
// 
// 	while (fgets( acBuffer, HISHREPLENGTH, fin2 ) != NULL) 
// 	{
// 	    //printf("Read from file: %s\n", acBuffer)
// 	    int len = strlen(acBuffer) -1;
// 	    if( acBuffer[len] == '\n')
// 		acBuffer[len] = '\0';
// 	    strcpy(strings[n], acBuffer);
// 	    n++;
// 	}
// 	
// 	strcpy(seq, strings[n-1]);
// 	strcpy(struc1, strings[0]);
// 	strcpy(struc2, strings[1]);
// 	n_node = n-1;
// 	
// 	/*
// 	for (i = 0; i < n_node; ++i)
// 	{
// 		for (j = 0; j < n_node; ++j)
// 		{
// 		      distances[i][j] = -INFINITY;
// 		}
// 	}*/
// 
// 	s = 0;
// 	t = 1;
// 	//printf("s=%d, t=%d\n", s, t);
// 	
// 	      
// 	int anchor_no = dijkstra(s,t,maxkeep);
// 	
//         /*
// 	printf("Hipath: ");
// 	for (i = anchor_no-1; i >= 0; i--) {
// 	    printf("%d\t", anchors[i]);
// 	}
// 	printf("\n");*/
// 	
// 	for (i = anchor_no-1; i >= 1; i--) {
// 	    int k_max=0;
// 	    float en_max=-INFINITY;
// 	    char struc_max[HISHREPLENGTH];
// 	    route=get_path(seq, strings[anchors[i]], strings[anchors[i-1]], maxkeep);
// 	    dist=getBasePairDistance(strings[anchors[i]], strings[anchors[i-1]]);
// 	    for (k=0;k<dist+1;k++){
// 		printf("%d: %g %s\n", k, route[k].en, route[k].s);
// 	    }
// 	}
// 
// 	// TODO: 1. 3 different outputs into 3 different files
// 	//       2. output saddle structure
// 	//       3. change the ruby script main_b.rb to fit the change of this file
// 	//       4. consider whether calculation of saddle ss
// 	
// 	//       5. debug the 2 scripts, change main_b.rb and don't need folding.rb any more
// 
// 
// 	myTime=clock()-myTime;
// 	FreeStringArray(strings,HISHREPSIZE);
// 	return 0;
// }



//############################ private methods #########################
static void mark_endgaps(char *seq, char egap) {
  int i,n;
  n = strlen(seq);
  for (i=0; i<n && (seq[i]=='-'); i++) {
    seq[i] = egap;
  }
  for (i=n-1; i>0 && (seq[i]=='-'); i--) {
    seq[i] = egap;
  }
}

static void print_pi(const pair_info pi, FILE *file) {
  const char *pname[8] = {"","CG","GC","GU","UG","AU","UA", "--"};
  int i;

  /* numbering starts with 1 in output */
  fprintf(file, "%5d %5d %2d %5.1f%% %7.3f",
	  pi.i, pi.j, pi.bp[0], 100.*pi.p, pi.ent);
  for (i=1; i<=7; i++)
    if (pi.bp[i]) fprintf(file, " %s:%-4d", pname[i], pi.bp[i]);
  if (!pi.comp) fprintf(file, " +");
  fprintf(file, "\n");
}

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

/*-------------------------------------------------------------------------*/

static char **annote(const char *structure, const char *AS[]) {
  /* produce annotation for colored drawings from PS_rna_plot_a() */
  char *ps, *colorps, **A;
  int i, n, s, pairings, maxl;
  short *ptable;
  char * colorMatrix[6][3] = {
    {"0.0 1", "0.0 0.6",  "0.0 0.2"},  /* red    */
    {"0.16 1","0.16 0.6", "0.16 0.2"}, /* ochre  */
    {"0.32 1","0.32 0.6", "0.32 0.2"}, /* turquoise */
    {"0.48 1","0.48 0.6", "0.48 0.2"}, /* green  */
    {"0.65 1","0.65 0.6", "0.65 0.2"}, /* blue   */
    {"0.81 1","0.81 0.6", "0.81 0.2"}  /* violet */
  };

  n = strlen(AS[0]);
  maxl = 1024;

  A = (char **) space(sizeof(char *)*2);
  ps = (char *) space(maxl);
  colorps = (char *) space(maxl);
  ptable = make_pair_table(structure);
  for (i=1; i<=n; i++) {
    char pps[64], ci='\0', cj='\0';
    int j, type, pfreq[8] = {0,0,0,0,0,0,0,0}, vi=0, vj=0;
    if ((j=ptable[i])<i) continue;
    for (s=0; AS[s]!=NULL; s++) {
      type = pair[encode_char(AS[s][i-1])][encode_char(AS[s][j-1])];
      pfreq[type]++;
      if (type) {
	if (AS[s][i-1] != ci) { ci = AS[s][i-1]; vi++;}
	if (AS[s][j-1] != cj) { cj = AS[s][j-1]; vj++;}
      }
    }
    for (pairings=0,s=1; s<=7; s++) {
      if (pfreq[s]) pairings++;
    }

    if ((maxl - strlen(ps) < 192) || ((maxl - strlen(colorps)) < 64)) {
      maxl *= 2;
      ps = realloc(ps, maxl);
      colorps = realloc(colorps, maxl);
      if ((ps==NULL) || (colorps == NULL))
	  nrerror("out of memory in realloc");
    }

    if (pfreq[0]<=2 && pairings>0) {
      snprintf(pps, 64, "%d %d %s colorpair\n",
	       i,j, colorMatrix[pairings-1][pfreq[0]]);
      strcat(colorps, pps);
    }

    if (pfreq[0]>0) {
      snprintf(pps, 64, "%d %d %d gmark\n", i, j, pfreq[0]);
      strcat(ps, pps);
    }
    if (vi>1) {
      snprintf(pps, 64, "%d cmark\n", i);
      strcat(ps, pps);
    }
    if (vj>1) {
      snprintf(pps, 64, "%d cmark\n", j);
      strcat(ps, pps);
    }
  }
  free(ptable);
  A[0]=colorps;
  A[1]=ps;
  return A;
}

/*-------------------------------------------------------------------------*/

#define PMIN 0.0008
static int compare_pair_info(const void *pi1, const void *pi2) {
  pair_info *p1, *p2;
  int  i, nc1, nc2;
  p1 = (pair_info *)pi1;  p2 = (pair_info *)pi2;
  for (nc1=nc2=0, i=1; i<=6; i++) {
    if (p1->bp[i]>0) nc1++;
    if (p2->bp[i]>0) nc2++;
  }
  /* sort mostly by probability, add
     epsilon * comp_mutations/(non-compatible+1) to break ties */
  return (p1->p + 0.01*nc1/(p1->bp[0]+1.)) <
	 (p2->p + 0.01*nc2/(p2->bp[0]+1.)) ? 1 : -1;
}

static void print_aliout(char **AS, plist *pl, int n_seq, char * mfe, FILE *aliout) {
  int i, j, k, n, num_p=0, max_p = 64;
  pair_info *pi;
  double *duck, p;
  short *ptable;
  for (n=0; pl[n].i>0; n++);

  max_p = 64; pi = space(max_p*sizeof(pair_info));
  duck =  (double *) space((strlen(mfe)+1)*sizeof(double));
  ptable = make_pair_table(mfe);

  for (k=0; k<n; k++) {
    int s, type;
    p = pl[k].p; i=pl[k].i; j = pl[k].j;
    duck[i] -=  p * log(p);
    duck[j] -=  p * log(p);

    if (p<PMIN) continue;

    pi[num_p].i = i;
    pi[num_p].j = j;
    pi[num_p].p = p;
    pi[num_p].ent =  duck[i]+duck[j]-p*log(p);
    for (type=0; type<8; type++) pi[num_p].bp[type]=0;
    for (s=0; s<n_seq; s++) {
      int a,b;
      a=encode_char(toupper(AS[s][i-1]));
      b=encode_char(toupper(AS[s][j-1]));
      type = pair[a][b];
      if ((AS[s][i-1] == '-')||(AS[s][j-1] == '-')) type = 7;
      if ((AS[s][i-1] == '~')||(AS[s][j-1] == '~')) type = 7;
      pi[num_p].bp[type]++;
      pi[num_p].comp = (ptable[i] == j) ? 1:0;
    }
    num_p++;
    if (num_p>=max_p) {
      max_p *= 2;
      pi = xrealloc(pi, max_p * sizeof(pair_info));
    }
  }
  free(duck);
  pi[num_p].i=0;
  qsort(pi, num_p, sizeof(pair_info), compare_pair_info );

  /* print it */
  fprintf(aliout, "%d sequence; length of alignment %d\n",
	  n_seq, (int) strlen(AS[0]));
  fprintf(aliout, "alifold output\n");

  for (k=0; pi[k].i>0; k++) {
    pi[k].comp = (ptable[pi[k].i] == pi[k].j) ? 1:0;
    print_pi(pi[k], aliout);
  }
  fprintf(aliout, "%s\n", mfe);
  free(ptable);
  free(pi);
}


static cpair *make_color_pinfo(char **sequences, plist *pl, int n_seq, bondT *mfe) {
  /* produce info for PS_color_dot_plot */
  cpair *cp;
  int i, n,s, a, b,z,t,j, c;
  int pfreq[7];
  for (n=0; pl[n].i>0; n++);
  c=0;
  cp = (cpair *) space(sizeof(cpair)*(n+1));
  for (i=0; i<n; i++) {
    int ncomp=0;
    if(pl[i].p>PMIN) {
      cp[c].i = pl[i].i;
      cp[c].j = pl[i].j;
      cp[c].p = pl[i].p;
      for (z=0; z<7; z++) pfreq[z]=0;
      for (s=0; s<n_seq; s++) {
	a=encode_char(toupper(sequences[s][cp[c].i-1]));
	b=encode_char(toupper(sequences[s][cp[c].j-1]));
	if ((sequences[s][cp[c].j-1]=='~')||(sequences[s][cp[c].i-1] == '~')) continue;
	pfreq[pair[a][b]]++;
      }
      for (z=1; z<7; z++) {
	if (pfreq[z]>0) {
	  ncomp++;
	}}
      cp[c].hue = (ncomp-1.0)/6.2;   /* hue<6/6.9 (hue=1 ==  hue=0) */
      cp[c].sat = 1 - MIN2( 1.0, (float) (pfreq[0]*2./*pi[i].bp[0]*//(n_seq)));
      c++;
    }
  }
  for (t=1; t<=mfe[0].i; t++) {
    int nofound=1;
      for (j=0; j<c; j++) {
	if ((cp[j].i==mfe[t].i)&&(cp[j].j==mfe[t].j)) {
	  cp[j].mfe=1;
	  nofound=0;
	  break;
	}
      }
      if(nofound) {
	fprintf(stderr,"mfe base pair with very low prob in pf: %d %d\n",mfe[t].i,mfe[t].j);
	cp = (cpair *) realloc(cp,sizeof(cpair)*(c+1));
	cp[c].i = mfe[t].i;
	cp[c].j = mfe[t].j;
	cp[c].p = 0.;
	cp[c].mfe=1;
	c++;
      }
    }
  return cp;
}

/*-------------------------------------------------------------------------*/

static void usage(void)
{
  nrerror("usage:\n"
	  "RNAalifold [-cv float] [-nc float] [-E] [-old] [-r] [-R ribosum]\n"
	  "        [-mis] [-aln] [-color] [-circ] [-s num] [-se num]\n"
	  "        [-p[0]] [-C] [-T temp] [-4] [-d] [-noGU] [-noCloseGU]\n"
	  "        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale]"
	  );
}
