#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <grace_np.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <signal.h>
#include <cash2003.h>
#include <cash2.h>
#include <mersenne.h>
#include <cash2-s.h>

static TYPE2** Medium;
static TYPE2** Cost;
static TYPE2** ColorMap;
static TYPE2** Diffusion_plane;

static TYPE2 empty={0,0,0,0,0,0.,0.,0.,0.,0.};



double growthrate=0.6; // Speed of growth, chance to invade empty space
double cost_s1 = 0.1;
double cost_s2 = 0.1; // Cost of antibiotic production
float death=0.07; // Death rate
float diffusion=0.2;
double decayrate = 0.03;
double rearrangement_rate_s1=0.015;
double rearrangement_rate_s2=0.015;
double growthrate_competitor = 0.6;
double production = 0.3;
int killing_time = 10000;
int number_sensitives = 1000;
int killing_period = 0;
double nastiness = 2.0;
int wt_seeds = 15;




char *name;
char buf[400]; //<-- temp buffer for string
char folder[400]; //<-- another temp buffer for string
char gridfolder[400]; //<-- another temp buffer for string
char costsfolder[400];
char productionfolder[400];
// char name[30] = "data";
void Initial(void)
{
  display = 0;
  char* readOut; // Stores arguments for loop
  for(int i = 0; i < (int)argc_g; i++) // Loops over the words on the command line
  {
    readOut = (char*)argv_g[i]; // Stores current word            // atof voor charater-to-float (of double)
    if(strcmp(readOut, "-n") == 0) nastiness = atof(argv_g[i+1]);         // atoi voor character-to-integer
    if(strcmp(readOut, "-d") == 0) diffusion = atof(argv_g[i+1]);
    if(strcmp(readOut, "-dec") == 0) decayrate = atof(argv_g[i+1]);
    if(strcmp(readOut, "-X") == 0) display = 0;
    if(strcmp(readOut, "-c1") == 0) cost_s1 = atof(argv_g[i+1]);
    if(strcmp(readOut, "-c2") == 0) cost_s2 = atof(argv_g[i+1]);
    if(strcmp(readOut, "-pr") == 0) production = atof(argv_g[i+1]);
    if(strcmp(readOut, "-name") == 0) name = argv_g[i+1];
    if(strcmp(readOut, "-rs") == 0) ulseedG = atoi(argv_g[i+1]);
    if(strcmp(readOut, "-gr") == 0) growthrate = atof(argv_g[i+1]);
    if(strcmp(readOut, "-gc") == 0) growthrate_competitor = atof(argv_g[i+1]);
    if(strcmp(readOut, "-ns") == 0) number_sensitives = atoi(argv_g[i+1]);
    if(strcmp(readOut, "-r1") == 0) rearrangement_rate_s1 = atof(argv_g[i+1]);
    if(strcmp(readOut, "-r2") == 0) rearrangement_rate_s2 = atof(argv_g[i+1]);
    if(strcmp(readOut, "-wt") == 0) wt_seeds = atoi(argv_g[i+1]);
  }

  sprintf(folder,"/linuxhome/tmp/bramve/%s/",name);
  sprintf(buf,"mkdir %s -p",folder);
  int success = system(buf);
  if(success == -1) {
    printf("Error in system call to make dir...\n"); exit(0);
  } //Exits if fails.

  init_genrand(ulseedG); // random seed moet weer geinitialiseerd worden volgens Bram

  sprintf(gridfolder,"/linuxhome/tmp/bramve/%s/grid",name);
  sprintf(buf,"mkdir %s -p",gridfolder);
  success = system(buf);
  if(success == -1) {
    printf("Error in system call to make grid dir...\n"); exit(0);
  } //Exits if fails.

  sprintf(productionfolder,"/linuxhome/tmp/bramve/%s/grid",name);
  sprintf(buf,"mkdir %s -p",gridfolder);
  success = system(buf);
  if(success == -1) {
    printf("Error in system call to make grid dir...\n"); exit(0);
  } //Exits if fails.

  sprintf(costsfolder,"/linuxhome/tmp/bramve/%s/grid",name);
  sprintf(buf,"mkdir %s -p", costsfolder);
  success = system(buf);
  if(success == -1){
    printf("Error in system call to make costs dir...\n"); exit(0);}

	MaxTime = 100000;
	nrow = 400;
	ncol = 400;
	nplane =4;
	scale = 1;
	boundary = WRAP;
	//ulseedG = 57;

	boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.};
}

void InitialPlane(void)
{


	MakePlane(&Medium, &Cost, &Diffusion_plane, &ColorMap);
		TYPE2 bg, s1, s2, s3, s4, s5;
		bg.val = 0;
		s1.val = 1;
		s2.val = 2;
    s3.val = 3;
    s4.val = 4;
    s5.val = 5;

    double chance = 0.0000625; // very scientific
    printf("%f", chance);
		InitialSetS(Medium, 3, bg, s1, chance, s3, chance, s4, chance);
    InitialSet(Diffusion_plane, 0, 0);
}

 void NextState(int row, int col)
 {
   TYPE2 neighbor = empty;
   neighbor = RandomMooreS8(Medium, row, col);

   if(Medium[row][col].val == 0 && neighbor.val == 1){
     if(genrand_real1() < growthrate){
       if(genrand_real1() > rearrangement_rate_s1){ // WT reproduction
         Medium[row][col] = neighbor;
      }

			else {
				Medium[row][col].val = 2;
      }
    }
  }

  if(Medium[row][col].val == 0 && neighbor.val == 4){
    if(genrand_real1() < growthrate){
      if(genrand_real1() > rearrangement_rate_s2){ // WT reproduction
        Medium[row][col] = neighbor;
     }

     else {
       Medium[row][col].val = 5;
     }
   }
 }
  // Mutant reproduction
  if(Medium[row][col].val == 0 && neighbor.val == 2){
  	if(genrand_real1() < (growthrate - cost_s1)){
  		Medium[row][col] = neighbor;
    }
  }


  if(Medium[row][col].val == 1 || Medium[row][col].val == 2 || Medium[row][col].val == 4 || Medium[row][col].val == 5){

  	if(genrand_real1() < death){
  		Medium[row][col] = empty;
    }


    if(Medium[row][col].val == 2 || Medium[row][col].val == 5){

  		Diffusion_plane[row][col].fval += production;
    }
  }


if(Medium[row][col].val == 0 && neighbor.val == 3){
  if(genrand_real1() < growthrate_competitor){
    Medium[row][col].val = 3;}}

if(Medium[row][col].val == 3){
  if(genrand_real1() < death + nastiness*Diffusion_plane[row][col].fval){
    Medium[row][col].val = 0;}
  }

if(Diffusion_plane[row][col].fval > 0.0){
	Diffusion_plane[row][col].fval = Diffusion_plane[row][col].fval *(1 - decayrate);
}

}

void Update(void)
{
  int i, j;
  DiffusionFVAL(Diffusion_plane,diffusion,1);

  //Killing
  if(Time > 0 && Time%killing_time == 0){
    PerfectMix(Medium);


    int counter1 = 0;
    int arr_x1[wt_seeds];
    int arr_y1[wt_seeds];
    int arr_type1[wt_seeds];
    while(counter1 < wt_seeds){
      int x = genrand_int(1, nrow);
      int y = genrand_int(1, nrow);
      if(Medium[x][y].val == 1 || Medium[x][y].val == 2 ){
        arr_x1[counter1] = x;
        arr_y1[counter1] = y;
        arr_type1[counter1] = Medium[x][y].val;
        counter1 +=1;
      }
    }

    int counter2 = 0;
    int arr_x2[wt_seeds];
    int arr_y2[wt_seeds];
    int arr_type2[wt_seeds];
    while(counter2 < wt_seeds){
      int x = genrand_int(1, nrow);
      int y = genrand_int(1, nrow);
      if(Medium[x][y].val == 4 || Medium[x][y].val == 5){
        arr_x2[counter2] = x;
        arr_y2[counter2] = y;
        arr_type2[counter2] = Medium[x][y].val;
        counter2 +=1;
      }
    }

    int i, j;
		for(i=1;i<=nrow;i++)
		for(j=1;j<=ncol;j++){
      Diffusion_plane[i][j].fval = 0;
      if(Medium[i][j].val == 5) Medium[i][j] = empty;
      if(Medium[i][j].val == 4) Medium[i][j] = empty;
      if(Medium[i][j].val == 3) Medium[i][j] = empty;
      if(Medium[i][j].val == 2) Medium[i][j] = empty;
      if(Medium[i][j].val == 1) Medium[i][j] = empty;
    }
    for(i=0;i<=wt_seeds-1;i++){
      if(arr_type1[i] == 1){
        Medium[arr_x1[i]][arr_y1[i]].val = 1;
      }
    }
    for(i=0;i<=wt_seeds-1;i++){
      if(arr_type2[i] == 4){
        Medium[arr_x2[i]][arr_y2[i]].val = 4;
      }
    }

    PerfectMix(Medium);
    int additor = 0;
    while(additor < number_sensitives){
      int i = genrand_int(1, nrow);
      int j = genrand_int(1, ncol);
      if(Medium[i][j].val == 0){
        Medium[i][j].val = 3;
        Medium[i][j].fval = 0;
        additor++;
         }
	     }
     }


  FILE *populationfile;
  if(Time%1000 == 0){
    char populationbuffer[400];
    sprintf(populationbuffer,"%s/popsize2.txt",folder);    //T=1000 "/linuxhome/tmp/brem/Hierzo/grid_at_time_1000.dat"
    populationfile = fopen(populationbuffer, "a");
    if(Time == 0){
      fprintf(populationfile, "Time\tPopulation_wildtype\tPopulation_mutants\tsensitives\tpop_WT2\tpop_mutant2");
      fprintf(populationfile, "\n");}
    int pop_WT1 = countGlobal(Medium, 1);
    int pop_mutant1 = countGlobal(Medium, 2);
    int pop_sensitive = countGlobal(Medium, 3);
    int pop_WT2 = countGlobal(Medium, 4);
    int pop_mutant2 = countGlobal(Medium, 5);
    fprintf(populationfile, "%d\t%d\t%d\t%d\t%d\t%d\n", Time, pop_WT1, pop_mutant1, pop_sensitive, pop_WT2, pop_mutant2);
    fclose(populationfile);
  }


  if(Time%killing_time == 0){
    killing_period += killing_time;
  }

    FILE *valfile;
    if(Time%2500 == 0){
      char valbuffer[400];
      sprintf(valbuffer,"%s/grid_type%d.txt", gridfolder, Time);
      valfile = fopen(valbuffer, "a");
      for(i=1;i<=nrow;i++)
      {
        for(j=1;j<=ncol;j++)
        {
          fprintf(valfile, "%d\t%d\t%d\n", Medium[i][j].val, i, j);
        }
      }
      fprintf(valfile, "\n");
      fclose(valfile);
    }

    FILE *runfile;
    char commandfile[400];
    sprintf(commandfile,"%s/command%s.txt",folder, name);
    runfile = fopen(commandfile, "w");
    for(int i = 0; i < (int)argc_g; i++){
      fprintf(runfile,"%s ", (char*)argv_g[i]);}
    fclose(runfile);

    if(Time%50 ==0)
    {

      Synchronous(1, Medium);
      Display(Medium);
    }

}
