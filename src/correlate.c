/*
 * correlate.c    1.0  2001/05/16
 *
 * Copyright (c) 2001 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke and J. Daniel Gezelter, "A Random Sequential 
 *    Adsorption model for the differential coverage of Gold (111)
 *    surfaces by two related Silicon phthalocyanines,"
 *    J. Phys. Chem. B 105, pp. 6515-6519 (2001). 
 *    DOI: 10.1021/jp010985m
 * 
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. ALL EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND
 * WARRANTIES, INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE OR NON-INFRINGEMENT, ARE HEREBY
 * EXCLUDED.  THE UNIVERSITY OF NOTRE DAME AND ITS LICENSORS SHALL NOT
 * BE LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THE SOFTWARE OR ITS
 * DERIVATIVES. IN NO EVENT WILL THE UNIVERSITY OF NOTRE DAME OR ITS
 * LICENSORS BE LIABLE FOR ANY LOST REVENUE, PROFIT OR DATA, OR FOR
 * DIRECT, INDIRECT, SPECIAL, CONSEQUENTIAL, INCIDENTAL OR PUNITIVE
 * DAMAGES, HOWEVER CAUSED AND REGARDLESS OF THE THEORY OF LIABILITY,
 * ARISING OUT OF THE USE OF OR INABILITY TO USE SOFTWARE, EVEN IF THE
 * UNIVERSITY OF NOTRE DAME HAS BEEN ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGES.
 *
 * This software is not designed or intended for use in on-line
 * control of aircraft, air traffic, aircraft navigation or aircraft
 * communications; or in the design, construction, operation or
 * maintenance of any nuclear facility. Licensee represents and
 * warrants that it will not use or redistribute the Software for such
 * purposes.  
 */

/***************************************************************************
 * This program is used to read the output of the bent_umbrella and
 * octopus simulations and compute the radial distribution of the
 * objects on the surface, the radial correlation of the umbrellas,
 * and the orientation correlation parameter.
 *
 * Usage:
 *
 *  correlate <orientation file1> <orientation file2> ... 
 * 
 * NOTE: You may want to only do one at a time, because these guys can
 * take a while. However if multiple files are given, the final output
 * file will be an average of all the correlations.
 ****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parameters.h"
#include "circle_intersect.h"

struct linked_hist{ /* a linked list of histograms */
  struct linked_hist *next;
  double hist[max_bin];
};
struct linked_hist *first_corr_g = NULL; /* the first hist in the
					    correlated g(r) linked list*/
struct linked_hist *first_g = NULL; /* the first hist in the normal
				       g(r) linked hist. */



int main(argc, argv)
     int argc;
     char *argv[];
{
  double map_x(double); /*function for mapping the periodic box. */
  double map_y(double); /*function for mapping the periodic box. */

  int n_atoms; /*the number of particles, umbrellas, atoms, whatever.*/
  struct coords *locations; /*an array to hold all of the landed locations.*/
  const double hist_length = (NND * MAXLATTICE) / 2.0; /*the max hist_length */
  double dx; /* temp variable for the change in x */
  double dy; /* temp variable for the change in y */
  double rjk_sqr; /*the distance between the two points squared. */
  double rjk; /*the distance between the two points. */
  double r_lower; /* The inner shell radius when calculating the ideal 
		    density in the g(r) calculation.*/
  double nideal; /*the ideal number of particles in the bin */
  double r_upper; /* The outer shell radius when calculating the ideal
		     density in the g(r) calculation. */
  int bin; /*keeps track of which bin in the histogram the object falls in. */
  const double delr = hist_length / max_bin; /*the width of each bin. */
  struct linked_hist *g_of_r; /*the g(r) */
  struct linked_hist *corr_g_of_r; /* the correlated g(r) */
  struct linked_hist *current_g; /*used in calculating the total g(r) */
  struct linked_hist *current_corr_g; /*used in calc'ing the total corr_g(r) */
  double tot_g_of_r[max_bin]; /*the total g(r) */
  double tot_corr_g_of_r[max_bin]; /*the total correlated g(r) */
  double correlation = 0.0; /*the correlation parameter. */
  double total_correlation = 0.0; /*the correlation of all files. */
  double *corr_hist; /*histogram of the correlations. */
  double corr; /*the correlation of the two particles.*/
  int *hist; /*the histogram of the radial distribution. */
  double area; /* temporary variable to store the area of the surface. */
  double area_const; /*a constant used in the calculation of the bin area*/
  double dens; /* the average density of particles on the surface. */
  int i,j,k; /*these are my loop counters. */
  int excess; /*keeps track of the points which lie outside our measurement. */
  
  FILE *in_file; /* the current input file */
  FILE *out_file; /* the output file */

  char *out_name = "correlation.dat"; /*the output filename */
  char *in_name; /* the current input filename */
  char read_line[256]; /* read line for the file input. */
  

  
  /*check for proper usage */
  
  if(argc < 2){
    (void)fprintf(stderr,
		  "Proper usage:\n"
		  "\n"
		  "     %s filename1 filename2 ...\n\n",
		  argv[0]);
    exit(8);
  }
 
  /* allocate memory off the heap for the radial arrays. */
  
  hist = (int *)calloc(max_bin, sizeof(int));
  corr_hist = (double *)calloc(max_bin, sizeof(double));


  /* loop over the files and calculate the g(r)'s*/
  
  argc--; /* so we don't count the program name as a file. */

  for(i = 0; i < argc; i++){
    
    /* read the next file */
    
    in_name = argv[i+1];
    in_file = fopen(in_name, "r");
    if(in_file == NULL){
      (void)fprintf(stderr,"Could not open %s\n", in_name);
      exit(8);
    }
    

    /*read the first line and toss it, then read the number of
      particles and the area. */
    
    (void)fgets(read_line, sizeof(read_line), in_file);

    (void)fgets(read_line, sizeof(read_line), in_file);

    (void)sscanf(read_line, "%d\t%lf", &n_atoms, &area);    
  

    /* allocate memory for the locations array and the linked lists. */
    
    g_of_r = (struct linked_hist *)malloc(sizeof(struct linked_hist));
    corr_g_of_r = (struct linked_hist *)malloc(sizeof(struct linked_hist));
    locations = (struct coords *)calloc(n_atoms, sizeof(struct coords));


    /* read in the locations and orientations */
    
    (void)fgets(read_line, sizeof(read_line), in_file); /* toss this line */
    
    for(j = 0; j < n_atoms; j++){
      (void)fgets(read_line, sizeof(read_line), in_file);
      (void)sscanf(read_line, "%lf\t%lf\t%lf\n", &locations[j].x, 
		   &locations[j].y, &locations[j].theta);
    }
    
    /* close the input file. */
    
    (void)fclose(in_file);

    /* initialize all the histograms */
    
    for(j = 0; j < max_bin; j++){
      g_of_r->hist[j] = 0.0;
      corr_g_of_r->hist[j] = 0.0;
      hist[j] = 0;
      corr_hist[j] = 0.0;
    }
    excess = 0;



    /* loop over all pairs. */
    
    for(j = 0; j < n_atoms - 1; j++){
      
      for(k = j+1; k < n_atoms; k++){
	
	/* find the distance between the pair. */
	
	dx = locations[j].x - locations[k].x;
	dy = locations[j].y - locations[k].y;
	
	/* map the distances into the periodic box. */
	
	dx = map_x(dx);
	dy = map_y(dy);
	
	/*distance formula */

	rjk_sqr = dx*dx + dy*dy;
	rjk = sqrt(rjk_sqr);
	
	/* place the location in a bin */
	
	bin = (int)(rjk / delr);
	
	/* calculate the correlation term for the bin */
	
	corr = cos(locations[j].theta - locations[k].theta);
	correlation += 2.0 * corr; /*factor of 2, because each can see the 
				     other. */
	
	/* place the values in their respective histograms. */
	
	if(bin < max_bin){ /* make sure we don't get an out of bounds index*/
	  hist[bin] += 2; /*factor of 2, because each can see the other. */
	  corr_hist[bin] += 2.0 * corr;
	}
	else{
	  excess += 2; /* keep track of how many were outside our measurement*/
	} 
      }
    }
    
    /* calculate the density */
    
    dens = ((double)n_atoms) / area;
    area_const = dens * M_PI;
    

    /* calculate the g(r)'s */
    
    for(bin = 0; bin < max_bin; bin++){
      r_lower = ((double)bin) * delr;
      r_upper = r_lower + delr;
      nideal = area_const * (r_upper * r_upper - r_lower * r_lower);
      g_of_r->hist[bin] = ((double)hist[bin]) / ((double)n_atoms*nideal);
      corr_g_of_r->hist[bin] = corr_hist[bin] / ((double)n_atoms*nideal);

      /* the next step separates the radial distribution function out
	 of the angular distribution function. */

      if(g_of_r->hist[bin] != 0){
	corr_g_of_r->hist[bin] = corr_g_of_r->hist[bin] / g_of_r->hist[bin];
      }
      else{
	corr_g_of_r->hist[bin] = 0.0;
      }
    }
    
    /*calculate the correlation parameter */
    
    correlation /= ((double)n_atoms)*((double)(n_atoms - 1));
    total_correlation += correlation;


    /* take care of the linked lists. */

    g_of_r->next = first_g;
    corr_g_of_r->next = first_corr_g;
    first_g = g_of_r;
    first_corr_g = corr_g_of_r;
    
  }


  /* calculate the total correlation. */

  total_correlation /= (double)argc;

  /*calculate the total g(r)'s */

  for(i = 0; i < max_bin; i++){
    tot_g_of_r[i] = 0.0;
    tot_corr_g_of_r[i] = 0.0;
  }
  
  /*loop over all bins and collect the averages. */
  
  for(bin = 0; bin < max_bin; bin++){
    current_g = first_g;
    current_corr_g = first_corr_g;
    
    while(current_g != NULL){
      tot_g_of_r[bin] += current_g->hist[bin];
      tot_corr_g_of_r[bin] += current_corr_g->hist[bin];
      current_g = current_g->next;
      current_corr_g = current_corr_g->next;
    }
    tot_g_of_r[bin] /= (double)argc;
    tot_corr_g_of_r[bin] /= (double)argc;
  }

  /* write the output */
  
  out_file = fopen(out_name, "w");
  (void)fprintf(out_file, "#correlation parameter = %lf\n",
		total_correlation);
  (void)fprintf(out_file, "#radius\t g(r)\t corr g(r)\n");
  
  for(i = 0; i < max_bin; i++){
    (void)fprintf(out_file, "%lf\t%lf\t%lf\n", (delr * (double)(i + 1)),
		  tot_g_of_r[i], tot_corr_g_of_r[i]);
  }
  
  (void)fclose(out_file);
  return(0);
}


/* function maps the x direction of the periodic box */

double map_x(unmapped)
     double unmapped;
{ 
  const double length_X = NND * MAXLATTICE;/* used by the mapper */
  double mapped;

  if(unmapped < 0){
    mapped = unmapped - length_X * (double)((int)((unmapped / length_X)
						  - 0.5));
  }
  else{
    mapped = unmapped - length_X * (double)((int)((unmapped / length_X)
						  + 0.5));
  }
  return(mapped);
}


/* responsible for mapping the y direction in the periodic box. */

double map_y(unmapped)
     double unmapped;
{
  const double sqrt3 = 1.732050808; /* square root of three, a commonly
				       occurring constant in an hcp
				       lattice.  (the underlying lattice in 
				       the simulation. */
  
  const double length_Y = sqrt3 * NND * MAXLATTICE;/*used by the periodic
						     mapper . */
  double mapped;
  
  if(unmapped < 0){
    mapped = unmapped - length_Y * (double)((int)((unmapped / length_Y) 
						  - 0.5));
  }
  else{
    mapped = unmapped - length_Y * (double)((int)((unmapped / length_Y)
						  + 0.5));
  }
  return(mapped);
}
