/*
 * octopus.c    1.0  2001/05/16
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
 *    J. Phys. Chem. A XX, XXXX-YYYY (2001)). 
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


/****************************************************************************
 * This program is based on the bent_umbrella code. Modifications have
 * been made such that the object is now an "octopus". Namely, a
 * circle that lies parallel to the attachment plane and can attach
 * via a specified number of legs distributed radially around the
 * perimeter of the circle. Here the test for overlap is a simple
 * distance formula.
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "parameters.h"
#include "circle_intersect.h"




int debug = 0; /*boolean for the debugging output option.*/
int lattice = 0; /*boolean for turning on the underlying lattice. */
char *program_name; /*the name of the program in case we need it*/

int main(argc, argv)
     int argc;
     char *argv[];
{
  double map_x(double); /*function for mapping the periodic box. */
  double map_y(double); /*function for mapping the periodic box. */

  void usage(void); /*little function to print out the usage of the command
		      line arguments*/
   const double sqrt3 = 1.732050808; /* square root of three, a commonly 
				       occurring constant in an hcp lattice. */
  
  /* lengths of the X and Y sides in angstroms. */
  const double sidelengthX = NND * MAXLATTICE;
  const double sidelengthY = sqrt3 * NND * MAXLATTICE;
  const double epSqr = EPSILON * EPSILON; /*the Square of Epsilon*/
  const double diameter = RADIUS * 2.0; /*the diameter of the Umbrella*/
  const double diameter_sqr = diameter * diameter; /* square of the diameter*/
  const double rsqr = RADIUS * RADIUS; /*the square of the radius */
  const double circleArea = M_PI * rsqr; /*the area of the umbrella*/
  const double squareArea = sidelengthX * sidelengthY; /* the area of the
							  surface */
  const int binsize = sidelengthX / (2.5 * RADIUS); /*the # of bins */
  const int mid_bin = (int)(((double)binsize) / 2.0); /*the middle bin */
  const double dBinX = (sidelengthX + NND) / binsize; /* the x bin width*/
  const double dBinY = (sidelengthY + NND)/ binsize; /* the y bin width*/
  /* note the extra little bit helps catch anyone on the edges. */

  time_t seed; /* seed value for the random number generator. */
  int i, j, k, l; /*my loop counter variables*/
  double mapX = 0.0; /*the x position mapped to the nearest lattice site*/
  double mapY = 0.0; /*the y position mapped to the nearest lattice site*/
  int n = 0; /*the index value for the x coordinate of the lattice site*/
  int m = 0; /*the index value for the y coordinate of the lattice site*/
  double mapDsqr = 0.0; /*the square of the mapped distance to the lattice 
			  site*/
  int stick = 0; /*boolean to let us know if the umbrella sticks to the 
		   surface*/
  double dx = 0.0; /*the difference between x coordinates*/
  double dy = 0.0; /*the difference between y coordinates*/
  double d_sqr; /*the distance between two points squared. */
  struct coords *locations; /*an array of MAXN elements that contain the 
			      coordinates of the attached umbrellas*/
  struct coords temp = {0.0, 0.0, 0.0}; /*temporary storage variable for the 
					  random x, y, and theta drops. */
  int Xbin = 0; /*divides the x locations into binsize bins*/
  int Ybin = 0; /*divides the y locations into binsize bins*/
  int mod_Xbin; /*keeps track of the bin, including periodic boundary
		  conditions */
  int mod_Ybin; /*same as mod_Xbin, only for the y axis. */
  int nLanded[binsize][binsize]; /* keeps track of the success for each bin*/
  long int attempts = 0; /*keeps track of the number of attempts*/
  int success = 0; /*keeps track of the number of successes*/
  int overlap = 0; /*boolean to keep track of any overlaps.*/
  double coverage = 0.0; /*the percent coverage of the surface*/
  double coverageSum = 0.0; /*the sum of the coverages, used to compute the 
			      averages*/
  double avgCoverageOld = 0.0; /*the old average coverage*/
  double avgCoverageNew = -10.0; /*the new average coverage*/
  double threshold = THRESHOLD;
  struct coords spaces[6]; /*these are the coordinates of the bind-able spaces
			     surrounding an hcp lattice site*/

  int visualize = 0; /* boolean to tell if the visualization output is on */
  double nx, ny, nz; /* The normals for the visualization output */
  int visual_count = 0; /*lets the visualization program know how many 
			  particles to expect. */

  int n_legs = 8; /* the number of legs the octopus has. */
  struct coords legs[n_legs]; /* the locations of the legs of the "octopus"*/
  double leg_dx;  /* two variables to map the locations of the legs from the */
  double leg_dy;  /* centers. */

  int surface_scan = 0; /*boolean to tell us to do the surface scan. */
  int covered = 0; /* boolean to see if the scanning point is covered. */
  const double scan_res = 1.0; /*the resolution on a side of the scan point */
  int scan_nx = (int)(sidelengthX / scan_res); /*the number of x scan points*/
  int scan_ny = (int)(sidelengthY / scan_res); /*the number of y scan points*/
  double scan_x;
  double scan_y; /*x and y coordinates of the scanning point. */
  const double scan_area = scan_res * scan_res; /*the area of a scan point */
  double covered_area_tot = 0.0; /*the cumulative covered area. */
  double coverage_true; /*the true percent coverage */

  char *outfile_name = COVERAGE_FILE; /*the name of the output file*/
  char *orient_name = ORIENTATION_FILE; /*the name of the orientation output 
					  file.*/
  char *visual_name = VISUAL_FILE; /*the name of the visualization
                                     output file */

  FILE *visual_file; /*the output file for the visualization locations. */
  FILE *output_file; /* the output file*/
  FILE *orient_file; /* the file of the orientation of the umbrella's*/

  /**************************************************************************
   * The next pointer array is created to keep track of which octopus is 
   * in which bin. The first index is the x coordinate of the bin. The second
   * index is the y coordinate of the bin. The third index is the number of the
   * particle in the bin (first, second, third, etc.)
   ***************************************************************************/
  
  /* a purposeful overestimate of the likely number of octopi per bin. */
  const int maxUPB = 2 * ((int)((dBinX * dBinY) / (circleArea / 8.0)));

  bool worked = true;
  size_t ii, jj;
  struct coords ***inBin = calloc( binsize, sizeof (struct coords) );
  if (!inBin) {
    (void)fprintf(stderr, "Bad calloc\n");
    usage();
  }

  for (ii = 0; worked && ii < binsize; ii++) {
    inBin[ii] = calloc( binsize, sizeof *inBin[ii] );
    worked = (inBin[ii] != NULL);
  }
  // If allocating any inBin[ii] failed, free all inBin[0] through
  // inBin[ii-1], *then* free inBin.  Freeing inBin alone won't free
  // each inBin[ii].
  if ( !worked ) {
    while ( ii-- ) {
      free( inBin[ii] );
    }
    free( inBin );
    (void)fprintf(stderr, "Bad allocation\n");
    usage();
  }
  
  // for each inBin[ii], allocate inBin[ii][jj].  
  for (ii = 0; worked && ii < binsize; ii++ ) {
    for (jj = 0; worked && jj < binsize; jj++ ) {
      inBin[ii][jj] = calloc( maxUPB, sizeof *inBin[ii][jj] );
      worked = (inBin[ii][jj] != NULL );
    }
  }

  // Same deal - if any inBin[ii][jj] allocation failed, free all
  // inBin[ii][0] through inBin[ii][jj-1], *then* free all inBin[0]
  // through inBin[ii], *then* free inBin.
  if ( !worked ) {
    do {
      while ( jj-- )
	free( inBin[ii][jj] );
      free( inBin[ii] );
    } while ( ii-- );
    free( inBin );
    (void)fprintf(stderr, "Bad allocation\n");
    usage();
  }

  /*struct coords inBin[binsize][binsize][maxUPB];*/
  
  /***************************************************************************
   * end declarations
   **************************************************************************/
  
  program_name = argv[0]; /*save the program name in case we need it*/
  
  for( i = 0; i < argc; i++){
    
    if(argv[i][0] =='-'){
      
      /* argv[i][1] is the actual option character */
      
      switch(argv[i][1]){
	
	/* -D => debug option */
	
      case 'D':
	debug = 1;
	break;
	
	/* -L => lattice turned on. */
	
      case 'L':
	lattice = 1;
	break;

	/* -c <name> => the coverage output file name
	 *     [i+1] actually starts the name
	 */
	
      case 'c':
	outfile_name = argv[i+1];
	break;
	
	/* -o <name> => the orientation file name
	 *     [i+1] actually starts the name
	 */
	
      case 'o':
	orient_name = argv[i+1];
	break;
	
	/* -v <name> => the visualization output file name
	 *    [i+1] is the name.
	 */

      case 'v':
	visual_name = argv[i+1];
	break;
	
	/* -V => turns on the visualization output. */

      case 'V':
	visualize = 1;
	break;
	
	/* -S => turns on the surface probe scan. */

      case 'S':
	surface_scan = 1;
	break;

      default:
	(void)fprintf(stderr, "Bad option %s\n", argv[i]);
	usage();
      }
    }
  }
  
  
  /* allocate memory off the heap for locations*/
  
  locations = (struct coords *)calloc(MAXN, sizeof(struct coords));
  
  if(debug){
    (void)fprintf(stderr, "Memory allocation for the large arrays is a go.\n");
  }
  
  
  /*initialize the seed for the random number generator drand48()*/
  
  (void)time(&seed);

  srand48((int)seed);

  /* Initialize the output file*/

  output_file = fopen(outfile_name, "w");
  (void)fprintf(output_file, 
		"#N_Attempts\tPercent_Coverage\tAcceptance_Ratio\n");

  if(debug){
    (void)fprintf(stderr, "wrote the headers to storage file.\n");
  }
  
  /* initialize the arrays*/

  for(i = 0; i < binsize; i++){
    locations[i].x = 0.0;
    locations[i].y = 0.0;
    locations[i].theta = 0.0;
  }

  for(i = 0; i < binsize; i++){
    for(j = 0; j < binsize; j++){
      nLanded[i][j] = 0;
    }
  }

  /* initialize the gap locations */  

  spaces[0].x = 0.0;
  spaces[0].y = NND / sqrt3;
  
  spaces[1].x = NND / 2.0;
  spaces[1].y = NND / (2.0 * sqrt3);
  
  spaces[2].x = NND / 2.0;
  spaces[2].y = -NND / (2.0 * sqrt3);

  spaces[3].x = 0.0;
  spaces[3].y = -NND / sqrt3;
  
  spaces[4].x = -NND / 2.0;
  spaces[4].y = -NND / (2.0 * sqrt3);
  
  spaces[5].x = -NND / 2.0;
  spaces[5].y = NND / (2.0 * sqrt3);
 
  /* initialize the leg locations */

  for( i=0; i < n_legs; i++){
    legs[i].x = RADIUS * cos((double)i * 2.0 * M_PI / (double)n_legs);
    legs[i].y = RADIUS * sin((double)i * 2.0 * M_PI / (double)n_legs);
  }

  if(debug){
    (void)fprintf(stderr, "Initialization of arrays and gap locations, GO.\n");
  }
  
  if(lattice){
    threshold = threshold / 1000.0; /*adjusts the threshold to be harsher 
				      for lattice simulations */
  }


  while(success < MAXN){ /*break loop if we exceed maxn success*/
    
    /*output every 10000 attempts*/
    
    if(attempts % 10000 == 0){
      (void)fprintf(output_file, "%ld \t %lf \t %lf \n", attempts, coverage,
	     ((double)success / (double)attempts));

      if(debug){
	(void)fprintf(stderr, "wrote: %ld\t%lf\t%lf\n", attempts, coverage,
	     ((double)success / (double)attempts));
      }
    }

    /*reset the booleans*/
    
    overlap = 0;
    stick = 0;

    /*pick the first random point and orientation*/
    
    temp.x = sidelengthX * drand48();
    temp.y = sidelengthY * drand48();
    temp.theta = 2.0 * M_PI * drand48();

    

    /*********************************************************************
     * This next section entails the check for the octopus landing
     * onto one of the lattice sites. Since the underlying lattice is
     * hcp, there are two atoms located within each replicated
     * cell. the x coordinates of the first site are as follows:
     * 
     * x1 = n * nnd: where n  is an integer.
     *        
     * y1 = m * sqrt(3) * nnd: where m is an integer.
     *
     * the coordinates of the second site are as follows:
     *
     * x2 = nnd * (2n+1) / 2: again, n is an integer.
     *
     * y2 = nnd * sqrt(3) * (2m+1) / 2: m is an integer.
     *    
     * solving the the equations for n and m respectively, the 
     * following are obtained:
     *       
     * n = x1 / nnd :for site one.;
     *
     * m = y1 / (sqrt(3)*nnd) :for site one.
     *
     * n = (x2 / nnd) - 1/2 :for site two.
     *
     * m = (y2 /(sqrt(3) * nnd)) - 1/2 :for site two.
     * 
     * With these equations, we can see if one of the attachment point
     * of the octopus lies within epsilon of a lattice site. However,
     * in practice, once we calculate n and m with actual data points,
     * they will be double precision numbers. To these decimal
     * numbers, we will add 0.5 before casting them as integers. This
     * will place the nearest lattice point at the origin, and the
     * attachment point at some distance away which will be less than
     * half the distance to the next lattice point.
     *
     **********************************************************************
     *
     * NOTE: in this variation of the program, all of the above
     * applies, however with the exception that now we will test for
     * attachment to a gap. The gap locations are tracked with the
     * spaces array, and we simply place the template array at each
     * lattice site.
     * 
     **********************************************************************/

    if(lattice){
      
      /*check each leg */
      
      for( j = 0; j<n_legs && stick == 0; j++){
	
	/* rotate the leg by the particle's theta. */

	/******************************************************************
	 *
	 * When Cartesian coordinates are used, the 2D rotational
	 * matrix is used to rotate the tracked points before they can
	 * be mapped to the reference frame. The matrix looks as
	 * follows:
	 * 
	 * [ x' ] = [  cos@    sin@ ]  [ x ]
	 * [ y' ]   [ -sin@    cos@ ]  [ y ] 
	 * 
	 * x' and y' are the rotated points, and @ is the angle
	 * through which they will be rotated. The matrix
	 * multiplication can be carried out to give the following:
	 * 
	 * x' = x * cos@ + y * sin@
	 *
	 * y' = y * cos@ - x * sin@
	 *
	 ******************************************************************/
	
	leg_dx = legs[j].x * cos(temp.theta) + legs[j].y * sin(temp.theta);
	leg_dy = legs[j].y * cos(temp.theta) - legs[j].x * sin(temp.theta);
	
	/* Check lattice site 1 */
	
	n = (int)(((temp.x + leg_dx) / NND) + 0.5);
	m = (int)(((temp.y + leg_dy) / (sqrt3 * NND)) + 0.5);
	
	mapX = (temp.x + leg_dx) - n * NND;
	mapY = (temp.y + leg_dy) - m * sqrt3 * NND;
	
	/* cycle over the gaps around the lattice site*/
	
	for(i = 0; i < 6 && stick == 0; i++){
	  dx = mapX - spaces[i].x;
	  dy = mapY - spaces[i].y;
	  mapDsqr = dx * dx + dy * dy;
	  
	  /* check to see if the point is within epsilon of the gap*/
	  
	  if(mapDsqr < epSqr){
	    
	    stick = 1; /* flip the boolean, to be true*/
	    
	    /*set the landing site to be the gap site*/
	    
	    temp.x = (n * NND + spaces[i].x) - leg_dx;
	    temp.y = (m * sqrt3 * NND + spaces[i].y) - leg_dy;
	    
	    /* map the point into the periodic box if neccassary. */
	    
	    if(temp.x > sidelengthX){
	      temp.x = temp.x - sidelengthX;
	    }
	    if(temp.y > sidelengthY){
	      temp.y = temp.y - sidelengthY;
	    }
	    if(temp.x < 0){
	      temp.x = sidelengthX + temp.x;
	    }
	    if(temp.y < 0){
	      temp.y = sidelengthY + temp.y;
	    }
	  }
	}	    
	
	/* check lattice site 2 */
	
	if(stick == 0){
	 
	  n = (int)((temp.x + leg_dx) / NND);
	  m = (int)((temp.y + leg_dy) / (NND * sqrt3));
	  
	  mapX = (temp.x + leg_dx) - (NND * (2 * n + 1)) / 2;
	  mapY = (temp.y + leg_dy) - (NND * sqrt3 * (2 * m + 1)) / 2;
	  
	  /* cycle over the gap points */
	  
	  for(i = 0; i < 6 && stick == 0; i++){
	    dx = mapX - spaces[i].x;
	    dy = mapY - spaces[i].y;
	    mapDsqr = dx * dx + dy * dy;
	    
	    /*see if the point is within epsilon of the gap */
	    
	    if(mapDsqr < epSqr){
	      stick = 1; /*flip the boolean to true*/
	    
	      /* set the landing site to be the gap site */

	      temp.x = (((NND * (2 * n + 1)) / 2) + spaces[i].x) - leg_dx;
	      temp.y = (((NND * sqrt3 * (2 * m + 1)) / 2) + spaces[i].y)
		- leg_dy;
	      
	      /* map the point into the periodic box if neccassary. */
	      
	      if(temp.x > sidelengthX){
		temp.x = temp.x - sidelengthX;
	      }
	      if(temp.y > sidelengthY){
		temp.y = temp.y - sidelengthY;
	      }
	      if(temp.x < 0){
		temp.x = sidelengthX + temp.x;
	      }
	      if(temp.y < 0){
		temp.y = sidelengthY + temp.y;
	      }
	    }
	  }
	}
      }
    }
    
    /*everybody automatically sticks if the lattice is turned off*/
    
    else{
      stick = 1;
    }
    
    /* if stick is true determine Xbin and Ybin for the overlap test */
    
    if(stick){
            
      /* find in which bin our point lies */
      
      Xbin = (int)(temp.x / dBinX);
      Ybin = (int)(temp.y / dBinY);
    }
    
    /* To see if the umbrella overlaps with any others, a 3 x 3 grid of bins,
       centered on the bin in which our point lies, is checked. */
    
    for(i = -1; i < 2 && !overlap && stick; i++){
      for(j = -1; j < 2 && !overlap; j++){
	
	mod_Xbin = (Xbin + i);
	mod_Ybin = (Ybin + j);
	
	/* wrap the bins back into the periodic box. */
	
	if(mod_Xbin < 0){
	  mod_Xbin = binsize + mod_Xbin;
	}
	if(mod_Ybin < 0){
	  mod_Ybin = binsize + mod_Ybin;
	}
	if(mod_Xbin >= binsize){
	  mod_Xbin = mod_Xbin - binsize;
	}	
	if(mod_Ybin >= binsize){
	  mod_Ybin = mod_Ybin - binsize;
	}	

	
	for(k = 0; k < nLanded[mod_Xbin][mod_Ybin] && !overlap; k++){
	  
	  dx = temp.x - inBin[mod_Xbin][mod_Ybin][k].x;
	  dy = temp.y - inBin[mod_Xbin][mod_Ybin][k].y;
	  
	  /* map the positions into the periodic box. */

	  dx = map_x(dx);
	  dy = map_y(dy);
	  d_sqr = dx * dx + dy * dy;
	  
	  /* if they overlap flip the boolean and get on with life. */

	  if(d_sqr < diameter_sqr){
	    overlap = 1;
	  }
	  
	  
	}
      }
    }
    
    
    /* If there was no overlap and the the umbrella stuck, the location
       gets added to the array. */
    
    if(!overlap && stick){
      locations[success].x = temp.x;
      locations[success].y = temp.y;
      locations[success].theta = temp.theta;
      inBin[Xbin][Ybin][nLanded[Xbin][Ybin]] = locations[success];
      nLanded[Xbin][Ybin]++;
      success++;
      
      /* If the number of particles in the bin exceed the max number, then 
	 exit prematurely, and let the user know that maxUPB needs to be 
	 bigger.*/
      
      if(nLanded[Xbin][Ybin] >= maxUPB){
	(void)fprintf(stderr, "maxUPB exceeded in Xbin: %d; Ybin: %d\n", 
		      Xbin, Ybin);
	goto exit;
      }
      
    }
    
    attempts++;
    coverage = (((double)success) * circleArea) / squareArea;
    coverageSum += coverage;
    
    /* check the average coverage every 10,000 steps. */
    
    if(attempts % 10000 == 0){
      avgCoverageOld = avgCoverageNew;
      avgCoverageNew = coverageSum / ((double)attempts);
      
      /*If the average coverage has not changed significantly, then quit. */
      
      if((avgCoverageNew - avgCoverageOld) < threshold){
	goto exit;
      }
    }
  }
    
  
  
  
  
 exit: /* gives us a direct way out of all of the nested loops */
  
  
  if(surface_scan){

    for(l = 0; l < scan_nx; l++){
      
      scan_x = scan_res * (double)l;
      
      for(m = 0; m < scan_ny; m++){
	
	scan_y = scan_res * (double)m;
	
	temp.x = scan_x;
	temp.y = scan_y;
	Xbin = (int)(temp.x / dBinX);
	Ybin = (int)(temp.y / dBinY);
	
	/*reset the boolean */
	
	covered = 0;
	
	for(i = -1; i < 2 && !covered; i++){
	  for(j = -1; j < 2 && !covered; j++){
	    
	    mod_Xbin = (Xbin + i);
	    mod_Ybin = (Ybin + j);
	    
	    /* wrap the bins back into the periodic box. */
	    
	    if(mod_Xbin < 0){
	      mod_Xbin = binsize + mod_Xbin;
	    }
	    if(mod_Ybin < 0){
	      mod_Ybin = binsize + mod_Ybin;
	    }
	    if(mod_Xbin >= binsize){
	      mod_Xbin = mod_Xbin - binsize;
	    }	
	    if(mod_Ybin >= binsize){
	      mod_Ybin = mod_Ybin - binsize;
	    }	
	    
	    for(k = 0; k < nLanded[mod_Xbin][mod_Ybin] && !covered; k++){
	      
	      dx = temp.x - inBin[mod_Xbin][mod_Ybin][k].x;
	      dy = temp.y - inBin[mod_Xbin][mod_Ybin][k].y;
	      
	      /* map the positions into the periodic box. */
	      
	      dx = map_x(dx);
	      dy = map_y(dy);
	      d_sqr = dx * dx + dy * dy;
	      
	      /* if they overlap flip the boolean and get on with life. */
	      
	      if(d_sqr < rsqr){
		covered = 1;
		covered_area_tot += scan_area;
	      }
	    }
	  }
	}
      }
    }
    
    coverage_true = covered_area_tot / squareArea;
  }
  
  /* if we didn't surface scan, then coverage true is just the normal
     coverage*/
    
  else{
    coverage_true = coverage;
  }


  /* write out the final locations */
  
  orient_file = fopen(orient_name, "w");
  (void)fprintf(orient_file,"#successes\tSquare_Area\t total coverage = %lf\n",
		coverage_true);
  (void)fprintf(orient_file,"%d\t%lf\n",success,squareArea);
  (void)fprintf(orient_file, 
		"#x_coord.\t y_Coord.\t Theta\n");
  for(i = 0; i < success; i++){
    (void)fprintf(orient_file,"%lf\t%lf\t%lf\n", locations[i].x, 
		  locations[i].y, locations[i].theta);
  }
  
  /* output the visualization locations if desired. */
  
  if(visualize){
  
    /*the output for visualization */
    
    visual_file = fopen(visual_name, "w");
    
    visual_count = 0;
    
    /*  check out the region around bin: (mid_bin, mid_bin) */
    
    Xbin = mid_bin;
    Ybin = mid_bin;
    
    for(i = -5; i < 6; i++){
      for(j = -5; j < 6; j++){
	
	mod_Xbin = (Xbin + i);
	mod_Ybin = (Ybin + j);
	
	/* wrap the bins back into the periodic box. */
	
	if(mod_Xbin < 0){
	  mod_Xbin = binsize + mod_Xbin;
	}
	if(mod_Ybin < 0){
	  mod_Ybin = binsize + mod_Ybin;
	}
	if(mod_Xbin >= binsize){
	  mod_Xbin = mod_Xbin - binsize;
	}	
	if(mod_Ybin >= binsize){
	  mod_Ybin = mod_Ybin - binsize;
	}	
	
	
	for(k = 0; k < nLanded[mod_Xbin][mod_Ybin]; k++){
	
	  /* the lazy man's way of counting the number of particles that will
	     be written. */
	  
	  visual_count++;
	}
      }
    }
	
    (void)fprintf(visual_file,
		  "#n_particles:\t%d\n",
		  visual_count);

    for(i = -5; i < 6; i++){
      for(j = -5; j < 6; j++){
	
	mod_Xbin = (Xbin + i);
	mod_Ybin = (Ybin + j);
	
	/* wrap the bins back into the periodic box. */
	
	if(mod_Xbin < 0){
	  mod_Xbin = binsize + mod_Xbin;
	}
	if(mod_Ybin < 0){
	  mod_Ybin = binsize + mod_Ybin;
	}
	if(mod_Xbin >= binsize){
	  mod_Xbin = mod_Xbin - binsize;
	}	
	if(mod_Ybin >= binsize){
	  mod_Ybin = mod_Ybin - binsize;
	}	
	
	
	for(k = 0; k < nLanded[mod_Xbin][mod_Ybin]; k++){

	  /* Note: the octopi all face up, so their normals should all
             be the same */
	  
	  nx = 0.0;
	  ny = 0.0;
	  nz = 1.0;
	  
	  /* output the locations */
	  
	  (void)fprintf(visual_file,
			"p:\t%lf\t%lf\t%lf\n"
			"n:\t%lf\t%lf\t%lf\n"
			"r:\t%lf"
			"\n",
			inBin[mod_Xbin][mod_Ybin][k].x,
			inBin[mod_Xbin][mod_Ybin][k].y,
			HANDLE_LENGTH,
			nx, ny, nz,
			RADIUS);
	  
	}
      }
    }

    (void)fclose(visual_file);
  
  }
  
  
  /* close the output files. */

  (void)fclose(orient_file);
  (void)fclose(output_file);
  
  /* free up the allocated memory */
  
  free(locations);
  // free in reverse order of allocation
  for ( ii = 0; i < binsize; ii++ ) {
    for ( jj = 0; j < binsize; jj++ )
      free( inBin[ii][jj] );
    free( inBin[ii] );
  }
  free( inBin );
   
  return (0);
}


/***************************************************************************
 * prints out the usage for the command line arguments, then exits.
 ***************************************************************************/

void usage(){
  
  (void)fprintf(stderr, 
		"The proper usage is: %s [options]\n\n"
		"Options:\n"
		"   -D         turn on the debugging output\n"
		"   -o <name>  the name of the orientation output file\n"
		"   -c <name>  the name of the coverage output file\n"
		"   -L         turns on the lattice constraint\n"
		"   -V         generate output locations and orientations\n"
		"              for visualization\n"
		"   -v <name>  the name for the visualization file\n"
		"   -S         turn on the surface probe scan for coverage\n",
		program_name);
  exit(8);
}



/***********************************************************************
 * maps the x direction into the periodic box.
 **********************************************************************/

double map_x(unmapped)
     double unmapped;
{ 
  const double length_X = NND * MAXLATTICE;/* length of the x side of the 
					      box.*/
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


/*************************************************************************
 *  responsible for mapping the y direction in the periodic box. 
 *************************************************************************/

double map_y(unmapped)
     double unmapped;
{
  const double sqrt3 = 1.732050808; /* square root of three, a commonly 
				       occurring constant in an hcp 
				       lattice. 
				       (the underlying lattice in the
				       simulation) */
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
