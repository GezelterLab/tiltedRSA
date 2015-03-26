/*
 * visualize.c    1.0  2001/05/16
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



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parameters.h"

struct particle {
  
  /* the cartesian coords */

  double x;
  double y;
  double z;

  /* the surface normal */

  double nx;
  double ny;
  double nz;

  /* the radius */
  
  double radius;
};

char *program_name; /*gives us the name that the program was run as*/

main(argc, argv)
     int argc;
     char *argv[];
{

  void usage(void); /*little function to print out the usage of the command
		      line arguments*/

  int i; /*loop counters */
  struct particle *umbrellas; /*the array of umbrellas. */
  FILE *input_file; /*the input file */
  FILE *povray_script; /*the povray output script */
  FILE *pov_header; /*the pov_header file. */
  double dx, dy, d2; /*variables used in the distance function */
  int n_particles; /* the number of particles */
  const double max_r = 500.0; /* the maximum radius to observe */
  const double rsqr = (max_r * max_r); /*r squared. */
  const double sqrt3 = 1.732050808; /* square root of three, a commonly 
				       occuring constant in an hcc lattice. */
  /* the next two are the center coordinates of the particle cluster */
  const double center_x = NND * MAXLATTICE / 2.0;
  const double center_y = NND * sqrt3 * MAXLATTICE / 2.0;
  double theta; /*the angle of the normal on thge xy plane */
  double gradient = 0.0; /* will be used to determine the color gradient*/
  double r, g, b; /* the red, green, blue components of color. */
  double rn, gn, bn;/* adjusted colors. */
  int neg = 0; /*boolean to know if the angle is negative */
  const double pi_3 = M_PI / 3.0; /* pi / 3 */
  const double two_pi_3 = 2.0 * pi_3; /* 2*pi / 3 */
  double nx, ny, n2; /* modified normals. */

  int generate_header = 0; /* boolean tells if the header file should 
			  be generated*/
  double height; /* used in the pov_header to place lights and the camera. */
  int n_lights = 3; /*specifies the number of lights in the header. */
  double light_theta; /*the angle by which the lights are dispersed */

  int color_fun = 0; /* boolean to turn on mr. happy gradient. */
  double brightness = 0.750; /* 0 -> 1 tones down mr. happy gradient. */
  char *visual_name = VISUAL_FILE; /* the name of the file containing
				      the visualization locations*/
  char *pov_name = POV_NAME; /* The name of the povray output script. */

  
  umbrellas = (struct particle *)calloc(MAX_VISUALIZE, 
					sizeof(struct particle));

  /**************************************************************
   *
   * allows the user to specify input and output files. 
   *
   * Usage:
   *        visualize -I <input file> -O <output script>
   *
   ***************************************************************/
  
  program_name = argv[0]; /*save the program name for the usage blurb. */

  for(i = 0; i < argc; i++){
    
    /* see if we have any option charecters. */

    if(argv[i][0] == '-'){
      
      switch(argv[i][1]){
	
	/* -I <input file> => the input file */
	
      case 'I':
	visual_name = argv[i+1];
	break;

	/* -O <output file> => the output script. */
	
      case 'O':
	pov_name = argv[i+1];
	break;

	/* turn on the pov_header generation. */
	
      case 'g':
	generate_header = 1;
	break;

	/*turn on mr. happy color gradient. */

      case 'c':
	color_fun = 1;
	break;

      default:
	(void)fprintf(stderr, "bad option %s\n", argv[i]);
	usage();
      }
    }
  }
    
  input_file = fopen(visual_name, "r");
  
  (void)fscanf(input_file, "#n_particles:\t%d\n", &n_particles);
  
  for(i = 0; i < n_particles; i++){
    
    (void)fscanf(input_file,
		 "p:\t%lf\t%lf\t%lf\n"
		 "n:\t%lf\t%lf\t%lf\n"
		 "r:\t%lf"
		 "\n",
		 &umbrellas[i].x,
		 &umbrellas[i].y,
		 &umbrellas[i].z,
		 &umbrellas[i].nx,
		 &umbrellas[i].ny,
		 &umbrellas[i].nz,
		 &umbrellas[i].radius);
  }

  (void)fclose(input_file);

  povray_script = fopen(pov_name, "w");

  (void)fprintf(povray_script, "#include \"pov_header.pov\"\n\n");

  for(i=0; i< n_particles; i++){

    dx = umbrellas[i].x - center_x;
    dy = umbrellas[i].y - center_y;
    
    d2 = dx * dx + dy * dy;
    
    if(d2 < rsqr){
      
      /* mr. happy color fun. */
      
      if(color_fun){
	
	/* calculate the new normals */
	
	nx = umbrellas[i].nx;
	ny = umbrellas[i].ny;
	
	n2 = nx * nx + ny * ny;
	
	nx = nx / sqrt(n2);
	ny = ny / sqrt(n2); 
	
	/*initialize the components and the boolean */

	r = 0.0;
	g = 0.0;
	b = 0.0;
	neg = 0;
	
	/* the good stuff */

	theta = acos(nx);
	if(ny < 0){
	  neg = 1;
	}
	 
	/*if the angle is on the bottom of the circle */
	
	if(neg){
	  
	  theta = two_pi_3 - theta;
	  
	  if(theta < 0){
	    gradient = -(theta / pi_3);
	    r = 1.0 - gradient;
	    g = gradient;
	  }

	  else{
	    gradient = (theta / two_pi_3);
	    r = (1.0 - gradient);
	    b = gradient;
	  }
	}

	/* here the top of the circle is considered. */
	
	else{
	  theta = two_pi_3 - theta;
	  
	  if(theta < 0){
	    gradient = -(theta / pi_3);
	    g = (1.0 - gradient);
	    r = gradient;
	  }
	  
	  else{
	    gradient = theta / two_pi_3;
	    g = 1.0 - gradient;
	    b = gradient;
	  }
	}
      }
      
      /* no mr. happy color fun. :(  */
      
      else{
	r = 0.0;
	g = 0.5;
	b = 1.0;
      }

      /* adjust the brightness */
      
      r = r * brightness;
      g = g * brightness;
      b = b * brightness;

      /* adjust the diferences in the brightness of color components. */

      rn = 3.240479 * r - 1.53715 * g - 0.498535 * b;
      gn = -0.969256 * r + 1.875992 * g + 0.041556 * b;
      bn = 0.055648 * r - 0.204043 * g + 1.057311 * b;
      
      (void)fprintf(povray_script,
		    "disc {\n"
		    "   < %lf, %lf, %lf>,\n"
		    "   < %lf, %lf, %lf>,\n"
		    "   %lf\n"
		    "   texture{\n"
		    "     pigment { color rgb < %lf, %lf, %lf> }\n"
		    "   }\n"
		    "}\n"
		    "\n",
		    dx, dy, umbrellas[i].z,
		    umbrellas[i].nx, umbrellas[i].ny, umbrellas[i].nz,
		    umbrellas[i].radius,
		    rn, gn, bn);
    }
  }

  (void)fclose(povray_script);

  /* Generate the pov_header.pov file */

  if(generate_header){

    pov_header = fopen("pov_header.pov", "w");
    
    height = 15.0 * RADIUS;
    light_theta = M_PI * 2.0 / (double)n_lights;
        
    (void)fprintf(pov_header,
		  "#include \"shapes.inc\"\n"
		  "#include \"colors.inc\"\n"
		  "\n"
		  "camera {\n"
		  "  location <0, 0, %lf>\n"
		  "  look_at <0, 0, 0>\n"
		  "}\n"
		  "\n",
		  height);
    
    for(i = 0; i < n_lights; i++){

     (void)fprintf(pov_header,
		   "light_source {\n"
		   "  <%lf, %lf, %lf>\n"
		   "  color rgb <1.0, 1.0, 1.0>\n"
		   "}\n",
		   2.0 * height * cos((double)i * light_theta),
		   2.0 * height * sin((double)i * light_theta),
		   height);
    }
    
    (void)fprintf(pov_header,
		  "\n"
		  "plane {\n"
		  "  z, 0.0\n"
		  "  texture {\n"
		  "    pigment { checker color Gray color Black scale 5}\n"
		  "  }\n"
		  "}\n");
  
    (void)fclose(pov_header);
  }

}

/***************************************************************************
 * prints out the usage for the command line arguments, then exits.
 ***************************************************************************/

void usage(){
  
  (void)fprintf(stderr, 
		"The proper usage is: %s [options]\n\n"
		"Options:\n"
		"   -I <input name>   the name of input locations file\n"
		"   -O <output name>  the name of the output script\n"
		"   -g                generate the pov_header.pov file\n"
		"                       -contains info on lighting and\n"
		"                        camera locations\n"
		"   -c                Turns on the color gradient.\n"
		"                       -note: will work only for the\n"
		"                              tilted umbrellas.\n",
		program_name);
  exit(8);
}
