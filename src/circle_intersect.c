/*
 * circle_intersect.c    1.0  2001/05/16
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


/**************************************************************************
 * Collected here are the functions used to perform the more complex
 * tests for overlap in the umbrella simulations.
 **************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parameters.h"
#include "circle_intersect.h"


struct line_element { /* used to store the x and y elements of the parametric
			 lines */
  double cnst; /* the constant part */
  double coeff; /*the coefficient of the parameter. */
};

struct roots { /*used to pass the roots of the parametric equation 
		 back and forth between functions. */
  double root1;
  double root2;
};

struct point { /*structure to hold the coordinates of a point in 3 space. */
  double x;
  double y;
  double z;
};


/**************************************************************************
 * This function is used to determine whether the handle of a test
 * umbrella lies within the elliptical projection of a previously
 * attached umbrella.
 **************************************************************************/

int handle_compare(set_ptr, test_ptr)
     struct coords *set_ptr;
     struct coords *test_ptr;
{
  double map_x(double); /*function for mapping the periodic box. */
  double map_y(double); /*function for mapping the periodic box. */

  double inner_angle = M_PI - PHI; /*the angle between the umbrella
				     top and the handle */
  double minor; /* the radius of the minor axis. (the radius of the major axis
		   is simply the radius of the umbrella.) */
  double focus_pos; /*the positive focus of the elliptical projection */
  double focus_neg; /*the negative focus of the elliptical projection.*/
  double x_rot; /*the test umbrella's x position rotated into the projection's 
		  reference frame.*/
  double y_rot; /*the test umbrella's y position rotated into the projection's 
		  reference frame.*/
  double dx; /*temp variable for storing change in x*/
  double dy; /*temp variable for storing the change in y*/
  double distance_sum = 0.0; /*the sum of the distances from the point to the 
			       two foci. */

  

  /* first calculate the minor axis, and in turn calculate the foci */

  minor = RADIUS * sin(inner_angle);

  focus_pos = sqrt(RADIUS * RADIUS - minor * minor);
  focus_neg = -sqrt(RADIUS * RADIUS - minor * minor);

  /* scale the test umbrella's coordinates to share the same origin as the 
     set umbrella*/

  dx = test_ptr->x - set_ptr->x;
  dy = test_ptr->y - set_ptr->y;

  /* evoke the periodic box. */
  
  dx = map_x(dx);
  dy = map_y(dy);

  /* rotate the new x and y into the projection reference frame */
  
  x_rot = dx * cos(set_ptr->theta) + dy * sin(set_ptr->theta);
  y_rot = dy * cos(set_ptr->theta) - dx * sin(set_ptr->theta);

  /* sum the distances from the test point to the two foci*/

  dx = x_rot - focus_pos;
  dy = y_rot;
  
  distance_sum += sqrt(dx*dx + dy*dy);
  
  dx = x_rot - focus_neg;
  dy = y_rot;
  
  distance_sum += sqrt(dx*dx + dy*dy);

  /*if the sum of the distance is less than 2*Radius then the point lies within
    the ellipse. */

  if(distance_sum < (2 * RADIUS)){
    return(1);
  }
  else{
    return(0);
  }
}

/****************************************************************************
 * The following function will test the tops of two given umbrellas
 * for intersection in three space. The algorithm firsts computes the
 * parametric line that arises from the intersection of two planes
 * (the tops of the umbrellas). This line is then tested for
 * intersection with the two umbrella tops. If roots do exist, then
 * the final test is whether the line sequentially enters and leaves
 * one umbrella top before it intersects the second.
 ****************************************************************************/

int circle_compare(site1_ptr, site2_ptr)
     struct coords *site1_ptr;
     struct coords *site2_ptr;
{
  double map_x(double); /*function for mapping the periodic box. */
  double map_y(double); /*function for mapping the periodic box. */

  /* function to find the roots of intersection between a parametric line,
     and a sphere. */
  int find_roots(struct line_element*, struct line_element*, 
		 struct line_element*, struct point*, struct roots*); 

  double n1x,n1y,n1z; /*temporary variables to hold the x, y, and z 
			coefficients of the vector normal to the plane of 
			the first umbrella's top.*/
  double n2x,n2y,n2z; /*temporary variables to hold the x, y, and z 
			coefficients of the vector normal to the plane of 
			the second umbrella's top.*/
  double p1x,p1y,p1z; /*temporary variables to hold the positions of the 
			center of the first umbrella */
  double p2x,p2y,p2z; /*temporary variables to hold the positions of the 
			center of the first umbrella */
  double denominator; /*temporary value used to calculate the the parametric
			line. */
  double dx; /*temp variable for storing change in x*/
  double dy; /*temp variable for storing the change in y*/
  double d1; /* the offset term of the plane of the 1st umbrella top. */
  double d2; /* the offset term of the plane of the 2nd umbrella top. */
  struct line_element x_element; /* the x element. */
  struct line_element y_element; /* the y element. */
  struct line_element z_element; /* the z element. */
  struct point sphere_center; /*the center of the sphere to be tested*/
  struct roots intersect_1; /*the points of intersection on the 1st sphere. */
  struct roots intersect_2; /*the points of intersection on the 2nd sphere. */
  int roots_exist = 0; /*tells us if roots exist */
  double a,b,c,d; /*temporary variables for testing overlap of the circles */

  /* this is the square of the diameter, also the closest approach for two flat
     circles (phi = pi/2 ). */

 
  /* set the points */
  
  p1x = site1_ptr->x;
  p1y = site1_ptr->y;
  p1z = HANDLE_LENGTH;
  
  p2x = site2_ptr->x;
  p2y = site2_ptr->y;
  p2z = HANDLE_LENGTH;

  /*evoke periodic boundary conditions. */

  dx = p2x - p1x;
  dy = p2y - p1y;

  p1x = 0.0;
  p2x = map_x(dx);

  p1y = 0.0;
  p2y = map_y(dy);

  
  /*set the normals*/
  
  n1x = cos(M_PI - PHI) * cos(M_PI_2 + site1_ptr->theta);
  n1y = cos(M_PI - PHI) * sin(M_PI_2 + site1_ptr->theta);
  n1z = sin(M_PI - PHI);
  
  n2x = cos(M_PI - PHI) * cos(M_PI_2 + site2_ptr->theta);
  n2y = cos(M_PI - PHI) * sin(M_PI_2 + site2_ptr->theta);
  n2z = sin(M_PI - PHI);
  
  /* set the offset terms and the denominator of the parametric line */
  
  d1 = n1x * p1x + n1y * p1y + n1z * p1z;
  d2 = n2x * p2x + n2y * p2y + n2z * p2z;
  
  denominator = n1x*n2y - n1y*n2x + n1z*n2x - n1x*n2z + n1y*n2z - n1z*n2y;
  
  /* set the elements of the parametric line */
  
  x_element.cnst = ((n2y - n2z) * d1 / denominator) + 
    ((n1z - n1y) * d2 / denominator);
  x_element.coeff = n1y*n2z - n1z*n2y;
  
  y_element.cnst = ((n2z - n2x) * d1 / denominator) + 
    ((n1x - n1z) * d2 / denominator);
  y_element.coeff = n1z*n2x - n1x*n2z;
  
  z_element.cnst = ((n2x - n2y) * d1 / denominator) +
    ((n1y - n1x) * d2 / denominator);
  z_element.coeff = n1x*n2y - n1y*n2x;
  
  /*set the 1st sphere's center */
  
  sphere_center.x = p1x;
  sphere_center.y = p1y;
  sphere_center.z = p1z;
  
  /* find the points of intersection with the first circle and the line*/
  
  roots_exist = find_roots(&x_element, &y_element, &z_element, 
			   &sphere_center, &intersect_1);
  
  /* if the circle was indeed intersected, then check the intersection of 
     the second circle */
  
  if(roots_exist){
    
    /*set the 2nd sphere's center */
    
    sphere_center.x = p2x;
    sphere_center.y = p2y;
    sphere_center.z = p2z;
    
    /* find the points of intersection of the second circle with the line */
    
    roots_exist = find_roots(&x_element, &y_element, &z_element, 
			     &sphere_center, &intersect_2);
    
    /* if the second circle was intersected, then test for overlap of the 
       two circles */
    
    if(roots_exist){
      
      a = intersect_1.root1;
      b = intersect_1.root2;
      c = intersect_2.root1;
      d = intersect_2.root2;
      
      /*test to make sure neither set of points fall within the other
	set note: the ! sends only the cases that *DO* overlap into
	the if loop.*/
      
      if(!((a<c && b<c && a<d && b<d) || (a>c && b>c && a>d && b>d))){
	
	return(1); /*overlap true returned. */
      }
    }
  }
  
  return(0); /* overlap false returned */
}

/***************************************************************************
 * This function calculates the intersection points of the line with
 * the two umbrella tops. Note that the tops are thought of as
 * sphere's. This is because the line is already constrained to lie in
 * the plane of the umbrella's top. This then allows for the
 * relatively simple test of a line in three space intersecting with a
 * sphere. The roots that are calculated in this function are the
 * actual parameter values for the line at intersection. This makes
 * the sequential test easier in the circle_compare algorithm.
 ****************************************************************************/

static int find_roots(x_element_ptr, y_element_ptr, z_element_ptr,
		      center_ptr, intersect_ptr)
     struct line_element *x_element_ptr; /* the pointer to the x element
					    of the parametric line. */
     struct line_element *y_element_ptr;  /* the pointer to the y element
					    of the parametric line. */
     struct line_element *z_element_ptr;  /* the pointer to the x element
					    of the parametric line. */
     struct point *center_ptr; /*the pointer to the center of the sphere*/
     struct roots *intersect_ptr; /*the pointer to the roots of intersection*/
{
  
  double mx, my, mz; /*the coeff.'s for each line_element */
  double bx, by, bz; /*the cnst's for each line_element */
  double x0, y0, z0; /*the center of the sphere */
  double a, b, c; /*the coefficients of quadratic equation 
		    (a*x^2 + b*x + c = 0) */

  /* set the coeff.'s of the parametric line */
  
  mx = x_element_ptr->coeff;
  my = y_element_ptr->coeff;
  mz = z_element_ptr->coeff;
  
  /*set the cnst's of the parametric line*/
  
  bx = x_element_ptr->cnst;
  by = y_element_ptr->cnst;
  bz = z_element_ptr->cnst;

  /* set the center of the sphere*/
  
  x0 = center_ptr->x;
  y0 = center_ptr->y;
  z0 = center_ptr->z;
  
  /* set the coefficients of the quadratic equation. */

  a = mx*mx + my*my + mz*mz;
  
  b = 2.0*mx*(bx - x0) + 2.0*my*(by - y0) + 2.0*mz*(bz - z0);

  c = bx*bx + by*by + bz*bz + x0*x0 + y0*y0 + z0*z0 - 2.0*x0*bx - 2.0*y0*by -
    2.0*z0*bz - RADIUS*RADIUS;

  /* test to make sure that there are roots */
  
  if((b*b - 4.0*a*c) > 0){
    
    intersect_ptr->root1 = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);
    intersect_ptr->root2 = (-b - sqrt(b*b - 4.0*a*c)) / (2.0*a);

    return(1); /*return that roots were found */
  }
  
  return(0); /*return that roots did not exist. */
}


static double map_x(unmapped)
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


/* responsible for mapping the y direction in the periodic box. */

static double map_y(unmapped)
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
