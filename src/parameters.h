/*
 * paramters.h    1.0  2001/05/16
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
 * The following is a list of parameters for the RSA simulations.
 ****************************************************************************/


#define MAXN 5000000 /* The maximum number of successful landings */
#define MAXLATTICE 4000 /* The maximum number of lattice points in the
			    x or y direction */
#define NND 2.3 /* The nearest neighbor distance of the lattice */
#define EPSILON 0.1 /* the tolerance of accepting a bond to a lattice
                       space */
#define RADIUS 14.0 /* The radius of the dropped objects */
#define THRESHOLD 0.000004 /*the threshold for stopping the
                              simulation */
#define POV_NAME "visual.pov" /* The visualization pov-ray script*/
#define VISUAL_FILE "locations.dat" /*the output file for the
                                      visualization of the particles
                                      on the surface. */
#define COVERAGE_FILE "coverage.dat" /*the output file containing the
				      coverage data*/
#define ORIENTATION_FILE "orientation.dat" /*the output file containing
					      the orientation of the
					      particles*/
#define PHI 1.911135531 /*109.5 degrees in radians, the angle the
			  umbrella handle is bent. note: keep between
			  pi and pi/2 */
#define HANDLE_LENGTH 5.0 /*the length of the handle*/
#define max_bin 10000 /*the max number of bins in the histogram. */
#define MAX_VISUALIZE 10000 /* the max number of particles to
                               visualize*/
