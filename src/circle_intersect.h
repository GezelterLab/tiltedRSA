/*
 * circle_intersect.h    1.0  2001/05/16
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



/*structure to keep coordinates together*/

struct coords {
  double x; /*the x coordinate */
  double y; /*the y coordinate */
  double theta; /*orientation of the umbrella to the plane 
		  NOTE: the angle is measured relative to the y axis.
		  ie. 0 degrees has the face point straight up. (on the xy 
		  plane) */
};

/**************************************************************************
 * function to tell if the handle of the umbrella is overlapping another's 
 * top.
 *  
 * usage: handle_compare(&set_location, &test_location)
 *       
 *        struct coords set_location => the pre-existing location to which to 
 *                                      check against.
 *        struct coords test_location => The location to check.
 *
 *
 * returns: 1 => if overlap is true.
 *          0 => if overlap is false.
 *************************************************************************/

extern int handle_compare(struct coords, struct coords);

/**************************************************************************
 * Function to check whether the top of one umbrella overlaps with that 
 * of another umbrella.
 * 
 * usage: circle_compare(&location1, &location2)
 *  
 *        struct coords location1 => the first umbrella.
 *        struct coords location2 => the second umbrella.
 *
 * 
 * returns: 1 => if overlap is true
 *          0 => if overlap is false
 ***************************************************************************/

extern int circle_compare(struct coords, struct coords);


