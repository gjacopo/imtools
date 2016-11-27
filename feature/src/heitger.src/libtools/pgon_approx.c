/* 
 * $Log: pgon_approx.c,v $
 * Revision 1.2  1995/10/02  12:25:57  ohenri
 * *** empty log message ***
 *
 * Revision 1.1  1995/06/21  12:58:19  ohenri
 * Initial revision
 *
 * Revision 1.1  1995/06/21  12:58:19  ohenri
 * Initial revision
 * 
 *
 */

/* pgon_approx.c  -  approximate a curve given as a point list by a polygon.

    Copyright (C) 1992  Christian Brechbuehler

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 1, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    Christian Brechbuehler			brech@vision.ethz.ch
    IKT / Bildwissenschaften
    ETH Zentrum, Gloriastrasse 35
    CH-8092 Zurich


    11 Oct 88	creation				brech@vision.ethz.ch
*/

#include <stdio.h>
#include <stdlib.h>
#include "tools.h"

/* Approximiert die als Punktliste der Laenge npt gegebene Kurve curve durch
 * einen Polygonzug  mit maximaler normaler Abweichung toleranz.
 * Gibt die Anzahl der Linien  zurueck.
 */

int polygon_approx(Punkt  *curve,int npt,float toleranz,float exponent,
		   Punkt *polygon)
/* An interface that makes sure the end of the last span, which is not the
 * start of any other span, will be written to _polygon_. It is the caller's
 * responsability to provide sufficient storage space in _polygon_.
 */
{
  int	n_spans= polygon_approx_i(curve, npt, toleranz, exponent, polygon);
  polygon[n_spans]= curve[npt-1];
  return n_spans;
}

float Tolerance(float d,float e)
{
  if (d==0.0) fprintf(stderr,"kuken ...\n");
  return flog(10.0*fpow(d,e));
}

int polygon_approx_i(Punkt *curve,int npt,float toleranz,float exponent, 
		     Punkt *polygon)
/* Recursive procedure for subdividing the curve into sections that are well
 * enough (i.e., inside _toleranz_) represented by a straight line (a span of
 * a polygon. The starting points of all spans are written to _polygon_, and
 * the ending point is implicitly given as the starting point of the next
 * span. If some application (e.g., with closed curves) does not need the end
 * point of the last span, which is always the last point of _curve_, it may
 * call this function. Normally, "polygon_approx" should be used instead.
 */
{
  Punkt *ende, *p_max, *p;
  int nx, ny, dq_max, d, dq, diff, teil1, index;
  float tol;

  ende= &curve[npt-1]; 
  nx= ende ->y - curve->y;
  ny= curve->x - ende ->x;

  if (toleranz == 0.0) {
    if (nx == 0 && ny == 0) tol=Tolerance((float) npt,exponent);
    else tol=Tolerance(fsqrt(nx*nx+ny*ny),exponent);
  }
  else tol=toleranz;

  if (nx || ny){				  /* open curve */
    dq_max= (nx*nx + ny*ny) * tol*tol;
    diff= ende->x*nx+ ende->y*ny;
    p_max= 0;
    for (p= curve; p != ende; p++){
      d= p->x*nx + p->y*ny - diff;
      if ( (dq= d*d) > dq_max) { dq_max= dq; p_max= p;}
    }
  }
  else {                                          /* closed curve or 1 point */
    int dx, dy;
    dq_max= tol * tol;
    p_max= 0;
    for (p= curve; p != ende; p++) {
      dx= p->x - curve->x; dy= p->y - curve->y;
      dq= dx*dx + dy*dy;
      if (dq > dq_max) { dq_max= dq; p_max= p;}
    }
  }
  
  if (p_max){                                     /* above tol -> split */
    index= p_max - curve;
    teil1= polygon_approx_i(curve, index+1, toleranz, exponent, polygon);
    return teil1 + polygon_approx_i(p_max, npt-index, toleranz, exponent, polygon+teil1);
  }
  /* else */                                      /* toleranz was observed */
  *polygon= *curve;
  return 1;
} /* END polygon_approx */
