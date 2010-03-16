/* Burr Solver
 * Copyright (C) 2003-2009  Andreas Röver
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include "stl.h"

#include <math.h>
#include <string.h>

/** \page STL Surface Tessellation Language
 *
 * STL is a file format allowing to define 3D-Objects using just triangles. The
 * format is used for 3D rapid prototyping machines.
 *
 * See http://en.wikipedia.org/wiki/STL_(file_format) for details in the file format.
 *
 * BurrTools supports text mode STL as well as binary mode STL file export. The working
 * is again similar to all. You have the base class and a derived class for concrete
 * grid types (e.g. one for bricks, one for spheres...)
 *
 * The concrete classes do the grid dependend stuff and add lots of triangles to the file
 */

#if defined(WIN32) || defined(__APPLE__)
const char * basename(const char * name) {
  const char * res1 = strchr(name, '/');
  const char * res2 = strchr(name, '\\');

  const char * res = res1>res2 ? res1 : res2;

  if (res == 0)
    res = name;
  else
    res++;

  return res;
}
#endif

void stlExporter_c::open(const char * name) {

  char fname[1000];
  snprintf(fname, 1000, "%s.stl", name);

  const char * title = basename(name);

  if (binaryMode) {

    f = fopen(fname,"wb");

    if (!f) throw new stlException_c("Could not open file");

    int pos = 0;

    for (int i = 0; i < 84; i++) {
      if (fwrite(title+pos, 1, 1, f) != 1)
        throw stlException_c("Could not write to file");
      if (title[pos]) pos++;
    }

    triangleCount = 0;

  } else {

    f = fopen(fname,"w");

    if (!f) throw new stlException_c("Could not open file");

    fprintf(f, "solid %s\n", title);

  }
}

void stlExporter_c::close(void) {

  if (binaryMode) {

    // write out the triangle count into the header
    fseek(f, 80, SEEK_SET);
    if (fwrite(&triangleCount, 1, 4, f) != 1)
      throw stlException_c("Could not write to file");

  } else {

    fprintf(f, "endsolid\n");

  }

  fclose(f);
}

#define Epsilon 1.0e-5

void stlExporter_c::outTriangle(
        /* the 3 vertexes of the triangle */
        double x1, double y1, double z1,
        double x2, double y2, double z2,
        double x3, double y3, double z3,
        /* a point that is on the _INSIDE_ of the plane of thr triangle
         * this is used to calculate the proper normal and reorientate the
         * points anti clockwise
         */
        double xp, double yp, double zp) {

  double sx = x1;
  double sy = y1;
  double sz = z1;
  double v1x = x2-x1;
  double v1y = y2-y1;
  double v1z = z2-z1;
  double v2x = x3-x1;
  double v2y = y3-y1;
  double v2z = z3-z1;

  float nx = v1y*v2z - v1z*v2y;
  float ny = v1z*v2x - v1x*v2z;
  float nz = v1x*v2y - v1y*v2x;

  float l = sqrt(nx*nx+ny*ny+nz*nz);

  if (l < Epsilon)
    return;

  nx /= l;
  ny /= l;
  nz /= l;

  if (binaryMode) {

    float d;

    bool wrOk = true;

    if (nx*(xp-sx)+ny*(yp-sy)+nz*(zp-sz) > 0) {

      d = -nx; wrOk &= fwrite(&d, 1, 4, f) == 1; d = -ny; wrOk &= fwrite(&d, 1, 4, f) ==1; d = -nz; wrOk &= fwrite(&d, 1, 4, f) == 1;
      d = x1; wrOk &= fwrite(&d, 1, 4, f); d = y1; wrOk &= fwrite(&d, 1, 4, f); d = z1; wrOk &= fwrite(&d, 1, 4, f);
      d = x3; wrOk &= fwrite(&d, 1, 4, f); d = y3; wrOk &= fwrite(&d, 1, 4, f); d = z3; wrOk &= fwrite(&d, 1, 4, f);
      d = x2; wrOk &= fwrite(&d, 1, 4, f); d = y2; wrOk &= fwrite(&d, 1, 4, f); d = z2; wrOk &= fwrite(&d, 1, 4, f);

    } else {

      d = nx; wrOk &= fwrite(&d, 1, 4, f); d = ny; wrOk &= fwrite(&d, 1, 4, f); d = nz; wrOk &= fwrite(&d, 1, 4, f);
      d = x1; wrOk &= fwrite(&d, 1, 4, f); d = y1; wrOk &= fwrite(&d, 1, 4, f); d = z1; wrOk &= fwrite(&d, 1, 4, f);
      d = x2; wrOk &= fwrite(&d, 1, 4, f); d = y2; wrOk &= fwrite(&d, 1, 4, f); d = z2; wrOk &= fwrite(&d, 1, 4, f);
      d = x3; wrOk &= fwrite(&d, 1, 4, f); d = y3; wrOk &= fwrite(&d, 1, 4, f); d = z3; wrOk &= fwrite(&d, 1, 4, f);

    }

    if (!wrOk)
      throw stlException_c("Could not write to file");

    // attribute
    int i = 0;
    if (fwrite(&i, 1, 2, f) != 1)
      throw stlException_c("Could not write to file");

    triangleCount++;

  } else {

    if (nx*(xp-sx)+ny*(yp-sy)+nz*(zp-sz) > 0) {

      fprintf(f,"  facet normal %9.4e %9.4e %9.4e\n", -nx, -ny, -nz);
      fprintf(f,"    outer loop\n");
      fprintf(f,"      vertex %9.4e %9.4e %9.4e\n", x1, y1, z1);
      fprintf(f,"      vertex %9.4e %9.4e %9.4e\n", x3, y3, z3);
      fprintf(f,"      vertex %9.4e %9.4e %9.4e\n", x2, y2, z2);
      fprintf(f,"    endloop\n");
      fprintf(f,"  endfacet\n");

    } else {

      fprintf(f,"  facet normal %9.4e %9.4e %9.4e\n", nx, ny, nz);
      fprintf(f,"    outer loop\n");
      fprintf(f,"      vertex %9.4e %9.4e %9.4e\n", x1, y1, z1);
      fprintf(f,"      vertex %9.4e %9.4e %9.4e\n", x2, y2, z2);
      fprintf(f,"      vertex %9.4e %9.4e %9.4e\n", x3, y3, z3);
      fprintf(f,"    endloop\n");
      fprintf(f,"  endfacet\n");
    }
  }
}

