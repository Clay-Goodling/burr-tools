/* Burr Solver
 * Copyright (C) 2003-2007  Andreas R�ver
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
#ifndef __SEPARATOR_H__
#define __SEPARATOR_H__

#include "Layouter.h"

#include "FL/Fl_Group.H"

// a widget to separate 2 groups
class LSeparator_c : public Fl_Group,  public layoutable_c {

public:

  LSeparator_c(int x, int y, int w, int h, const char * label, bool button);

  virtual void getMinSize(int *width, int *height) const {
    *width = 10;
    *height = 10;
  }
};

#endif