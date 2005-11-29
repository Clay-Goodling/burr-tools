/* Burr Solver
 * Copyright (C) 2003-2005  Andreas Röver
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
#ifndef __LAYOUTER_H__
#define __LAYOUTER_H__

#include <FL/Fl_Widget.h>
#include <FL/Fl_Group.h>

#include <FL/fl_draw.H>

#include <FL/Fl_Button.h>
#include <FL/Fl_Slider.h>
#include <FL/Fl_Box.h>
#include <FL/Fl_Check_Button.h>
#include <FL/Fl_Round_Button.h>
#include <FL/Fl_Input.h>
#include <FL/Fl_Double_Window.h>

#include <vector>

class layoutable_c {

public:

  typedef enum {
    STRETCH_MINUS, STRETCH_PLUS, STRETCH_MIDDLE, STRETCH, STRETCH_SQUARE
  } stretch;

private:

  /* these are the grid positions of the widgets for the layouter */
  unsigned int _gridX, _gridY, _gridW, _gridH;

  /* what to do with the widget, if the space is bigger than the widget
  * requested: _MINUS put it on the upper or left edge, _PLUS put it on
  * the lower or right edge, _MIDDLE put it in the middle, STRETCH stretch
  * it to the available size
  * stratch _SQUARE takes the smaller size in the two dimensions and applies
  * this size to both (width and height)
  */
  stretch _stretchX, _stretchY;

  /* surroung the widget by pitch pixel */
  unsigned char _pitch;

  /* if space is left the available pixels will be distributed weighted by
   * these values
   */
  unsigned char _weightX, _weightY;

  /* if you want to enforce a size for the widget independent of the size
   * the widget itself requests, use these 2 values
   */
  int _minWidth, _minHeight;

public:

  layoutable_c(int gx = 0, int gy = 0, int gw = 1, int gh = 1) : _gridX(gx), _gridY(gy), _gridW(gw), _gridH(gh), _stretchX(STRETCH), _stretchY(STRETCH), _pitch(0), _weightX(0), _weightY(0), _minWidth(0), _minHeight(0) {}

  int getMinWidth(void) const { return _minWidth; }
  int getMinHeight(void) const { return _minHeight; }
  unsigned char getPitch(void) const { return _pitch; }
  unsigned char getWeightX(void) const { return _weightX; }
  unsigned char getWeightY(void) const { return _weightY; }
  stretch getStretchX(void) const { return _stretchX; }
  stretch getStretchY(void) const { return _stretchY; }

  /* the following function is for the layouter to find out how much space
   * the widget requires at least
   */
  virtual void getMinSize(int *width, int *height) const = 0;

  /* called by the layouter to find out the position
   * of the widget
   */
  void getGridValues(unsigned int *x, unsigned int *y, unsigned int *w, unsigned int *h) const {
    *x = _gridX;
    *y = _gridY;
    *w = _gridW;
    *h = _gridH;
  }

public:

  /* sets the layouting position for this widget. This is necessary
   * for e.g. automatic layout in tables. The caller creates unlayouted
   * widgets that get placed by the container
   */
  void setGridValues(unsigned int x, unsigned int y, unsigned int w, unsigned int h) {
    _gridX = x;
    _gridY = y;
    _gridW = w;
    _gridH = h;
  }

  /* sets the layouting parameter, this function will
   * automatically enable layouting for this widget, because
   * without layouting these values do not make any sense
   */
  void setLayoutParameter(stretch stretchX, stretch stretchY, unsigned char pitch, unsigned char weightX, unsigned char weightY) {
    _stretchX = stretchX;
    _stretchY = stretchY;
    _pitch = pitch;
    _weightX = weightX;
    _weightY = weightY;
  }

  void pitch(unsigned char pitch) { _pitch = pitch; }
  void stretchLeft(void) { _stretchX = STRETCH_MINUS; }
  void stretchRight(void) { _stretchX = STRETCH_PLUS; }
  void stretchTop(void) { _stretchY = STRETCH_MINUS; }
  void stretchBottom(void) { _stretchY = STRETCH_PLUS; }
  void stretchSquare(void) { _stretchX = _stretchY = STRETCH_SQUARE; }
  void stretchHCenter(void) { _stretchX = STRETCH_MIDDLE; }
  void stretchVCenter(void) { _stretchY = STRETCH_MIDDLE; }
  void stretchCenter(void) { _stretchX = _stretchY = STRETCH_MIDDLE; }
  void weight(unsigned char x, unsigned char y) { _weightX = x; _weightY = y; }

  /* sets the minimum size a widget should have */
  void setMinimumSize(unsigned int minWidth, unsigned int minHeight) {
    _minWidth = minWidth;
    _minHeight = minHeight;
  }
};

class layouter_c : public Fl_Group, public layoutable_c {

  public:

  virtual void getMinSize(int *width, int *height) const;

  void calcLayout(int task, std::vector<int> *widths, std::vector<int> *heights,
                  std::vector<int> *widgetsW, std::vector<int> *widgetsH, int targetW = 0, int tagetH = 0) const;

  virtual void resize(int x, int y, int w, int h);

  layouter_c(int x = 0, int y = 0, int w = 1, int h = 1) : Fl_Group(0, 0, 0, 0), layoutable_c(x, y, w, h) {}
};

/* now some basic widgets made layoutable */

class LFl_Box : public Fl_Box, public layoutable_c {

  public:

  LFl_Box(const char *txt, int x = 0, int y = 0, int w = 1, int h = 1) : Fl_Box(0, 0, 0, 0, txt), layoutable_c(x, y, w, h) {}
  LFl_Box(int x = 0, int y = 0, int w = 1, int h = 1) : Fl_Box(0, 0, 0, 0), layoutable_c(x, y, w, h) {}

  virtual void getMinSize(int *width, int *height) const {
    *width = 0;
    fl_font(labelfont(), labelsize());
    fl_measure(label(), *width, *height);
  }
};

class LFl_Frame : public layouter_c {

  layouter_c *l;

  public:

  LFl_Frame(int x = 0, int y = 0, int w = 1, int h = 1) : layouter_c(x, y, w, h) {
    l = new layouter_c();
    l->pitch(5);
    box(FL_ENGRAVED_BOX);
    pitch(1);
  }

  ~LFl_Frame() {
    delete l;
  }
};

class LFl_Button : public Fl_Button, public layoutable_c {

  public:

  LFl_Button(const char * text, int x = 0, int y = 0, int w = 1, int h = 1) : Fl_Button(0, 0, 0, 0, text), layoutable_c(x, y, w, h) {}

  virtual void getMinSize(int *width, int *height) const {
    *width = 0;
    fl_font(labelfont(), labelsize());
    fl_measure(label(), *width, *height);
    *width += 20;
    *height += 10;
  }
};

class LFl_Check_Button : public Fl_Check_Button, public layoutable_c {

  public:

  LFl_Check_Button(const char * text, int x = 0, int y = 0, int w = 1, int h = 1) : Fl_Check_Button(0, 0, 0, 0, text), layoutable_c(x, y, w, h) {
    stretchVCenter();
  }

  virtual void getMinSize(int *width, int *height) const {
    *width = 0;
    fl_font(labelfont(), labelsize());
    fl_measure(label(), *width, *height);
    *width += 18;
    *height += 4;
  }
};

class LFl_Round_Button : public Fl_Round_Button, public layoutable_c {

  public:

  LFl_Round_Button(const char *text, int x = 0, int y = 0, int w = 1, int h = 1) : Fl_Round_Button(0, 0, 0, 0, text), layoutable_c(x, y, w, h) {
    stretchVCenter();
  }

  virtual void getMinSize(int *width, int *height) const {
    *width = 0;
    fl_font(labelfont(), labelsize());
    fl_measure(label(), *width, *height);
    *width += 18;
    *height += 4;
  }
};

class LFl_Slider : public Fl_Slider, public layoutable_c {

  public:

  LFl_Slider(int x = 0, int y = 0, int w = 1, int h = 1) : Fl_Slider(0, 0, 0, 0), layoutable_c(x, y, w, h) {
  }

  virtual void getMinSize(int *width, int *height) const {
    switch (type()) {
    case FL_VERTICAL:
    case FL_VERT_FILL_SLIDER:
    case FL_VERT_NICE_SLIDER:
      *width = 15;
      *height = 30;
      break;
    case FL_HORIZONTAL:
    case FL_HOR_FILL_SLIDER:
    case FL_HOR_NICE_SLIDER:
      *width = 30;
      *height = 15;
      break;
    }
  }
};

class LFl_Input : public Fl_Input, public layoutable_c {

  public:

  LFl_Input(int x = 0, int y = 0, int w = 1, int h = 1) : Fl_Input(0, 0, 0, 0), layoutable_c(x, y, w, h) {
    stretchVCenter();
  }

  virtual void getMinSize(int *width, int *height) const {
    *width = 30;
    *height = 20;
  }

  // sets width so that the given text will fit into the input line
  void setMinWidth(const char *) {
  }
};

class LFl_Double_Window : public Fl_Double_Window {

  layouter_c * lay;

  public:

  LFl_Double_Window(void) : Fl_Double_Window(10, 10) {
    lay = new layouter_c();
    lay->resize(0, 0, 10, 10);
    resizable(lay);
  }

  void show(void) {
    int wy, hy;
    lay->getMinSize(&wy, &hy);

    size_range(wy, hy, wy, hy);

    size(wy, hy);

    Fl_Double_Window::show();
  }
};

#endif