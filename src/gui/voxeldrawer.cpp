/* Burr Solver
 * Copyright (C) 2003-2006  Andreas R�ver
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
#include "voxeldrawer.h"

#include "pieceColor.h"

#include "../lib/voxel.h"
#include "../lib/puzzle.h"
#include "../lib/assembly.h"

#include <math.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

voxelDrawer_c::voxelDrawer_c(int x,int y,int w,int h) :
  VoxelView(x, y, w, h), markerType(-1),
  colors(pieceColor), _useLightning(true),
  _gtChanged(true)
{
};

void voxelDrawer_c::setTransformationMatrix(GLfloat m[16]) {
  for (unsigned int i = 0; i < 16; i++)
    transformMatrix[i] = m[i];
}

void voxelDrawer_c::gridTypeChanged(void) {
  GLfloat m[16] = {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1};

  setTransformationMatrix(m);

  _gtChanged = false;
  redraw();
}

void voxelDrawer_c::drawVoxelSpace() {

  glShadeModel(GL_FLAT);

  for (unsigned int run = 0; run < 2; run++) {
    for (unsigned int piece = 0; piece < shapes.size(); piece++) {

      if (shapes[piece].a == 0)
        continue;

      // in run 0 we only paint opaque objects and in run 1 only transparent ones
      // this lets the transparent objects be always in front of the others
      if ((run == 0) && (shapes[piece].a != 1)) continue;
      if ((run == 1) && (shapes[piece].a == 1)) continue;

      glPushMatrix();

      switch(trans) {
      case ScaleRotateTranslate:
        glTranslatef(shapes[piece].x - shapes[piece].shape->getHx(),
                     shapes[piece].y - shapes[piece].shape->getHy(),
                     shapes[piece].z - shapes[piece].shape->getHz());
        glScalef(shapes[piece].scale, shapes[piece].scale, shapes[piece].scale);
        addRotationTransformation();
        glMultMatrixf(transformMatrix);
        {
          float cx, cy, cz;
          calculateSize(shapes[piece].shape, &cx, &cy, &cz);
          glTranslatef(-0.5*cx, -0.5*cy, -0.5*cz);
        }
        break;
      case TranslateRoateScale:
        addRotationTransformation();
        glTranslatef(shapes[piece].x - shapes[piece].shape->getHx(),
                     shapes[piece].y - shapes[piece].shape->getHy(),
                     shapes[piece].z - shapes[piece].shape->getHz());
        glMultMatrixf(transformMatrix);
        {
          float cx, cy, cz;
          calculateSize(shapes[piece].shape, &cx, &cy, &cz);
          glTranslatef(-0.5*cx, -0.5*cy, -0.5*cz);
        }
        glScalef(shapes[piece].scale, shapes[piece].scale, shapes[piece].scale);
        break;
      case CenterTranslateRoateScale:
        addRotationTransformation();
        glMultMatrixf(transformMatrix);
        glTranslatef(shapes[piece].x - shapes[piece].shape->getHx(),
                     shapes[piece].y - shapes[piece].shape->getHy(),
                     shapes[piece].z - shapes[piece].shape->getHz());
        glTranslatef(-centerX, -centerY, -centerZ);
        glScalef(shapes[piece].scale, shapes[piece].scale, shapes[piece].scale);
        break;
      default:
        break;
      }

      if (_showCoordinateSystem) {
        if (_useLightning) glDisable(GL_LIGHTING);
        glDisable(GL_BLEND);
        glBegin(GL_LINES);
        float cx, cy, cz;
        calculateSize(shapes[piece].shape, &cx, &cy, &cz);

        glColor3f(1, 0,    0); glVertex3f(-1, -1, -1); glVertex3f(cx+1, -1, -1);
        glColor3f(0, 0.75, 0); glVertex3f(-1, -1, -1); glVertex3f(-1, cy+1, -1);
        glColor3f(0, 0,    1); glVertex3f(-1, -1, -1); glVertex3f(-1, -1, cz+1);
        glEnd();
        if (_useLightning) glEnable(GL_LIGHTING);
        glEnable(GL_BLEND);
      }

      const shapeInfo * shape = &shapes[piece];

      for (unsigned int x = 0; x < shape->shape->getX(); x++)
        for (unsigned int y = 0; y < shape->shape->getY(); y++)
          for (unsigned int z = 0; z < shape->shape->getZ(); z++) {

            if (shape->shape->isEmpty(x, y , z))
              continue;

            float cr, cg, cb, ca;
            cr = cg = cb = 0;
            ca = 1;

            switch (colors) {
              case pieceColor:
                if ((x+y+z) & 1) {
                  cr = lightPieceColor(shape->r);
                  cg = lightPieceColor(shape->g);
                  cb = lightPieceColor(shape->b);
                  ca = shape->a;
                } else {
                  cr = darkPieceColor(shape->r);
                  cg = darkPieceColor(shape->g);
                  cb = darkPieceColor(shape->b);
                  ca = shape->a;
                }
                break;
              case paletteColor:
                unsigned int color = shape->shape->getColor(x, y, z);
                if ((color == 0) || (color - 1 >= palette.size())) {
                  if ((x+y+z) & 1) {
                    cr = lightPieceColor(shape->r);
                    cg = lightPieceColor(shape->g);
                    cb = lightPieceColor(shape->b);
                    ca = shape->a;
                  } else {
                    cr = darkPieceColor(shape->r);
                    cg = darkPieceColor(shape->g);
                    cb = darkPieceColor(shape->b);
                    ca = shape->a;
                  }
                } else {
                  cr = palette[color-1].r;
                  cg = palette[color-1].g;
                  cb = palette[color-1].b;
                  ca = shape->a;
                }
            }

            if (shape->dim) {
              cr = 1 - (1 - cr) * 0.2;
              cg = 1 - (1 - cg) * 0.2;
              cb = 1 - (1 - cb) * 0.2;
            }

            glColor4f(cr, cg, cb, ca);

            switch (shape->mode) {
              case normal:
                if (shape->shape->getState(x, y , z) == voxel_c::VX_VARIABLE) {
                  drawNormalVoxel(shape->shape, x, y, z, shape->a, shape->dim ? 0 : 0.05);
                  glColor4f(0, 0, 0, shape->a);
                  drawVariableMarkers(shape->shape, x, y, z);
                } else
                  drawNormalVoxel(shape->shape, x, y, z, shape->a, shape->dim ? 0 : 0.05);
                break;
              case gridline:
                drawFrame(shape->shape, x, y, z, 0.05);
                break;
              case invisible:
                break;
            }
          }

      // the marker should be only active, when only one shape is there
      // otherwise it's drawn for every shape
      if ((markerType >= 0) && (mX1 <= mX2) && (mY1 <= mY2)) {

        if (_useLightning) glDisable(GL_LIGHTING);
        glDisable(GL_BLEND);

        glColor4f(0, 0, 0, 1);

        drawCursor(shapes[piece].shape->getX(), shapes[piece].shape->getY(), shapes[piece].shape->getZ());

        if (_useLightning) glEnable(GL_LIGHTING);
        glEnable(GL_BLEND);
      }

      glPopMatrix();
    }

    glDepthMask(GL_FALSE);
  }

  glDepthMask(GL_TRUE);
}

void voxelDrawer_c::drawData(void) {

  if (_gtChanged)
    gridTypeChanged();

  glShadeModel(GL_FLAT);
  glEnable(GL_DEPTH_TEST);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glLineWidth(3);

  if (_useLightning)
    glEnable(GL_LIGHTING);
  else
    glDisable(GL_LIGHTING);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();

  drawVoxelSpace();
  glPopMatrix();
}

unsigned int voxelDrawer_c::addSpace(const voxel_c * vx) {
  shapeInfo i;

  i.r = i.g = i.b = 1;
  i.a = 1;
  i.shape = vx;

  i.mode = normal;

  i.x = i.y = i.z = 0;
  i.scale = 1;

  i.dim = false;

  shapes.push_back(i);

  updateRequired();

  return shapes.size()-1;
}

void voxelDrawer_c::clearSpaces(void) {

  for (unsigned int i = 0; i < shapes.size(); i++)
    delete shapes[i].shape;

  shapes.clear();
  updateRequired();
}

unsigned int voxelDrawer_c::spaceNumber(void) {
  return shapes.size();
}

void voxelDrawer_c::setSpaceColor(unsigned int nr, float r, float g, float b, float a) {
  shapes[nr].r = r;
  shapes[nr].g = g;
  shapes[nr].b = b;
  shapes[nr].a = a;

  updateRequired();
}

void voxelDrawer_c::setSpaceColor(unsigned int nr, float a) {
  shapes[nr].a = a;

  updateRequired();
}

void voxelDrawer_c::setSpacePosition(unsigned int nr, float x, float y, float z, float scale) {
  shapes[nr].x = x;
  shapes[nr].y = y;
  shapes[nr].z = z;
  shapes[nr].scale = scale;

  updateRequired();
}

void voxelDrawer_c::setSpaceDim(unsigned int nr, bool dim) {

  shapes[nr].dim = dim;

  updateRequired();
}


void voxelDrawer_c::setDrawingMode(unsigned int nr, drawingMode mode) {
  shapes[nr].mode = mode;

  updateRequired();
}

void voxelDrawer_c::setColorMode(colorMode color) {
  colors = color;

  updateRequired();
}

void voxelDrawer_c::setTransformationType(transformationType type) {
  trans = type;

  updateRequired();
}

void voxelDrawer_c::addPaletteEntry(float r, float g, float b) {

  colorInfo ci;

  ci.r = r;
  ci.g = g;
  ci.b = b;

  palette.push_back(ci);
}

void voxelDrawer_c::setMarker(int x1, int y1, int x2, int y2, int z, int mT) {
  markerType = mT;
  mX1 = x1;
  mY1 = y1;
  mX2 = x2;
  mY2 = y2;
  mZ = z;
}

void voxelDrawer_c::hideMarker(void) {
  markerType = -1;
}

void voxelDrawer_c::showSingleShape(const puzzle_c * puz, unsigned int shapeNum) {

  hideMarker();
  clearSpaces();
  unsigned int num = addSpace(puz->getGridType()->getVoxel(puz->getShape(shapeNum)));

  setSpaceColor(num, pieceColorR(shapeNum), pieceColorG(shapeNum), pieceColorB(shapeNum), 1);

  setTransformationType(TranslateRoateScale);
  showCoordinateSystem(true);
}

void voxelDrawer_c::showProblem(const puzzle_c * puz, unsigned int probNum, unsigned int selShape) {

  hideMarker();
  clearSpaces();

  if (probNum < puz->problemNumber()) {

    // first find out how to arrange the pieces:
    unsigned int square = 3;
    while (square * (square-2) < puz->probShapeNumber(probNum)) square++;

    unsigned int num;

    float diagonal = 0;

    // now find a scaling factor, so that all pieces fit into their square
    if (puz->probGetResult(probNum) < puz->shapeNumber()) {

      if (puz->probGetResultShape(probNum)->getDiagonal() > diagonal)
        diagonal = puz->probGetResultShape(probNum)->getDiagonal();
    }

    // check the selected shape
    if (selShape < puz->shapeNumber()) {

      if (puz->getShape(selShape)->getDiagonal() > diagonal)
        diagonal = puz->getShape(selShape)->getDiagonal();
    }

    for (unsigned int p = 0; p < puz->probShapeNumber(probNum); p++)
      if (puz->probGetShapeShape(probNum, p)->getDiagonal() > diagonal)
        diagonal = puz->probGetShapeShape(probNum, p)->getDiagonal();

    diagonal = sqrt(diagonal)/1.5;

    // now place the result shape
    if (puz->probGetResult(probNum) < puz->shapeNumber()) {

      num = addSpace(puz->getGridType()->getVoxel(puz->probGetResultShape(probNum)));
      setSpaceColor(num,
                            pieceColorR(puz->probGetResult(probNum)),
                            pieceColorG(puz->probGetResult(probNum)),
                            pieceColorB(puz->probGetResult(probNum)), 1);
      setSpacePosition(num,
                               0.5* (square*diagonal) * (1.0/square - 0.5),
                               0.5* (square*diagonal) * (0.5 - 1.0/square), -20, 1.0);
    }

    // now place the selected shape
    if (selShape < puz->shapeNumber()) {

      num = addSpace(puz->getGridType()->getVoxel(puz->getShape(selShape)));
      setSpaceColor(num,
                            pieceColorR(selShape),
                            pieceColorG(selShape),
                            pieceColorB(selShape), 1);
      setSpacePosition(num,
                               0.5* (square*diagonal) * (0.5 - 0.5/square),
                               0.5* (square*diagonal) * (0.5 - 0.5/square), -20, 0.5);
    }

    // and now the shapes
    int unsigned line = 2;
    int unsigned col = 0;
    for (unsigned int p = 0; p < puz->probShapeNumber(probNum); p++) {
      num = addSpace(puz->getGridType()->getVoxel(puz->probGetShapeShape(probNum, p)));

      setSpaceColor(num,
                            pieceColorR(puz->probGetShape(probNum, p)),
                            pieceColorG(puz->probGetShape(probNum, p)),
                            pieceColorB(puz->probGetShape(probNum, p)), 1);

      setSpacePosition(num,
                               0.5* (square*diagonal) * ((col+0.5)/square - 0.5),
                               0.5* (square*diagonal) * (0.5 - (line+0.5)/square),
                               -20, 0.5);

      col++;
      if (col == square) {
        col = 0;
        line++;
      }
    }

    setTransformationType(ScaleRotateTranslate);
    showCoordinateSystem(false);
  }
}

void voxelDrawer_c::showColors(const puzzle_c * puz, bool show) {

  if (show) {

    clearPalette();
    for (unsigned int i = 0; i < puz->colorNumber(); i++) {
      unsigned char r, g, b;
      puz->getColor(i, &r, &g, &b);
      addPaletteEntry(r/255.0, g/255.0, b/255.0);
    }
    setColorMode(paletteColor);

  } else
    setColorMode(pieceColor);

}

void voxelDrawer_c::showAssembly(const puzzle_c * puz, unsigned int probNum, unsigned int solNum) {

  hideMarker();
  clearSpaces();

  if ((probNum < puz->problemNumber()) &&
      (solNum < puz->probSolutionNumber(probNum))) {

    unsigned int num;

    const assembly_c * assm = puz->probGetAssembly(probNum, solNum);

    unsigned int piece = 0;

    // and now the shapes
    for (unsigned int p = 0; p < puz->probShapeNumber(probNum); p++)
      for (unsigned int q = 0; q < puz->probGetShapeCount(probNum, p); q++) {

        num = addSpace(puz->getGridType()->getVoxel(puz->probGetShapeShape(probNum, p), assm->getTransformation(piece)));

        setSpacePosition(num, assm->getX(piece), assm->getY(piece), assm->getZ(piece), 1);

        setSpaceColor(num,
                              pieceColorR(puz->probGetShape(probNum, p), q),
                              pieceColorG(puz->probGetShape(probNum, p), q),
                              pieceColorB(puz->probGetShape(probNum, p), q), 1);

        piece++;
      }

    float cx, cy, cz;
    calculateSize(puz->probGetResultShape(probNum), &cx, &cy, &cz);
    setCenter(cx*0.5, cy*0.5, cz*0.5);
    setTransformationType(CenterTranslateRoateScale);
    showCoordinateSystem(false);
  }
}

void voxelDrawer_c::showAssemblerState(const puzzle_c * puz, unsigned int probNum, const assembly_c * assm) {

  hideMarker();
  clearSpaces();

  if (probNum < puz->problemNumber()) {

    unsigned int num;
    unsigned int piece = 0;

    // and now the shapes
    for (unsigned int p = 0; p < puz->probShapeNumber(probNum); p++)
      for (unsigned int q = 0; q < puz->probGetShapeCount(probNum, p); q++) {

        if (assm->getTransformation(piece) < 0xff) {

          num = addSpace(puz->getGridType()->getVoxel(puz->probGetShapeShape(probNum, p), assm->getTransformation(piece)));

          setSpacePosition(num, assm->getX(piece), assm->getY(piece), assm->getZ(piece), 1);

          setSpaceColor(num,
              pieceColorR(puz->probGetShape(probNum, p), q),
              pieceColorG(puz->probGetShape(probNum, p), q),
              pieceColorB(puz->probGetShape(probNum, p), q), 1);
        }

        piece++;
      }

    setCenter(0.5*puz->probGetResultShape(probNum)->getX(),
                      0.5*puz->probGetResultShape(probNum)->getY(),
                      0.5*puz->probGetResultShape(probNum)->getZ()
                     );
    setTransformationType(CenterTranslateRoateScale);
    showCoordinateSystem(false);
  }
}

void voxelDrawer_c::showPlacement(const puzzle_c * puz, unsigned int probNum, unsigned int piece, unsigned char trans, int x, int y, int z) {

  clearSpaces();
  hideMarker();
  setTransformationType(CenterTranslateRoateScale);
  showCoordinateSystem(false);
  setCenter(0.5*puz->probGetResultShape(probNum)->getX(),
                    0.5*puz->probGetResultShape(probNum)->getY(),
                    0.5*puz->probGetResultShape(probNum)->getZ()
                   );

  int num;

  if (trans < puz->getGridType()->getSymmetries()->getNumTransformationsMirror()) {

    int shape = 0;
    unsigned int p = piece;
    while (p >= puz->probGetShapeCount(probNum, shape)) {
      p -= puz->probGetShapeCount(probNum, shape);
      shape++;
    }

    num = addSpace(puz->getGridType()->getVoxel(puz->probGetShapeShape(probNum, shape), trans));
    setSpacePosition(num, x, y, z, 1);
    setSpaceColor(num,
                          pieceColorR(puz->probGetShape(probNum, shape), p),
                          pieceColorG(puz->probGetShape(probNum, shape), p),
                          pieceColorB(puz->probGetShape(probNum, shape), p), 1);
    setDrawingMode(num, normal);
  }

  num = addSpace(puz->getGridType()->getVoxel(puz->probGetResultShape(probNum)));
  setSpaceColor(num,
                        pieceColorR(puz->probGetResult(probNum)),
                        pieceColorG(puz->probGetResult(probNum)),
                        pieceColorB(puz->probGetResult(probNum)), 1);
  setDrawingMode(num, gridline);
}


void voxelDrawer_c::updatePositions(PiecePositions *shifting) {

  for (unsigned int p = 0; p < spaceNumber(); p++) {
    setSpacePosition(p, shifting->getX(p), shifting->getY(p), shifting->getZ(p), 1);
    setSpaceColor(p, shifting->getA(p));
  }
}

void voxelDrawer_c::dimStaticPieces(PiecePositions *shifting) {

  for (unsigned int p = 0; p < spaceNumber(); p++) {
    setSpaceDim(p, !shifting->moving(p));
  }
}

void voxelDrawer_c::updateVisibility(PieceVisibility * pcvis) {

  for (unsigned int p = 0; p < spaceNumber(); p++) {

    switch(pcvis->getVisibility(p)) {
    case 0:
      setDrawingMode(p, normal);
      break;
    case 1:
      setDrawingMode(p, gridline);
      break;
    case 2:
      setDrawingMode(p, invisible);
      break;
    }
  }
}

void voxelDrawer_c::drawGridRect(double x0, double y0, double z0,
                     double v1x, double v1y, double v1z,
                     double v2x, double v2y, double v2z, int diag) {

  glBegin(GL_LINES);
  glVertex3f(x0, y0, z0); glVertex3f(x0+v1x, y0+v1y, z0+v1z);
  glVertex3f(x0+v1x, y0+v1y, z0+v1z); glVertex3f(x0+v1x+v2x, y0+v1y+v2y, z0+v1z+v2z);
  glVertex3f(x0+v1x+v2x, y0+v1y+v2y, z0+v1z+v2z); glVertex3f(x0+v2x, y0+v2y, z0+v2z);
  glVertex3f(x0+v2x, y0+v2y, z0+v2z); glVertex3f(x0, y0, z0);

  int state1 = 0;
  int state2 = 0;


  float x1 = x0 + v1x;
  float y1 = y0 + v1y;
  float z1 = z0 + v1z;

  float x2 = x0 + v1x;
  float y2 = y0 + v1y;
  float z2 = z0 + v1z;

  float xe = x0 + v2x;
  float ye = y0 + v2y;
  float ze = z0 + v2z;

  while ((fabs(x1 - xe) > 0.01) || (fabs(y1 - ye) > 0.01) || (fabs(z1 - ze) > 0.01)) {
    // v1=(x1, y1, z1) first goes along vector1 =(v1x, v1y, v1z) and then along vector2

    if (state1 == 0) {
      x1 -= v1x/diag;
      y1 -= v1y/diag;
      z1 -= v1z/diag;

      if ((v1x) && (fabs(x1 - x0) < 0.01) ||
          (v1y) && (fabs(y1 - y0) < 0.01) ||
          (v1z) && (fabs(z1 - z0) < 0.01)) {
        state1 = 1;
      }
    } else {

      x1 += v2x/diag;
      y1 += v2y/diag;
      z1 += v2z/diag;
    }

    if (state2 == 0) {
      x2 += v2x/diag;
      y2 += v2y/diag;
      z2 += v2z/diag;

      if ((v2x) && (fabs(x2 - (x0+v2x+v1x)) < 0.01) ||
          (v2y) && (fabs(y2 - (y0+v2y+v1y)) < 0.01) ||
          (v2z) && (fabs(z2 - (z0+v2z+v1z)) < 0.01)) {
        state2 = 1;
      }
    } else {

      x2 -= v1x/diag;
      y2 -= v1y/diag;
      z2 -= v1z/diag;
    }

    glVertex3f(x1, y1, z1); glVertex3f(x2, y2, z2);
  }

  glEnd();
}


void voxelDrawer_c::drawGridTriangle(double x0, double y0, double z0,
                         double v1x, double v1y, double v1z,
                         double v2x, double v2y, double v2z, int diag) {

  glBegin(GL_LINES);
  glVertex3f(x0, y0, z0); glVertex3f(x0+v1x, y0+v1y, z0+v1z);
  glVertex3f(x0+v1x, y0+v1y, z0+v1z); glVertex3f(x0+v2x, y0+v2y, z0+v2z);
  glVertex3f(x0+v2x, y0+v2y, z0+v2z); glVertex3f(x0, y0, z0);

  float x1 = x0;
  float y1 = y0;
  float z1 = z0;

  float x2 = x0;
  float y2 = y0;
  float z2 = z0;

  float xe = x0 + v1x;
  float ye = y0 + v1y;
  float ze = z0 + v1z;

  while ((fabs(x1 - xe) > 0.01) || (fabs(y1 - ye) > 0.01) || (fabs(z1 - ze) > 0.01)) {

    x1 += v1x/diag;
    y1 += v1y/diag;
    z1 += v1z/diag;

    x2 += v2x/diag;
    y2 += v2y/diag;
    z2 += v2z/diag;

    glVertex3f(x1, y1, z1); glVertex3f(x2, y2, z2);
  }

  glEnd();
}

// this function finds out if a given square is inside the selected region
// this check includes the symmetric and comulmn edit modes
bool voxelDrawer_c::inRegion(int x, int y, int z, int x1, int x2, int y1, int y2, int z1, int z2, int sx, int sy, int sz, int mode) {

  if ((x < 0) || (y < 0) || (z < 0) || (x >= sx) || (y >= sy) || (z >= sz)) return false;

  if (mode == 0)
    return (x1 <= x) && (x <= x2) && (y1 <= y) && (y <= y2) && (z1 <= z) && (z <= z2);
  if (mode == voxelDrawer_c::TOOL_STACK_Y)
    return (x1 <= x) && (x <= x2) && (z1 <= z) && (z <= z2);
  if (mode == voxelDrawer_c::TOOL_STACK_X)
    return (y1 <= y) && (y <= y2) && (z1 <= z) && (z <= z2);
  if (mode == voxelDrawer_c::TOOL_STACK_Z)
    return (x1 <= x) && (x <= x2) && (y1 <= y) && (y <= y2);

  if (mode == voxelDrawer_c::TOOL_STACK_X + voxelDrawer_c::TOOL_STACK_Y)
    return ((x1 <= x) && (x <= x2) || (y1 <= y) && (y <= y2)) && ((z1 <= z) && (z <= z2));
  if (mode == voxelDrawer_c::TOOL_STACK_X + voxelDrawer_c::TOOL_STACK_Z)
    return ((x1 <= x) && (x <= x2) || (z1 <= z) && (z <= z2)) && ((y1 <= y) && (y <= y2));
  if (mode == voxelDrawer_c::TOOL_STACK_Y + voxelDrawer_c::TOOL_STACK_Z)
    return ((y1 <= y) && (y <= y2) || (y1 <= y) && (y <= y2)) && ((x1 <= x) && (x <= x2));

  if (mode == voxelDrawer_c::TOOL_STACK_X + voxelDrawer_c::TOOL_STACK_Y + voxelDrawer_c::TOOL_STACK_Z)
    return ((x1 <= x) && (x <= x2) && (y1 <= y) && (y <= y2) ||
        (x1 <= x) && (x <= x2) && (z1 <= z) && (z <= z2) ||
        (y1 <= y) && (y <= y2) && (z1 <= z) && (z <= z2));

  if (mode & voxelDrawer_c::TOOL_MIRROR_X)
    return inRegion(x, y, z, x1, x2, y1, y2, z1, z2, sx, sy, sz, mode & ~voxelDrawer_c::TOOL_MIRROR_X) ||
      inRegion(sx-x-1, y, z, x1, x2, y1, y2, z1, z2, sx, sy, sz, mode & ~voxelDrawer_c::TOOL_MIRROR_X);

  if (mode & voxelDrawer_c::TOOL_MIRROR_Y)
    return inRegion(x, y, z, x1, x2, y1, y2, z1, z2, sx, sy, sz, mode & ~voxelDrawer_c::TOOL_MIRROR_Y) ||
      inRegion(x, sy-y-1, z, x1, x2, y1, y2, z1, z2, sx, sy, sz, mode & ~voxelDrawer_c::TOOL_MIRROR_Y);

  if (mode & voxelDrawer_c::TOOL_MIRROR_Z)
    return inRegion(x, y, z, x1, x2, y1, y2, z1, z2, sx, sy, sz, mode & ~voxelDrawer_c::TOOL_MIRROR_Z) ||
      inRegion(x, y, sz-z-1, x1, x2, y1, y2, z1, z2, sx, sy, sz, mode & ~voxelDrawer_c::TOOL_MIRROR_Z);

  return false;
}

