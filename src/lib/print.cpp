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
#include "print.h"

#include "voxel.h"
#include "puzzle.h"
#include "disassembly.h"
#include "assembly.h"

void print(const voxel_c * v, char base) {
  for (unsigned int z = 0; z < v->getZ(); z++) {
    printf(" +");
    for (unsigned int x = 0; x < v->getX(); x++)
      printf("-");
    printf("+");
  }
  printf(" bx %i-%i by %i-%i bz %i-%i h: %i %i %i \n", v->boundX1(), v->boundX2(), v->boundY1(), v->boundY2(), v->boundZ1(), v->boundZ2(), v->getHx(), v->getHy(), v->getHz());

  for (unsigned int y = 0; y < v->getY(); y++) {
    for (unsigned int z = 0; z < v->getZ(); z++) {
      printf(" !");
      for (unsigned int x = 0; x < v->getX(); x++)
        if (v->get(x, y, z) != 0)
          printf("%c", base + v->get(x, y, z)-1);
        else
          printf(" ");
      printf("!");
    }
    printf("\n");
  }

  { for (unsigned int z = 0; z < v->getZ(); z++) {
      printf(" +");
      for (unsigned int x = 0; x < v->getX(); x++)
        printf("-");
      printf("+");
    }
  }
  printf("\n");
}

void print(const puzzle_c * p) {

  for (unsigned int s = 0; s < p->shapeNumber(); s++) {
    printf("shape %i:\n", s);
    print(p->getShape(s));
  }

  printf("=======================================================\n");

  for (unsigned int pr = 0; pr < p->problemNumber(); pr++) {
    printf("problem %i (%s):\n", pr, p->probGetName(pr).c_str());
    printf(" result shape: %i\n", p->probGetResult(pr));

    for (unsigned int sh = 0; sh < p->probShapeNumber(pr); sh++)
      if (p->probGetShapeCount(pr, sh) > 1)
        printf(" piece shape: %i times shape number %i\n", p->probGetShapeCount(pr, sh), p->probGetShape(pr, sh));
      else
        printf(" piece shape: %i\n", p->probGetShape(pr, sh));

    printf("-------------------------------------------------------\n");
  }
}

static void print_rec(const separation_c * s, voxel_c ** pieces, int sx, int sy, int sz, unsigned int * pieceNum) {

  for (unsigned int i = 0; i <= s->getMoves(); i++) {

    const state_c * st = s->getState(i);

    for (int z = -sz/2; z < sz+sz/2; z++) {
      printf(" +");
      for (int x = -sx/2; x < sx+sx/2; x++)
        printf("-");
      printf("+");
    }
    printf("\n");

    for (int y = -sy/2; y < sy+sy/2; y++) {
      for (int z = -sz/2; z < sz+sz/2; z++) {
        printf(" !");
        for (int x = -sx/2; x < sx+sx/2; x++) {
          char c = ' ';

          for (unsigned int pc = 0; pc < s->getPieceNumber(); pc++)
            if (pieces[pc]->isFilled2(x - (st->getX(pc) - pieces[pc]->getHx()),
                                      y - (st->getY(pc) - pieces[pc]->getHy()),
                                      z - (st->getZ(pc) - pieces[pc]->getHz()))) {
              c = 'a' + pieceNum[pc];
              break;
            }

          printf("%c", c);
        }
        printf("!");
      }
      printf("\n");
    }

    { for (int z = -sz/2; z < sz+sz/2; z++) {
        printf(" +");
        for (int x = -sx/2; x < sx+sx/2; x++)
          printf("-");
        printf("+");
      }
    }
    printf("\n");
  }

  printf(" puzzle separates into  2 parts\n");

  if (s->getRemoved()) {

    unsigned int * pieceNum2 = new unsigned int [s->getPieceNumber()];

    int pos = 0;
    for (unsigned int i = 0; i < s->getPieceNumber(); i++)
      if (s->getState(s->getMoves())->pieceRemoved(i))
        pieceNum2[pos++] = pieceNum[i];

    print_rec(s->getRemoved(), pieces, sx, sy, sz, pieceNum2);

    delete [] pieceNum2;
  }
  if (s->getLeft()) {

    unsigned int * pieceNum2 = new unsigned int [s->getPieceNumber()];

    int pos = 0;
    for (unsigned int i = 0; i < s->getPieceNumber(); i++)
      if (!(s->getState(s->getMoves())->pieceRemoved(i)))
        pieceNum2[pos++] = pieceNum[i];

    print_rec(s->getLeft(), pieces, sx, sy, sz, pieceNum2);

    delete [] pieceNum2;
  }
}

void print(const separation_c * s, const assembly_c * a, const puzzle_c * p, unsigned int prob) {

  const voxel_c * res = p->probGetResultShape(prob);

  voxel_c ** pieces = new voxel_c*[a->placementCount()];

  unsigned int pc = 0;

  for (unsigned int i = 0; i < p->probShapeNumber(prob); i++)
    for (unsigned int j = 0; j < p->probGetShapeCount(prob, i); j++) {

      pieces[pc] = p->getGridType()->getVoxel(p->probGetShapeShape(prob, i));
      bt_assert(pieces[pc]->transform(a->getTransformation(pc)));
      pc++;
    }

  unsigned int * pieceNum = new unsigned int [a->placementCount()];
  for (unsigned int i = 0; i < a->placementCount(); i++)
    pieceNum[i] = i;

  print_rec(s, pieces, res->getX(), res->getY(), res->getZ(), pieceNum);

  for (unsigned int pc = 0; pc < a->placementCount(); pc++)
    delete pieces[pc];

  delete [] pieceNum;
  delete [] pieces;
}

void print(const assembly_c * a, const puzzle_c * p, unsigned int prob) {

  const voxel_c * res = p->probGetResultShape(prob);


  voxel_c ** pieces = new voxel_c*[a->placementCount()];

  unsigned int pc = 0;

  for (unsigned int i = 0; i < p->probShapeNumber(prob); i++)
    for (unsigned int j = 0; j < p->probGetShapeCount(prob, i); j++) {

      pieces[pc] = p->getGridType()->getVoxel(p->probGetShapeShape(prob, i));
      bt_assert(pieces[pc]->transform(a->getTransformation(pc)));
      pc++;
    }

  for (unsigned int z = 0; z < res->getZ(); z++) {
    printf(" +");
    for (unsigned int x = 0; x < res->getX(); x++)
      printf("-");
    printf("+");
  }
  printf("\n");

  for (unsigned int y = 0; y < res->getY(); y++) {
    for (unsigned int z = 0; z < res->getZ(); z++) {
      printf(" !");
      for (unsigned int x = 0; x < res->getX(); x++) {
        char c = ' ';

        for (unsigned int pc = 0; pc < a->placementCount(); pc++)
          if (pieces[pc]->isFilled2(x - (a->getX(pc) - pieces[pc]->getHx()),
                                    y - (a->getY(pc) - pieces[pc]->getHy()),
                                    z - (a->getZ(pc) - pieces[pc]->getHz()))) {
            c = 'a' + pc;
            break;
          }

        printf("%c", c);
      }
      printf("!");
    }
    printf("\n");
  }

  { for (unsigned int z = 0; z < res->getZ(); z++) {
      printf(" +");
      for (unsigned int x = 0; x < res->getX(); x++)
        printf("-");
      printf("+");
    }
  }
  printf("\n");

  for (unsigned int pc = 0; pc < a->placementCount(); pc++)
    delete pieces[pc];

  delete [] pieces;
}
