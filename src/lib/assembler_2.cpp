/* BurrTools
 *
 * BurrTools is the legal property of its developers, whose
 * names are listed in the COPYRIGHT file, which is included
 * within the source distribution.
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
#include "assembler_2.h"

#include "bt_assert.h"
#include "problem.h"
#include "puzzle.h"
#include "voxel.h"
#include "assembly.h"
#include "gridtype.h"
#include "cadical/src/cadical.hpp"

#include "../tools/xml.h"

#include <cstdlib>
#include <cstring>
#include <memory>

#include "../config.h"

#ifdef WIN32
#define snprintf _snprintf
#endif

#define ASSEMBLER_VERSION "3.0"


assembler_2_c::assembler_2_c(const problem_c & prob) :
  assembler_c(),
  problem(prob),
  solver(0),
  reducePiece(0),
  avoidTransformedAssemblies(0), avoidTransformedMirror(0)
{
}

assembler_2_c::~assembler_2_c() {
  if (avoidTransformedMirror) delete avoidTransformedMirror;
  if (solver) delete solver;
}

/* add a piece to the cache, but only if it is not already there. If it is added return the
 * piece pointer otherwise return null
 */
static voxel_c * addToCache(voxel_c * cache[], unsigned int * fill, voxel_c * piece) {

  for (unsigned int i = 0; i < *fill; i++)
    if (cache[i]->identicalInBB(piece)) {
      delete piece;
      return 0;
    }

  cache[*fill] = piece;
  (*fill)++;
  return piece;
}


bool assembler_2_c::canPlace(const voxel_c * piece, int x, int y, int z) const {

  if (!piece->onGrid(x, y, z))
    return false;

  const voxel_c * result = getResultShape(problem);

  for (unsigned int pz = piece->boundZ1(); pz <= piece->boundZ2(); pz++)
    for (unsigned int py = piece->boundY1(); py <= piece->boundY2(); py++)
      for (unsigned int px = piece->boundX1(); px <= piece->boundX2(); px++)
        if (
            // the piece can not be place if the result is empty and the piece is filled at a given voxel
            ((piece->getState(px, py, pz) == voxel_c::VX_FILLED) &&
             (result->getState(x+px, y+py, z+pz) == voxel_c::VX_EMPTY)) ||

            // the piece can also not be placed when the colour constraints don't fit
            !problem.placementAllowed(piece->getColor(px, py, pz), result->getColor(x+px, y+py, z+pz))

           )
          return false;

  return true;
}


/**
 * this function prepares the matrix of nodes for the recursive function
 * I've done some additions to Knuths algorithm to implement variable
 * voxels (empty spaces in the solution) and multiple instances of the
 * same piece. Empty voxels in the result are done by removing columns
 * from the matrix. This will prevent the algorithm from filling the
 * corresponding voxels. But we need to have the constraints that these
 * columns place on the solution. This is done by adding these columns
 * to the matrix but behind the normal columns. These additional columns
 * wont be searched by the alg. if it looks for the next task to achieve.
 *
 * Multiple instances of the same piece is handles in a similar way. To
 * prevent finding the same solution again and again with just the
 * pieces swapping places we number the pieces and their possible
 * placements and disallow that the position number of piece n is lower
 * than the position number of piece n-1. This can be achieved by adding
 * more constraint columns. There need to be one column for each
 *
 * negative result show there is something wrong: the place -result has not
 * possible position inside the result
 */
int assembler_2_c::prepare(void) {

  const voxel_c * result = getResultShape(problem);

  if (solver) delete solver;
  solver = new CaDiCaL::Solver;
  solverState = CaDiCaL::READY;

  /* find the symmetry breaker
   *
   * OK, what idea is behind this: we try to find as few double solutions as possible
   * because we don't want to fist search them and later on discard them because they are
   * double, so what do we do to prevent double solutions?
   *
   * Select one piece and remove rotations from this piece so that we don't even try to
   * place this piece in all possible positions. But which rotations need to be removed?
   * This depends on the symmetries that are present in the result and the symmetries
   * that are present in the piece
   */
  symmetries_t resultSym = result->selfSymmetries();
  const gridType_c * gt = problem.getPuzzle().getGridType();
  const symmetries_c * sym = problem.getPuzzle().getGridType()->getSymmetries();
  unsigned int symBreakerShape = 0xFFFFFFFF;

  /* so, if we have just the self-symmetry in the result, everything needs to be tried
   * and not rotations can be removed
   */
  if (!unSymmetric(resultSym)) {

    /* now we try to find the most "suitable" piece for our rotation removal. What is
     * suitable? Suitable is the piece shape that has the least common symmetries with
     * the result and that has the fewest pieces
     *
     * FIXME: if there is more than one suitable piece, select the one with the most
     * placements, this will gain us a little (or even bigger) speed-up
     * as its a difference if we select a piece that has only one placement anyway
     * or select one with 400 placements of which 23/24th can be dropped
     */
    unsigned int bestFound = sym->countSymmetryIntersection(resultSym, problem.getPartShape(0)->selfSymmetries());
    symBreakerShape = 0;

    for (unsigned int i = 1; i < problem.getNumberOfParts(); i++) {

      unsigned int cnt = sym->countSymmetryIntersection(resultSym, problem.getPartShape(i)->selfSymmetries());

      if (cnt < bestFound) {
        bestFound = cnt;
        symBreakerShape = i;
      }
    }

    if (sym->symmetriesLeft(resultSym, problem.getPartShape(symBreakerShape)->selfSymmetries()))
      checkForTransformedAssemblies(symBreakerShape, 0);

    if (sym->symmetryContainsMirror(resultSym)) {
      /* we need to to the mirror check here, and initialise the mirror
       * structure, otherwise no mirror check will be done
       */

      /* so, we need to find out which case this puzzle is, depending on the pieces
       * 1) all pieces contain mirror symmetries -> check mirrors, but no pairs
       * 2) as least one piece has no mirror symmetries
       *   2a) all pieces with no mirror symmetries have a mirror partner -> check mirrors, find pairs
       *   2b) at least one piece with no mirror symmetries has no partner -> no mirror check
       */

      typedef struct {
        unsigned int shape;    // the shape of this piece
        unsigned int mirror;   // the mirror shape of this piece
        unsigned int trans;
      } mm;

      mm * mirror = new mm[piecenumber];

      // first initialize
      for (unsigned int i = 0; i < problem.getNumberOfParts(); i++) {
        mirror[i].shape = i;
        mirror[i].mirror = (unsigned int)-1;
        mirror[i].trans = 255;
      }

      bool mirrorCheck = true;

      // now go over all shapes
      for (unsigned int i = 0; i < piecenumber; i++) {

        // we have already found the mirror for this shape
        if (mirror[i].mirror < piecenumber)
          continue;

        if (!sym->symmetryContainsMirror(problem.getPartShape(mirror[i].shape)->selfSymmetries())) {
          /* this shape is not self mirroring, so we need to look out
           * for a shape that is the mirror of this shape
           */
          bool found = false;

          // now see if we can find another shape that is the mirror of the current shape
          for (unsigned int j = i+1; j < piecenumber; j++) {

            if (mirror[j].mirror < piecenumber)
              continue;

            unsigned int trans = problem.getPartShape(mirror[i].shape)->getMirrorTransform(
                problem.getPartShape(mirror[j].shape));

            if (trans > 0) {
              // found a mirror shape

              mirror[i].mirror = j;
              mirror[i].trans = trans;
              mirror[j].mirror = i;
              mirror[j].trans = problem.getPartShape(mirror[j].shape)->getMirrorTransform(
                  problem.getPartShape(mirror[i].shape));

              found = true;
              break;
            }
          }

          // when we could not find a mirror transformation for the non mirrorable piece
          // we can stop and we don't need to make mirror checks
          if (!found) {
            mirrorCheck = false;
            break;
          }
        }
      }

      if (mirrorCheck) {
        /* all the shapes are either self mirroring or have a mirror pair
         * so we create the mirror structure and we do the mirror check
         */
        mirrorInfo_c * mir = new mirrorInfo_c();

        for (unsigned int i = 0; i < piecenumber; i++)
          if (mirror[i].trans != 255)
            mir->addPieces(i, mirror[i].mirror, mirror[i].trans);

        checkForTransformedAssemblies(symBreakerShape, mir);
      }

      delete [] mirror;
    }
  }

  voxel_c ** cache = new voxel_c *[sym->getNumTransformationsMirror()];

  literalPositions.clear();
  literalPositions.assign(1, assembler_2_c::piecePosition(0,0,0,0,0));
  lit = 0;

  std::vector<unsigned int> piecePlacements [piecenumber];
  std::vector<unsigned int> voxelPlacements [result->getXYZ()];

  /* now we insert one shape after another */
  for (unsigned int pc = 0; pc < problem.getNumberOfParts(); pc++) {

    reducePiece = pc;

    /* this array contains all the pieces found so far, this will help us
     * to not add two times the same piece to the structure */
    unsigned int cachefill = 0;
    unsigned int placements = 0;

    /* go through all possible rotations of the piece
     * if shape is new to cache, add it to the cache and also
     * add the shape to the matrix, in all positions that it fits
     */
    for (unsigned int rot = 0; rot < sym->getNumTransformations(); rot++) {

      voxel_c * rotation = gt->getVoxel(problem.getPartShape(pc));
      if (!rotation->transform(rot)) {
        delete rotation;
        continue;
      }

      rotation = addToCache(cache, &cachefill, rotation);

      if (rotation) {
        for (int x = (int)result->boundX1()-(int)rotation->boundX1(); x <= (int)result->boundX2()-(int)rotation->boundX2(); x++)
          for (int y = (int)result->boundY1()-(int)rotation->boundY1(); y <= (int)result->boundY2()-(int)rotation->boundY2(); y++)
            for (int z = (int)result->boundZ1()-(int)rotation->boundZ1(); z <= (int)result->boundZ2()-(int)rotation->boundZ2(); z++)
              if (canPlace(rotation, x, y, z)) {
                literalPositions.push_back(assembler_2_c::piecePosition(x, y, z, rot, pc));
                lit++;
                piecePlacements[pc].push_back(lit);
                placements = 1;

                /* now add the used cubes of the piece */
                for (unsigned int pz = rotation->boundZ1(); pz <= rotation->boundZ2(); pz++)
                  for (unsigned int py = rotation->boundY1(); py <= rotation->boundY2(); py++)
                    for (unsigned int px = rotation->boundX1(); px <= rotation->boundX2(); px++)
                      if (rotation->getState(px, py, pz) == voxel_c::VX_FILLED)
                        voxelPlacements[result->getIndex(x+px, y+py, z+pz)].push_back(lit);
              }


        /* for the symmetry breaker piece we also add all symmetries of the box */
        if (pc == symBreakerShape)
          for (unsigned int r = 1; r < sym->getNumTransformations(); r++)
            if (sym->symmetrieContainsTransformation(resultSym, r)) {

              voxel_c * vx = gt->getVoxel(problem.getPartShape(pc));

              if (!vx->transform(rot) || !vx->transform(r)) {
                delete vx;
                continue;
              }

              addToCache(cache, &cachefill, vx);
            }
      }
    }

    for (unsigned int i = 0; i < cachefill; i++)  delete cache[i];

    /* check, if the current piece has at least one placement */
    if (placements == 0) {
      delete [] cache;
      return -problem.getShapeIdOfPart(pc);
    }
  }

  /* ensure each piece is placed and not placed more than once */
  for (unsigned int pc = 0; pc < problem.getNumberOfParts(); pc++) {

    /* at least one placement*/
    for (unsigned int i = 0; i < piecePlacements[pc].size(); i++)
      solver->add(piecePlacements[pc][i]);
    solver->add(0);

    /* no more than one placement */
    for (unsigned int i = 0; i < piecePlacements[pc].size(); i++)
      for (unsigned int j = i+1; j < piecePlacements[pc].size(); j++) {
        solver->add(-piecePlacements[pc][i]);
        solver->add(-piecePlacements[pc][j]);
        solver->add(0);
      }
  }

  /* ensure each voxel in the result is filled if required and not double filled*/
  for (unsigned int vox = 0; vox < result->getXYZ(); vox++) {

    /* at least one placement*/
    if (result->getState(vox) == voxel_c::VX_FILLED) {
      for (unsigned int i = 0; i < voxelPlacements[vox].size(); i++)
        solver->add(voxelPlacements[vox][i]);
      solver->add(0);
    }

    /* no more than one placement */
    for (unsigned int i = 0; i < voxelPlacements[vox].size(); i++)
      for (unsigned int j = i+1; j < voxelPlacements[vox].size(); j++) {
        solver->add(-voxelPlacements[vox][i]);
        solver->add(-voxelPlacements[vox][j]);
        solver->add(0);
      }
  }

  delete [] cache;

  return 1;
}


assembler_2_c::errState assembler_2_c::createMatrix(bool keepMirror, bool keepRotations, bool comp) {

  bt_assert(problem.resultValid());

  complete = comp;

  if (!canHandle(problem))
    return ERR_PUZZLE_UNHANDABLE;

  /* get and save piece number of puzzle */
  piecenumber = problem.getNumberOfPieces();

  /* count the filled and variable units */
  int res_vari = getResultShape(problem)->countState(voxel_c::VX_VARIABLE);
  int res_filled = getResultShape(problem)->countState(voxel_c::VX_FILLED) + res_vari;

  // check if number of voxels in pieces is not bigger than
  // number of voxel in result

  // check if number of filled voxels in result
  // is not bigger than number of voxels in pieces
  int h = res_filled;

  for (unsigned int j = 0; j < problem.getNumberOfParts(); j++)
    h -= problem.getPartShape(j)->countState(voxel_c::VX_FILLED);

  if (h < 0) {
    errorsState = ERR_TOO_MANY_UNITS;
    errorsParam = -h;
    return errorsState;
  }

  if (h > res_vari) {
    errorsState = ERR_TOO_FEW_UNITS;
    errorsParam = h-res_vari;
    return errorsState;
  }

  holes = h;

  /* fill the nodes arrays */
  int error = prepare();

  // check, if there is one piece not placeable
  if (error <= 0) {
    errorsState = ERR_CAN_NOT_PLACE;
    errorsParam = -error;
    return errorsState;
  }

  iterations = 0;

  if (keepMirror)
    avoidTransformedMirror = 0;

  if (keepRotations)
    avoidTransformedAssemblies = false;

  errorsState = ERR_NONE;
  return errorsState;
}

void assembler_2_c::reduce(void) {

  // TODO: cadical reduce
}

assembly_c * assembler_2_c::getAssembly(void) {

  assembly_c * assembly = new assembly_c(problem.getPuzzle().getGridType());

  // if no pieces are placed, or we finished return an empty assembly
  if (solverState != CaDiCaL::SATISFIED) {
    for (unsigned int i = 0; i < piecenumber; i++)
      assembly->addNonPlacement();
    return assembly;
  }

    /* first we need to find the order the piece are in */
  int * pieces = new int[piecenumber];
  unsigned char * trans = new unsigned char[piecenumber];
  int * xs = new int[piecenumber];
  int * ys = new int[piecenumber];
  int * zs = new int[piecenumber];

  /* fill the array with 0xff, so that we can distinguish between
   * placed and unplaced pieces
   */
  memset(pieces, -1, sizeof(int) * piecenumber);

  int usedPoses [piecenumber];

  for (unsigned int lit = 1; lit < literalPositions.size(); lit++) {
    if (solver->val(lit) > 0) {
      assembler_2_c::piecePosition pos = literalPositions[lit];
      int piece = pos.piece;

      usedPoses[piece] = lit;

      pieces[piece] = lit;
      trans[piece] = pos.transformation;
      xs[piece] = pos.x;
      ys[piece] = pos.y;
      zs[piece] = pos.z;
    }
  }

  for (unsigned int pc = 0; pc < piecenumber; pc++)
    solver->add(-usedPoses[pc]);
  solver->add(0);

  for (unsigned int i = 0; i < piecenumber; i++)
    if (pieces[i] < 0)
      assembly->addNonPlacement();
    else {
      fprintf(stderr, "piece: %d, trans: %d, xyz: (%d, %d, %d)\n", i, trans[i], xs[i], ys[i], zs[i]);
      assembly->addPlacement(trans[i], xs[i], ys[i], zs[i]);
    }

  delete [] pieces;
  delete [] trans;
  delete [] xs;
  delete [] ys;
  delete [] zs;

  // sort is not necessary because there is only one of each piece
  // assembly->sort(puzzle, problem);

  return assembly;
}

void assembler_2_c::checkForTransformedAssemblies(unsigned int pivot, mirrorInfo_c * mir) {
  avoidTransformedAssemblies = true;
  avoidTransformedPivot = pivot;
  avoidTransformedMirror = mir;
}

/* this function handles the assemblies found by the assembler engine
 */
void assembler_2_c::solution(void) {

  if (getCallback()) {

    assembly_c * assembly = getAssembly();

    if (avoidTransformedAssemblies && assembly->smallerRotationExists(problem, avoidTransformedPivot, avoidTransformedMirror, complete))
      delete assembly;
    else {
      getCallback()->assembly(assembly);
    }
  }
}

void assembler_2_c::iterativeSAT() {
  terminator = new terminator_c();
  solver->connect_terminator(terminator);

  int result;
  while (result = solver->solve(), result == 10) {
    solverState = CaDiCaL::SATISFIED;
    solution();
  }
  if (result) solverState = CaDiCaL::UNSATISFIED;
  else solverState = CaDiCaL::UNKNOWN;

  solver->disconnect_terminator();
  delete terminator;
}

void assembler_2_c::assemble(assembler_cb * callback) {

  debug = false;

  if (errorsState == ERR_NONE) {
    asm_bc = callback;
    iterativeSAT();
  }
}

float assembler_2_c::getFinished(void) const {
  return 0;
  //TODO: cadical finished
}

void assembler_2_c::stop(void) {
  if (terminator) terminator->abort = true;
}

assembler_c::errState assembler_2_c::setPosition([[maybe_unused]] const char * string, const char * version) {

  /* check for the right version */
  if (strcmp(version, ASSEMBLER_VERSION))
    return ERR_CAN_NOT_RESTORE_VERSION;

  // TODO: restore from cadical save

  return ERR_CAN_NOT_RESTORE_SYNTAX;
}

void assembler_2_c::save(xmlWriter_c & xml) const
{
  xml.newTag("assembler");
  xml.newAttrib("version", ASSEMBLER_VERSION);

  // TODO: save solver and literal poses

  xml.endTag("assembler");
}

void assembler_2_c::debug_step(unsigned long num) {
  debug = true;
  debug_loops = num;
  asm_bc = 0;
  // TODO: cadical limit and run
}

bool assembler_2_c::canHandle(const problem_c & p) {

  // we can not handle if there is one shape having not a counter of 1
  for (unsigned int s = 0; s < p.getNumberOfParts(); s++)
    if ((p.getPartMaximum(s) > 1) ||
        (p.getPartMaximum(s) != p.getPartMinimum(s)))

      return false;

  return true;
}

