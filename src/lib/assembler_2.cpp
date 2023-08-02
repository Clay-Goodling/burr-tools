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

int assembler_2_c::atMostOne(std::vector<unsigned int> lits, int i0, int i1) {

  bt_assert(i0 <= i1);

  int k = i1 - i0;
  if (k <= 4) {
    for (int i = i0; i < i1; i++)
      for (int j = i+1; j < i1; j++) {
        solver->add(-lits[i]);
        solver->add(-lits[j]);
        solver->add(0);
      }
    return (k * (k - 1)) >> 1;
  }

  int m = ++lit;
  for (int i = i0; i < i0 + (k >> 1); i++) {
    solver->add(-lits[i]);
    solver->add(-m);
    solver->add(0);
  }

  for (int i = i0 + (k >> 1); i < i1; i++) {
    solver->add(-lits[i]);
    solver->add(m);
    solver->add(0);
  }

  return k + atMostOne(lits, i0, i0 + (k >> 1)) + atMostOne(lits, i0 + (k >> 1), i1);
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
    unsigned int symBreakerPiece = 0;
    unsigned int pc = problem.getPartMaximum(0);
    unsigned int bestFound = sym->countSymmetryIntersection(resultSym, problem.getPartShape(0)->selfSymmetries());
    symBreakerShape = 0;

    for (unsigned int i = 1; i < problem.getNumberOfParts(); i++) {

      unsigned int cnt = sym->countSymmetryIntersection(resultSym, problem.getPartShape(i)->selfSymmetries());

      if ((problem.getPartMaximum(i) < problem.getPartMaximum(symBreakerShape)) ||
          ((problem.getPartMaximum(i) == problem.getPartMaximum(symBreakerShape)) && (cnt < bestFound))) {
        bestFound = cnt;
        symBreakerShape = i;
        symBreakerPiece = pc;
      }

      pc += problem.getPartMaximum(i);
    }

    bool tmp = sym->symmetriesLeft(resultSym, problem.getPartShape(symBreakerShape)->selfSymmetries());

    bool pieceRanges = false;
    for (unsigned int i = 0; i < problem.getNumberOfParts(); i++)
      if (problem.getPartMinimum(i) != problem.getPartMaximum(i)) {
        pieceRanges = true;
        break;
      }

    if (tmp || (problem.getPartMaximum(symBreakerShape) > 1) || pieceRanges) {

      // we can not use the symmetry breaker shape, if there is more than one piece
      // of this shape in the problem
      if (pieceRanges || problem.getPartMaximum(symBreakerShape) > 1) {
        symBreakerShape = 0xFFFFFFFF;
        symBreakerPiece = 0xFFFFFFFF;
      }

      checkForTransformedAssemblies(symBreakerPiece, 0);
    }


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

      pc = 0;
      // first initialize
      for (unsigned int i = 0; i < problem.getNumberOfParts(); i++)
        for (unsigned int p = 0; p < problem.getPartMaximum(i); p++) {
          mirror[pc].shape = i;
          mirror[pc].mirror = (unsigned int)-1;
          mirror[pc].trans = 255;
          pc++;
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

      if (mirrorCheck || pieceRanges) {
        /* all the shapes are either self mirroring or have a mirror pair
         * so we create the mirror structure and we do the mirror check
         * we also need to that when ranges are used because the final solution
         * might use only mirrorable pieces and then we need this information
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

  std::vector<unsigned int> piecePlacements [problem.getNumberOfParts()];
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
                literalPositions.push_back(assembler_2_c::piecePosition(x+rotation->getHx(), y+rotation->getHy(), z+rotation->getHz(), rot, pc));
                lit++;
                piecePlacements[pc].push_back(lit);
                placements++;

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

    /* check, if the current part has at least the minimum required placements */
    if (placements < problem.getPartMinimum(pc)) {
      delete [] cache;
      return -problem.getShapeIdOfPart(pc);
    }
  }

  int ppc = 0;
  /* ensure each part is placed an acceptable number of times */
  for (unsigned int pc = 0; pc < problem.getNumberOfParts(); pc++) {
    fprintf(stderr, "piece %d, placements %lu, min %d, max %d\n", pc, piecePlacements[pc].size(), problem.getPartMinimum(pc), problem.getPartMaximum(pc));
    ppc += cardinalityConstraint(piecePlacements[pc], problem.getPartMinimum(pc), problem.getPartMaximum(pc));
  }
  fprintf(stderr, "piece placement clauses: %d\n", ppc);

  int vuc = 0;
  /* ensure each voxel in the result is filled if required and not double filled*/
  for (unsigned int vox = 0; vox < result->getXYZ(); vox++) {
    vuc += cardinalityConstraint(voxelPlacements[vox], result->getState(vox) == voxel_c::VX_FILLED, 1);
  }
  fprintf(stderr, "voxel usage clauses: %d\n", vuc);

  delete [] cache;

  return 1;
}

int assembler_2_c::cardinalityConstraint(std::vector<unsigned int> lits, unsigned int min, unsigned int max) {
  bt_assert(min <= max);

  if (lits.size() == 0) return 0;
  if (max == 0) {
    for (unsigned int i = 0; i < lits.size(); i++) {
      solver->add(-lits[i]), solver->add(0);
    }
    return lits.size();
  }

  int counters = ++lit;
  #define COUNTER(row, column) (counters + (row)*max + (column))
  lit = COUNTER(lits.size()-1, max-1);

  // first lit iff first bit of first counter
  solver->add(-lits[0]), solver->add(COUNTER(0,0)), solver->add(0);
  solver->add(lits[0]), solver->add(-COUNTER(0,0)), solver->add(0);

  // rest of first counter is zero
  for (unsigned int i = 1; i < max; i++) {
    solver->add(-COUNTER(0, i)), solver->add(0);
  }
  int clauses = 1 + max;
  
  // next n-1 counters count correctly
  for (unsigned int i = 1; i < lits.size(); i++) {
    // first bit iff lit or first bit of last counter
    solver->add(-lits[i]), solver->add(COUNTER(i, 0)), solver->add(0);
    solver->add(-COUNTER(i-1, 0)), solver->add(COUNTER(i, 0)), solver->add(0);
    solver->add(-COUNTER(i, 0)), solver->add(lits[i]), solver->add(COUNTER(i-1, 0)), solver->add(0);
    clauses += 3;

    // next bits iff lit and last bit of last counter or same bit of last counter
    for (unsigned int j = 1; j < max; j++) {
      solver->add(-lits[i]), solver->add(-COUNTER(i-1, j-1)), solver->add(COUNTER(i, j)), solver->add(0);
      solver->add(-COUNTER(i-1, j)), solver->add(COUNTER(i, j)), solver->add(0);
      solver->add(-COUNTER(i, j)), solver->add(COUNTER(i-1, j)), solver->add(lits[i]), solver->add(0);
      solver->add(-COUNTER(i, j)), solver->add(COUNTER(i-1, j)), solver->add(COUNTER(i-1, j-1)), solver->add(0);
      clauses += 4;
    }

    // no overflow
    solver->add(-lits[i]), solver->add(-COUNTER(i-1, max-1)), solver->add(0); 
    clauses++;
  }

  // if min > 0, require min bit of last counter to be 1
  if (min > 0) {
    solver->add(COUNTER(lits.size()-1, min-1)), solver->add(0);
    clauses++;
  }

  return clauses;
}

assembler_2_c::errState assembler_2_c::createMatrix(bool keepMirror, bool keepRotations, bool comp) {

  bt_assert(problem.resultValid());

  complete = comp;

  if (!canHandle(problem))
    return ERR_PUZZLE_UNHANDABLE;

  /* get and save piece number of puzzle */
  piecenumber = problem.getNumberOfPieces();

  /* count the filled and variable units */
  unsigned int res_vari = getResultShape(problem)->countState(voxel_c::VX_VARIABLE);
  unsigned int res_filled = getResultShape(problem)->countState(voxel_c::VX_FILLED) + res_vari;

  // check if number of voxels in pieces is not bigger than
  // number of voxel in result

  // check if number of filled voxels in result
  // is not bigger than number of voxels in pieces
  unsigned int min = 0;
  unsigned int max = 0;

  for (unsigned int j = 0; j < problem.getNumberOfParts(); j++) {
    min += problem.getPartShape(j)->countState(voxel_c::VX_FILLED) * problem.getPartMinimum(j);
    max += problem.getPartShape(j)->countState(voxel_c::VX_FILLED) * problem.getPartMaximum(j);
  }

  if (min == max)
    holes = res_filled - min;
  else if (problem.maxHolesDefined())
    holes = problem.getMaxHoles();
  else
    holes = 0xFFFFFF;

  if (min > res_filled) {
    errorsState = ERR_TOO_MANY_UNITS;
    errorsParam = min-res_filled;
    return errorsState;
  }

  if (max < res_filled-res_vari) {
    errorsState = ERR_TOO_FEW_UNITS;
    errorsParam = res_filled-res_vari-max;
    return errorsState;
  }

  /* build the sat reduction */
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
  solver->optimize(2);
}

assembly_c * assembler_2_c::getAssembly(void) {

  assembly_c * assembly = new assembly_c(problem.getPuzzle().getGridType());

  // if no pieces are placed, or we finished return an empty assembly
  if (solverState != CaDiCaL::SATISFIED) {
    for (unsigned int i = 0; i < piecenumber; i++)
      assembly->addNonPlacement();
    return assembly;
  }

  std::vector<int> usedLits;
  std::vector<assembler_2_c::piecePosition> placements [problem.getNumberOfParts()];

  for (unsigned int placement = 1; placement < literalPositions.size(); placement++) {
    if (solver->val(placement) > 0) {
      assembler_2_c::piecePosition pos = literalPositions[placement];
      int piece = pos.piece;
      placements[piece].push_back(pos);
      usedLits.push_back(placement);
    }
  }

  for (unsigned int pc = 0; pc < usedLits.size(); pc++)
    solver->add(-usedLits[pc]);
  solver->add(0);

  for (unsigned int i = 0; i < problem.getNumberOfParts(); i++) {
    unsigned int j = 0;
    while (j < placements[i].size()) {
      assembly->addPlacement(
        placements[i][j].transformation,
        placements[i][j].x,
        placements[i][j].y,
        placements[i][j].z);
      j++;
    }

    while (j < problem.getPartMaximum(i)) {
      assembly->addNonPlacement();
      j++;
  }
  }

  assembly->sort(problem);
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
  return solverState && solverState == CaDiCaL::UNSATISFIED;
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

bool assembler_2_c::canHandle(const problem_c &) {
  return true;
}

