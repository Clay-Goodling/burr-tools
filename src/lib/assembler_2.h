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
#ifndef __ASSEMBLER_2_H__
#define __ASSEMBLER_2_H__

#include "assembler.h"

#include <vector>
#include <set>
#include <stack>
#include <memory>

#include "cadical/src/cadical.hpp"

class gridType_c;
class mirrorInfo_c;

/**
 * This is an assembler class.
 *
 * It is more or less identical to Don Knuths idea. Some changes have been done though
 * to provide for holes. This class can not handle ranges or multi-pieces.
 *
 * All involved pieces must be there exactly one time. But in that case it is a bit
 * faster than assembler_1.
 */
class assembler_2_c : public assembler_c {

protected:

  const problem_c & problem;

private:

  /* the sat solver */
  CaDiCaL::Solver * solver;
  CaDiCaL::State solverState;

  class terminator_c : public CaDiCaL::Terminator {
  public:
    bool abort;

    terminator_c(): abort(false) {};

    bool terminate () {
      return abort;
    }
  };
  terminator_c * terminator;

  /* used to abort the searching */
  bool abbort;

  /* used to save if the search is running */
  bool running;

  /* this function gets called whenever an assembly was found
   * when a call-back is available it will call getAssembly to
   * obtain the assembly for the found solution when the
   * field avoidTransformedAssemblies is true then the assembly
   * is checked, if it has been found before. The assembly
   * is normalised in inserted into a set of assemblies for
   * later reference
   */
  void solution(void);

  void iterativeSAT(void);

  /* this function checks, if the given piece can be placed
   * at the given position inside the result
   */
  bool canPlace(const voxel_c * piece, int x, int y, int z) const;

  /* this function creates the matrix for the search function
   * because we need to know how many nodes we need to allocate the
   * arrays with the right size, we add a parameter. If this is true
   * the function will not access the array but only count the number
   * of nodes used. This number is returned
   *
   * return error codes
   */
  int prepare(void);

  /* internal error state */
  errState errorsState;
  int errorsParam;

  /* number of iterations the assemble routine run */
  unsigned long iterations;

  /* the number of holes the assembles piece will have. Holes are
   * voxels in the variable voxel set that are not filled. The other
   * voxels are all filled
   */
  int holes;

  /* now this isn't hard to guess, is it? */
  unsigned int piecenumber;

  /* the message object that gets called with the solutions as param */
  assembler_cb * asm_bc;

  /* this value contains the piecenumber that the reduce procedure is currently working on
   * the value is only valid, when reduce is running
   */
  unsigned int reducePiece;
  
  /* this vector contains the placement (transformation and position) for
   * a piece in a row
   */
  class piecePosition {

  public:

    int x, y, z;
    unsigned char transformation;
    unsigned int piece;

    piecePosition(int x_, int y_, int z_, unsigned char transformation_, unsigned int pc): 
    x(x_), y(y_), z(z_), transformation(transformation_), piece(pc) {}
  };
  std::vector<piecePosition> literalPositions;

  /* the number of used literals*/
  int lit;

  /* the members for rotations rejection
   */
  bool avoidTransformedAssemblies;
  unsigned int avoidTransformedPivot;
  mirrorInfo_c * avoidTransformedMirror;

  /// set to true, when complete rotation analysis is requested
  bool complete;

  /* the variables for debugging assembling processes
   */
  bool debug;         // debugging enabled
  int debug_loops;    // how many loops to run ?

protected:

  /* finally after assembling a puzzle and creating something meaningful from the cover
   * information you need to call the callback of the user, use this function to get the
   * callback class
   */
  assembler_cb * getCallback(void) { return asm_bc; }

  unsigned int getPiecenumber(void) { return piecenumber; }

  /* call this function if you think that there might be
   * rotated assemblies found. Here a description of how the whole aspect of
   * rotation avoiding is supposed to work
   * the front end is supposed to initialise the assembler so that as few as
   * possible double assemblies are found by selecting one piece and not placing
   * this piece in all possible positions. But this will not always work, if
   * the front end is are not absolutely certain that it has avoided all possible
   * rotations it should call this function. This will then add an additional check
   * for each found assembly
   */
  void checkForTransformedAssemblies(unsigned int pivot, mirrorInfo_c * mir);

public:

  assembler_2_c(const problem_c & problem);
  ~assembler_2_c(void);

  /* functions that are overloaded from assembler_c, for comments see there */
  errState createMatrix(bool keepMirror, bool keepRotations, bool complete);
  void assemble(assembler_cb * callback);
  int getErrorsParam(void) { return errorsParam; }
  virtual float getFinished(void) const;
  virtual void stop(void);
  virtual bool stopped(void) const { return !terminator; }
  virtual errState setPosition(const char * string, const char * version);
  virtual void save(xmlWriter_c & xml) const;
  virtual void reduce(void);
  virtual unsigned int getReducePiece(void) const { return reducePiece; }
  virtual unsigned long getIterations(void) { return iterations; }

  /* some more special information to find out possible piece placements */
  bool getPiecePlacementSupported(void) const { return false; }
  unsigned int getPiecePlacement(
    [[maybe_unused]] unsigned int node,
    [[maybe_unused]] int delta,
    [[maybe_unused]] unsigned int piece,
    [[maybe_unused]] unsigned char *tran,
    [[maybe_unused]] int *x,
    [[maybe_unused]] int *y,
    [[maybe_unused]] int *z) const { return 0;};
  unsigned int getPiecePlacementCount([[maybe_unused]] unsigned int piece) const { return 0; };

  void debug_step(unsigned long num = 1);
  assembly_c * getAssembly(void);

  static bool canHandle(const problem_c & p);

private:

  // no copying and assigning
  assembler_2_c(const assembler_2_c&);
  void operator=(const assembler_2_c&);
};

#endif
