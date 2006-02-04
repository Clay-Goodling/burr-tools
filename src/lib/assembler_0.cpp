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


#include "assembler_0.h"

#include <xmlwrapp/attributes.h>

#ifdef WIN32
#define snprintf _snprintf
#endif

void assembler_0_c::GenerateFirstRow(unsigned int res_filled) {

  for (unsigned int i = 0; i < varivoxelStart; i++) {
    right.push_back(i+1);
    left.push_back(i-1);
    up.push_back(i);
    down.push_back(i);
    colCount.push_back(0);
  }

  /* make the linked list cyclic */
  left[0] = varivoxelStart - 1;
  right[varivoxelStart - 1] = 0;

  /* create column nodes for vari columns, these
   * are not inside the column header list and so
   * are not forced to be filled but otherwise behave
   * like normal columns
   */
  for (unsigned int j = varivoxelStart; j <= varivoxelEnd; j++) {
    left.push_back(j-1);
    right.push_back(j+1);
    up.push_back(j);
    down.push_back(j);
    colCount.push_back(0);
  }
  left[varivoxelStart] = varivoxelEnd;
  right[varivoxelEnd] = varivoxelStart;
}

int assembler_0_c::AddPieceNode(unsigned int piece, unsigned int rot, unsigned int x, unsigned int y, unsigned int z) {
  unsigned long piecenode = up.size();

  left.push_back(piecenode);
  right.push_back(piecenode);
  down.push_back(piece+1);
  up.push_back(up[piece+1]);

  down[up[piece+1]] = piecenode;
  up[piece+1] = piecenode;

  colCount.push_back(piece+1);
  colCount[piece+1]++;

  piecePositions.push_back(piecePosition(x, y, z, rot, piecenode));

  return piecenode;
}

void assembler_0_c::getPieceInformation(unsigned int node, unsigned char *tran, int *x, int *y, int *z) {

  for (int i = piecePositions.size()-1; i >= 0; i--)
    if (piecePositions[i].row <= node) {
      *tran = piecePositions[i].transformation;
      *x = piecePositions[i].x;
      *y = piecePositions[i].y;
      *z = piecePositions[i].z;

      return;
    }

  bt_assert(0);
}

void assembler_0_c::AddVoxelNode(unsigned int col, unsigned int piecenode) {
  unsigned long newnode = up.size();

  right.push_back(piecenode);
  left.push_back(left[piecenode]);
  right[left[piecenode]] = newnode;
  left[piecenode] = newnode;

  down.push_back(col);
  up.push_back(up[col]);
  down[up[col]] = newnode;
  up[col] = newnode;

  colCount.push_back(col);
  colCount[col]++;
}

void assembler_0_c::nextPiece(unsigned int piece, unsigned int count, unsigned int number) {
  multiPieceCount[piece] = count;
  multiPieceIndex[piece] = number;
  pieceStart[piece] = up.size();
}

assembler_0_c::assembler_0_c() :
  left(0), right(0), up(0), down(0), colCount(0),
  multiPieceCount(0), multiPieceIndex(0), pieceStart(0),
  pos(0), rows(0), columns(0), nodeF(0), numF(0),
  piece(0), searchState(0), addRows(0), avoidTransformedAssemblies(0)
{
}

assembler_0_c::~assembler_0_c() {
  if (rows) delete [] rows;
  if (columns) delete [] columns;
  if (multiPieceCount) delete [] multiPieceCount;
  if (pieceStart) delete [] pieceStart;
  if (multiPieceIndex) delete [] multiPieceIndex;
  if (nodeF) delete [] nodeF;
  if (numF) delete [] numF;
  if (piece) delete [] piece;
  if (searchState) delete [] searchState;
  if (addRows) delete [] addRows;
}

assembler_0_c::errState assembler_0_c::createMatrix(const puzzle_c * puz, unsigned int prob) {

  puzzle = puz;
  problem = prob;

  /* get and save piecenumber of puzzle */
  piecenumber = puz->probPieceNumber(prob);

  /* count the filled and variable units */
  int res_vari = puz->probGetResultShape(prob)->countState(voxel_c::VX_VARIABLE);
  int res_filled = puz->probGetResultShape(prob)->countState(voxel_c::VX_FILLED) + res_vari;

  for (unsigned int i = 0; i < puz->probShapeNumber(prob); i++)
    if (puz->probGetShapeShape(prob, i)->countState(voxel_c::VX_VARIABLE)) {
      errorsParam = puz->probGetShape(prob, i);
      errorsState = ERR_PIECE_WITH_VARICUBE;
      return errorsState;
    }

  varivoxelStart = 1 + piecenumber + res_filled - res_vari;
  varivoxelEnd = 1 + piecenumber + res_filled;

  // check if number of voxels in pieces is not bigger than
  // number of voxel in result

  // check if number of filled voxels in result
  // is not bigger than number of voxels in pieces
  int h = res_filled;

  for (unsigned int j = 0; j < puz->probShapeNumber(prob); j++)
    h -= puz->probGetShapeShape(prob, j)->countState(voxel_c::VX_FILLED) * puz->probGetShapeCount(prob, j);

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

  /* allocate all the required memory */
  rows = new unsigned int[piecenumber];
  columns = new unsigned int [piecenumber];
  nodeF = new unsigned int [piecenumber];
  numF = new unsigned int [piecenumber];
  piece = new unsigned int [piecenumber];
  addRows = new std::stack<unsigned int> [piecenumber];
  searchState = new unsigned int [piecenumber + 1];

  searchState[0] = 0;

  pieceStart = new unsigned int[piecenumber];
  multiPieceCount = new unsigned int[piecenumber];
  multiPieceIndex = new unsigned int[piecenumber];

  /* fill the nodes arrays */
  int error = prepare(puz, res_filled, res_vari, prob);

  // check, if there is one piece unplacable
  if (error <= 0) {
    errorsState = ERR_CAN_NOT_PLACE;
    errorsParam = -error;
    return errorsState;
  }

  memset(rows, 0, piecenumber * sizeof(int));
  memset(columns, 0, piecenumber * sizeof(int));
  pos = 0;
  iterations = 0;

  errorsState = ERR_NONE;
  return errorsState;
}


/* remove column from array, and also all the rows, where the column is one */
void assembler_0_c::cover(register unsigned int col)
{
  {
    register unsigned int l, r;

    l = left[col];
    r = right[col];

    left[r] = l;
    right[l] = r;
  }

  register unsigned int i = down[col];
  while (i != col) {

    register unsigned int j = right[i];
    while (j != i) {

      {
        register int u, d;

        u = up[j];
        d = down[j];

        up[d] = u;
        down[u] = d;
      }

      colCount[colCount[j]]--;

      j = right[j];
    }

    i = down[i];
  }
}

/* remove column from array, and also all the rows, where the column is one */
bool assembler_0_c::try_cover(register unsigned int col, unsigned int * columns)
{
  {
    register unsigned int l, r;

    l = left[col];
    r = right[col];

    left[r] = l;
    right[l] = r;
  }

  bool result = true;

  register unsigned int i = down[col];
  while (i != col) {

    register unsigned int j = right[i];
    while (j != i) {

      {
        register unsigned int u, d;

        u = up[j];
        d = down[j];

        up[d] = u;
        down[u] = d;
      }

      colCount[colCount[j]]--;

      if ((colCount[colCount[j]] == 0) && (columns[colCount[j]] == 0))
        result = false;

      j = right[j];
    }

    i = down[i];
  }

  return result;
}

void assembler_0_c::uncover(register unsigned int col) {

  register unsigned int i = up[col];
  while (i != col) {

    register unsigned int j = left[i];
    while (j != i) {

      colCount[colCount[j]]++;

      up[down[j]] = j;
      down[up[j]] = j;

      j = left[j];
    }

    i = up[i];
  }

  left[right[col]] = col;
  right[left[col]] = col;
}

/* remove all the columns from the matrix in which the given
 * row contains ones
 */
void assembler_0_c::cover_row(register unsigned int r) {
  register unsigned int j = right[r];
  while (j != r) {
    cover(colCount[j]);
    j = right[j];
  }
}

bool assembler_0_c::try_cover_row(register unsigned int r) {

  unsigned int *columns = new unsigned int[varivoxelEnd];
  memset(columns, 0, varivoxelEnd * sizeof(int));

  register unsigned int j = right[r];
  while (j != r) {
    columns[colCount[j]] = 1;
    j = right[j];
  }

  j = right[varivoxelEnd];
  while (j != varivoxelEnd) {
    columns[j] = 1;
    j = right[j];
  }


  j = right[r];
  while (j != r) {

    if (!try_cover(colCount[j], columns)) {
      do {
        uncover(colCount[j]);
        j = left[j];
      } while (j != r);

      delete [] columns;
      return false;
    }

    j = right[j];
  }

  delete [] columns;
  return true;
}

void assembler_0_c::uncover_row(register unsigned int r) {
  register unsigned int j = left[r];
  while (j != r) {
    uncover(colCount[j]);
    j = left[j];
  }
}

void assembler_0_c::remove_row(register unsigned int r) {
  register unsigned int j = r;
  do {
    register unsigned int u, d;

    colCount[colCount[j]]--;

    u = up[j];
    d = down[j];

    up[d] = u;
    down[u] = d;
    j = right[j];
  } while (j != r);
}

void assembler_0_c::reinsert_row(register unsigned int r) {
  register unsigned int j = r;
  do {
    up[down[j]] = j;
    down[up[j]] = j;

    colCount[colCount[j]]++;

    j = right[j];
  } while (j != r);
}


bool assembler_0_c::checkmatrix(unsigned int rec, unsigned int branch) {

  /* check the number of holes, if they are larger than allowed return */
  unsigned int j = right[varivoxelEnd];
  unsigned int count = holes;
  while (j != varivoxelEnd) {
    if (colCount[j] == 0) {
      if (!count)
        return true;
      count--;
    }
    j = right[j];
  }

  unsigned int col = right[0];

  while (col) {
    /* if there is a column that can not be covered
     * any longer, the whole problem can't be solved
     */
    if (colCount[col] == 0)
      return true;

    /* if one column must be covered in less than branch
     * level, check all these possibilities
     */
    if (rec && (colCount[col] <= branch)) {

      cover(col);

      unsigned int r = down[col];

      while (r != col) {

        if (!try_cover_row(r)) {

          uncover(col);
          return true;

        }

        bool ret = checkmatrix(rec - 1, branch);
        uncover_row(r);

        if (ret) {
          uncover(col);
          return true;
        }

        r = down[r];
      }

      uncover(col);
    }

    col = right[col];
  }

  return false;
}

void assembler_0_c::reduce(void) {

  unsigned int removed = 0;
  bool rem_sth;
  unsigned int iteration = 1;

  unsigned int rec = 0;
  unsigned int branch = 0;

  do {

    rem_sth = false;

    /* check all the pieces */
    for (unsigned int p = 0; p < piecenumber;) {

      reducePiece = p;

      // place the piece and check, if this leads to some
      // infillable holes or unplacable pieces or whaever
      // conditions that make a solution impossible

      cover(p+1);

      // we first collect all the rows that we finally want to
      // remove and only remove them after the complete check
      int * rowsToRemove = new int[colCount[p+1]];
      int countToRemove = 0;

      // go over all the placements of the piece and check, if
      // each for possibility
      unsigned int r = down[p+1];
      unsigned int row = 0;
      while (r != p+1) {

        // try to do this placement, if the placing goes
        // wrong already, we don't need to do the deep check
        if (!try_cover_row(r)) {
          rowsToRemove[countToRemove++] = row;
        } else {

          /* and check if that results in a dead end */
          if (checkmatrix(rec, branch))
            rowsToRemove[countToRemove++] = row;

          uncover_row(r);
        }

        r = down[r];
        row++;
      }

      uncover(p+1);

//      int count = multiPieceCount[p];

      for (int pc = 0; pc < 1/*count*/; pc++) {

        int r = down[p+1];
        int row = 0;

        for (int rem = 0; rem < countToRemove; rem++) {

          while (row < rowsToRemove[rem]) {
            r = down[r];
            row++;
          }

          int s = r;
          int rr = down[r];

          do {
            up[down[s]] = up[s];
            down[up[s]] = down[s];
            colCount[colCount[s]]--;
            s = right[s];
          } while (s != r);

          r = rr;
          row++;

          removed++;
        }

        p++;
      }

      if (countToRemove)
        rem_sth = true;

      delete [] rowsToRemove;
    }

//    return;

//    iteration++;
//    rec++;
//    if (rec > 2)
//      rec = 2;
//    branch = 1;

  } while (rem_sth);
}

assembly_c * assembler_0_c::getAssembly(void) {

  assembly_c * assembly = new assembly_c();

  // if no pieces are placed, or we finished return an empty assembly
  if (getPos() <= 0) return assembly;

  /* first we need to find the order the piece are in */
  unsigned int * pieces = new unsigned int[getPiecenumber()];

  /* fill the array with 0xff, so that we can distinguish between
   * placed and unplaced pieces
   */
  memset(pieces, 0xff, sizeof(unsigned int) * getPiecenumber());

  for (unsigned int i = 0; i < getPos(); i++) {
    bt_assert(getPiece(i) < getPiecenumber());
    pieces[getPiece(i)] = i;
  }

  for (unsigned int i = 0; i < getPiecenumber(); i++)
    if (pieces[i] >= getPos())
      assembly->addNonPlacement();
    else {
      unsigned char tran;
      int x, y, z;

      getPieceInformation(getRows(pieces[i]), &tran, &x, &y, &z);
      assembly->addPlacement(tran, x, y, z);
    }

  delete [] pieces;

  assembly->sort(puzzle, problem);

  return assembly;
}


void assembler_0_c::checkForTransformedAssemblies(unsigned int pivot) {
  avoidTransformedAssemblies = true;
  avoidTransformedPivot = pivot;
}


/* this function handles the assemblies found by the assembler engine
 */
void assembler_0_c::solution(void) {

  if (getCallback()) {

    assembly_c * assembly = getAssembly();

    if (avoidTransformedAssemblies && assembly->smallerRotationExists(puzzle, problem, avoidTransformedPivot))
      delete assembly;
    else
      getCallback()->assembly(assembly);
  }
}

/* to understand this function you need to first completely understand the
 * dancing link algorithm.
 * The next thing to know is how this alg. avoids finding the solutions
 * with identical pieces with all permutations. This is done by enforcing
 * an order of the pieces. The placements and the pieces are enumerated and
 * pieces with a smaller number must have a placement with a smaller number.
 * Then you need to know the layout of the matrix used here. It's exactly as
 * describes by Knuth. First columns for pieces. Followed by the columns for
 * the normal result voxels followed by the columnns for the variable
 * pieces. All placements for one piece are one below the other. One piece
 * follows after another. The nodes are enumerated line by line.
 * This fixed order allows for simple checks of the placement numbers against
 * each other (only relative comparisons are possible and necessary) by
 * knowing the piece this row belongs to and knowing the first node that
 * belongs to this piece and comparing the distances to this first node
 * for the two placements of the two pieces
 */
void assembler_0_c::iterativeMultiSearch(void) {

  abbort = false;
  running = true;

  // this variable is used to shor if we continue with our loop over
  // the rows or have finished
  bool cont;

  while (!abbort) {

    // we have finished if pos negative (or greater than piecenumber because of the
    // overflow
    if (pos > piecenumber)
      break;

    // check, if all pieces are placed and all voxels are filled
    // careful here abbort is also modified by another thread
    // i once used the expression
    //   abbort |= !solution()
    // and search halve a day why it didn't work. The value of abbort was read
    // then the function called then the new value calculated then the new value
    // written. Meanwhile abbort was pressed and abbort was changed. this new value got lost.
    if (!right[0])
      solution();

    // the debugger
    if (debug) {
      if (debug_loops <= 0)
        break;

      debug_loops --;

      if (debug_pos == pos)
        break;
    }

    // if all pieces are placed we can not go on even if there are
    // more columns that need attention, all of them should be
    // empty any ways, so we backtrack once and continue there
    if (pos == piecenumber)
      pos--;

    cont = false;
    iterations++;

    if (!rows[pos]) {

      // start with a new column
      if (searchState[pos] == 0) {
        /* search the best column for the next recursion step
         *
         * we actually do a bit more:
         * we find the column with the smallest column count. This column is selected
         * from the piece columns and the normal result columns
         *
         * we also look for piece and result columns that have a count of 0 that value
         * will lead to impossible arrangements
         *
         * FIXME: possible improvement
         * For some columns is is necessary that they are not at least one but even
         * bigger. Columns for multi pieces for example need to have at least the value
         * of the number of instances, e.g if there are 10 instances of a multi piece
         * each of that instance needs to fit at least 10 times
         */
        unsigned int c = right[0];
        unsigned int s = colCount[c];

        if (s) {
          register unsigned int j = right[c];

          while (j) {

            if (colCount[j] < s) {
              c = j;
              s = colCount[c];

              if (!s)
                break;
            }

            // multi pieces need to have at least num placements
            // as we are not currently in the process of placing a multi
            // piece either all of a group of multipieces have been placed and
            // the columns removed or none of them. So we can check the columns
            // that are still there and force the counters above the number
            // of pieces in the group of multipieces
            if ((j <= piecenumber) && (colCount[j] < multiPieceCount[j-1]) && (multiPieceIndex[j-1] > 0)) {
              s = 0;
              break;
            }


            j = right[j];
          }
        }

        // now check for the holes, only the variable columns can be unfillable
        // but there must not be more unfillable voxels than there are holes
        // if that is the case we can backtrack
        // sometimes this doesnt help much, but it also seems like
        // it doesn't cost a lot of time, so let's keep it in for the moment
        if (s) {
          unsigned int currentHoles = holes;
          register unsigned int j = right[varivoxelEnd];

          while (j != varivoxelEnd) {
            if (colCount[j] == 0) {
              if (currentHoles == 0) {
                s = 0;
                break;
              }
              currentHoles--;
            }
            j = right[j];
          }
        }

        if (s) {

          // we have found a valid column, start a search
          columns[pos] = c;
          rows[pos] = down[columns[pos]];

          if (c <= piecenumber)
            piece[pos] = c-1;
          else {
            // find out the piece for the current column
            // the piece can only change for non piece columns
            piece[pos] = 0;
            while ((piece[pos] < piecenumber-1) && (rows[pos] >= pieceStart[piece[pos]+1])) piece[pos]++;
          }

          bt_assert(multiPieceIndex[piece[pos]] == 0);

          cont = true;
        }
      } else {
        // this is the branch for the case we are forced to place certain pieces
        // this happens, if we placed a piece with a shape used by more than one piece
        // we do not need to search for a good column here as we already know wich
        // column we are going to place, but we must check, if there is one column
        // with a count of 0 so we can back track
        unsigned int j = right[0];

        do {
          if (!colCount[j]) {
            break;
          }

          // check that multi pieces can still be places the required number of times. En detail:
          // if we have a piece column that is for the current multi piece group and the
          // counter is smaller than the number of pieces still required to place, break
          if ((piece[pos]+1 <= j-1) && (piece[pos]+numF[pos] > j-1) && (colCount[j] < numF[pos])) {
            break;
          }

          // if we have a column for a piece and that piece is NOT the currently placed multi-piece
          // the number of possible placements needs to be at least the size of the multi piece group
          if ((j <= piecenumber) && ((piece[pos] > j-1) || (piece[pos]+numF[pos] <= j-1)) && (colCount[j] < multiPieceCount[j-1])) {
            break;
          }

          j = right[j];
        } while (j);

        // FIXME: we should also do the hole check in here

        if (!j) {
          // select the column.
          // we have a fixed column given by the last recursion step, so take this column
          columns[pos] = piece[pos] + 1;
          rows[pos] = down[columns[pos]];
          cont = true;
        }
      }

      // found no fitting row, or column with zero count
      if (!cont) {
        rows[pos] = 0;
        pos--;
        continue;
      }

      cover(columns[pos]);

    } else {


      // continue on a column we have already started, this is inside the loop in the
      // recourseve function, after we return from the recoursive call
      // we uncover our row, find the next one and continue, if there is a new row
      uncover_row(rows[pos]);
      cont = true;

      do {
        rows[pos] = down[rows[pos]];

        if (rows[pos] == columns[pos]) {
          cont = false;
          break;
        }

        // find out the piece for the current column
        // the piece can only change for non piece columns
        if (columns[pos] > piecenumber)
          while ((piece[pos] < piecenumber-1) && (rows[pos] >= pieceStart[piece[pos]+1])) piece[pos]++;

        if (searchState[pos] == 1) break;

      } while (multiPieceIndex[piece[pos]] > 0);
    }

    // check, if we continue, if we have found a new row to cover
    if (cont) {

      // cover the row

      // the next lines with the crazy if-nesting sets up the parameters for the next
      // recoursion step, this is different depending on the state we are in if there
      // are more pieces in front or behind us...

      // now depending on the properties of the current piece set up the variables for the
      // next recoursion step.
      if ((searchState[pos] == 0) && (multiPieceCount[piece[pos]] != 1)) {

        // this is the case when we are in normal mode but want to place a multipiece

        bt_assert(multiPieceIndex[piece[pos]] == 0);

        // is is a multi piece, so we have to at least go forward until we have placed all
        // the multi pieces

        if (columns[pos] != piece[pos]+1)
          nodeF[pos+1] = 0;
        else
          nodeF[pos+1] = rows[pos] - pieceStart[piece[pos]];

        numF[pos+1] = multiPieceCount[piece[pos]] - 1 - multiPieceIndex[piece[pos]];
        piece[pos+1] = piece[pos] + 1;
        searchState[pos+1] = 1;

      } else if ((searchState[pos] == 1) && (numF[pos] > 1)) {

        // this is the case that we are in moving forward mode and
        // there are more pieces to place, so let's do another forward step

        nodeF[pos+1] = rows[pos]-pieceStart[piece[pos]];
        numF[pos+1] = numF[pos] - 1;
        piece[pos+1] = piece[pos] + 1;
        searchState[pos+1] = 1;

      } else {

        // nothing special continue in normal mode
        searchState[pos+1] = 0;
      }

      if ((searchState[pos+1] == 1) && (nodeF[pos+1] > 0)) {
        /* we cover all rows in all multi pieces that are still available
         * whose row number is forbidden
         * problem: how to uncover?
         * for this we have the stack addRows, that contains all the additionally covered
         * rows
         */

        // for all pieces of this multi-piece set
        for (unsigned int n = piece[pos+1]; n < piece[pos+1]+numF[pos+1]; n++) {
          unsigned int r = down[n+1];

          // check all the rows
          while (r != n+1) {

            // remove the row, if it is forbidden
            if ((r - pieceStart[n]) < nodeF[pos+1]) {
              remove_row(r);
              addRows[pos].push(r);
            }

            r = down[r];
          }
        }
      }

      cover_row(rows[pos]);

      pos++;

    } else {

      // now ucover the removed rows
      while (addRows[pos].size()) {
        reinsert_row(addRows[pos].top());
        addRows[pos].pop();
      }

      // ok finished this column, uncover it and backtrack
      uncover(columns[pos]);

      rows[pos] = 0;
      pos--;
    }
  }

  running = false;
}

void assembler_0_c::assemble(assembler_cb * callback) {

  debug = false;

  if (errorsState == ERR_NONE) {
    asm_bc = callback;
    iterativeMultiSearch();
  }
}

float assembler_0_c::getFinished(void) {

  /* we don't need locking, as I hope that I have written the
   * code in a way that updated the data so, that it will never
   * be in an inconsistant state. The thing that will happen is that
   * the value may jump
   */

  if (!rows || !columns || !down.size())
    return 0;

  float erg = 0;

  if (pos > piecenumber)
    return 1;

  for (int i = pos - 1; i >= 0; i--) {

    unsigned int r = rows[i];
    unsigned int l = colCount[columns[i]];

    while (l && r && (r != down[columns[i]])) {
      erg += 1;
      r = up[r];
      l--;
    }

    erg /= colCount[columns[i]];
  }

  return erg;
}

static unsigned int getInt(const char * s, unsigned int * i) {

  char * s2;

  *i = strtol (s, &s2, 10);

  if (s2)
    return s2-s;
  else
    return 500000;
}

static unsigned int getLong(const char * s, unsigned long * i) {

  char * s2;

  *i = strtol (s, &s2, 10);

  if (s2)
    return s2-s;
  else
    return 500000;
}

assembler_c::errState assembler_0_c::setPosition(const char * string, const char * version) {

  /* we assert that the matrix is in the initial position
   * otherwise we would have to clean the stack
   */
  bt_assert(pos == 0);

  /* check for te right version */
  if (strcmp(version, "1.3"))
    return ERR_CAN_NOT_RESTORE;

  /* get the values from the string. For the moment we assume that
   * the string is correct and contains enough values, if not sad things
   * will happen
   * FIXME
   */
  string += getInt(string, &pos);
  string += getLong(string, &iterations);

  if (pos <= piecenumber)
    for (unsigned int i = 0; i <= pos; i++) {
      while (*string != '(') string++; string++;

      string += getInt(string, &searchState[i]);

      switch(searchState[i]) {
      case 0:
        string += getInt(string, &rows[i]);
        string += getInt(string, &columns[i]);
        string += getInt(string, &piece[i]);
        break;
      case 1:
        string += getInt(string, &rows[i]);
        string += getInt(string, &columns[i]);
        string += getInt(string, &nodeF[i]);
        string += getInt(string, &numF[i]);
        string += getInt(string, &piece[i]);
        break;
      }

      while (*string != ')') string++; string++;
    }

  /* here we need to get the matrix into this exact position as it has been, when we
   * saved the position that means we need to cover all rows and columns in the same
   * order as it happened in the original process
   */
  unsigned int p = 0;

  if (pos <= piecenumber) {

    while (p <= pos) {

      if ((p < piecenumber) && rows[p]) {
        cover(columns[p]);

        if ((p < pos-1) && (searchState[p+1] == 1) && (nodeF[p+1] > 0)) {
          /* we cover all rows in all multi pieces that are still available
           * whose row number is forbidden
           * problem: how to uncover?
           * for this we have the stack addRows, that contains all the additionally covered
           * rows
           */

          // for all pieces of this multi-piece set
          for (unsigned int n = piece[p+1]; n < piece[p+1]+numF[p+1]; n++) {
            unsigned int r = down[n+1];

            // check all the rows
            while (r != n+1) {

              // remove the row, if it is forbidden
              if ((r - pieceStart[n]) < nodeF[p+1]) {
                remove_row(r);
                addRows[p].push(r);
              }

              r = down[r];
            }
          }
        }

        cover_row(rows[p]);
      }

      p++;
    }
  }

  return ERR_NONE;
}

xml::node assembler_0_c::save(void) const {

  xml::node nd("assembler");

  nd.get_attributes().insert("version", "1.2");

  std::string cont;

  char tmp[100];

  snprintf(tmp, 100, "%i ", pos);
  cont += tmp;

  snprintf(tmp, 100, "%li ", iterations);
  cont += tmp;

  if (pos <= piecenumber)
    for (unsigned int j = 0; j <= pos; j++) {

      switch(searchState[j]) {
      case 0:
        snprintf(tmp, 100, "(%i %i %i %i)", searchState[j], rows[j], columns[j], piece[j]);
        break;
      case 1:
        snprintf(tmp, 100, "(%i %i %i %i %i %i)", searchState[j], rows[j], columns[j], nodeF[j], numF[j], piece[j]);
        break;
      }

      cont += tmp;

      if (j < pos)
        cont += " ";
    }

  nd.set_content(cont.c_str());

  return nd;
}

unsigned int assembler_0_c::getPiecePlacement(unsigned int node, int delta, unsigned int piece, unsigned char *tran, int *x, int *y, int *z) {

  if (!node)
    node = down[piece + 1];

  while (delta > 0) {
    node = down[node];
    delta--;
  }

  while (delta < 0) {
    node = up[node];
    delta++;
  }

  getPieceInformation(node, tran, x, y, z);
  return node;
}

unsigned int assembler_0_c::getPiecePlacementCount(unsigned int piece) {

  return colCount[piece+1];
}


void assembler_0_c::debug_step(unsigned long num) {
  debug = true;
  debug_loops = 1;
  debug_pos = -1;
  asm_bc = 0;
  iterativeMultiSearch();
}

void assembler_0_c::debug_run(unsigned int level) {
}

