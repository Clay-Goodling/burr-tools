#ifndef __BURRGROWER_H__
#define __BURRGROWER_H__

#include "puzzle.h"
#include "assembler.h"

#include <vector>

class puzzleSol_c : public assembler_cb {

private:

  puzzle_c * puzzle;

  unsigned long solutions;
  unsigned long maxMoves;
  unsigned long minMoves;
  unsigned long maxLevel;
  unsigned long minLevel;

public:

  puzzleSol_c(puzzle_c * p);
  puzzleSol_c(const puzzleSol_c * p);
  ~puzzleSol_c(void);

  double fitness(void);

  bool assembly(voxel_c * assm);

  bool nosol(void) { return solutions == 0; }

  const puzzle_c * getPuzzle(void) const { return puzzle; }

  unsigned long numSolutions(void) { return solutions; }
  unsigned long numMoves(void) { return maxMoves; }
  unsigned long numLevel(void) { return maxLevel; }

};



class burrGrower_c {

  const puzzle_c * base;
  unsigned int maxSetSize;

  // puzzles with unique solutions and highes first level
  std::vector<puzzleSol_c*> unique;

  // highest first level
  std::vector<puzzleSol_c*> highLevel;

  // highest overal moves
  std::vector<puzzleSol_c*> highMoves;

  // this function checks every found puzzle for some conditions
  // and saves the most interesting designs found
  void addToLists(puzzleSol_c * pz);

public:

  burrGrower_c(const puzzle_c *pz, unsigned int mss) : base(pz), maxSetSize(mss) {}

  void grow(std::vector<puzzleSol_c*> currentSet);
};

#endif
