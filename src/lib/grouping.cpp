#include "grouping.h"
#include <assert.h>


void grouping_c::addPieces(unsigned int pc, unsigned int group, unsigned int count) {

  if (group == 0)
    return;

  unsigned int i = findPiece(pc, group);

  if (i < pieces.size())
    pieces[i].count += count;

  else {

    struct piece s;

    s.piece = pc;
    s.group = group;
    s.count = count;

    pieces.push_back(s);

    if (group+1 > numGroups) numGroups = group+1;
  }
}

void grouping_c::reSet(void) {
  sets.clear();
}


void grouping_c::newSet(void) {

  struct set s;

  s.currentGroup = 1;

  sets.push_back(s);
}



bool grouping_c::addPieceToSet(unsigned int pc) {

  assert(sets.size() > 0);

  unsigned int set = sets.size()-1;

  // add the piece to the current set
  sets[set].pieces.push_back(pc);

  // try to add the set to the pieces
  unsigned int i = findPiece(pc, sets[set].currentGroup);

  // check, if the required group piece combination is available
  if ((i < pieces.size()) && (pieces[i].count > 0)) {

    pieces[i].count--;
    return true;

  } else {
    // not available -> remove set from group
    for (unsigned int p = 0; p < sets[set].pieces.size()-1; p++)
      pieces[findPiece(sets[set].pieces[p], sets[set].currentGroup)].count++;
  }

  // ok when we get here something doesn't fit, so
  // take the last set and increase the group
  // if we reaced the last group restart with the
  // first and increment the 2nd last group, until
  // the first group reaced the end
  // as soon as we find something suitable in here, stop

  // now increment the groups until we either find a fitting
  // assignment or tried everything
  do {

    // increment current Group
    sets[set].currentGroup++;

    // all groups have been tried for this set
    if (sets[set].currentGroup >= numGroups) {
      sets[set].currentGroup == 0;

      if (set == 0)
        return false;

      set--;

      // remove pieces of the current set from the current group
      for (unsigned int p = 0; p < sets[set].pieces.size(); p++)
        pieces[findPiece(sets[set].pieces[p], sets[set].currentGroup)].count++;

    } else {

      // try to place the pieces of the set into the current group

      for (unsigned int p = 0; p < sets[set].pieces.size(); p++) {
        unsigned int i = findPiece(sets[set].pieces[p], sets[set].currentGroup);

        // if the required group doesn't exists or not enough pieces are in there
        // try again with next group for set
        if ((i < pieces.size()) && (pieces[i].count > 0))
          pieces[i].count--;

        else {

          // remove already placed pieces from the group
          for (unsigned int p2 = 0; p2 < p; p++)
            pieces[findPiece(sets[set].pieces[p2], sets[set].currentGroup)].count++;

          set--;
          break;

        }
      }
      set++;
    }

  } while (set < sets.size());

  return true;
}
