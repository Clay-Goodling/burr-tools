/* Burr Solver
 * Copyright (C) 2003-2005  Andreas R�ver
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


#include "AssemblyCallbacks.h"

#ifdef WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

#include "../lib/disassembly.h"
#include "../lib/assm_0_frontend_0.h"
#include "../lib/disassembler_3.h"

#ifdef WIN32
unsigned long __stdcall start_th(void * c)
#else
void* start_th(void * c)
#endif
{
  assemblerThread * p = (assemblerThread*)c;

  assm_0_frontend_0_c * assm;

  if (!p->puzzle->probGetAssembler(p->prob)) {
    p->action = assemblerThread::ACT_PREPARATION;
    assm = new assm_0_frontend_0_c();
    p->puzzle->probSetAssembler(p->prob, assm);

    p->errState = assm->createMatrix(p->puzzle, p->prob);
    if (p->errState != assm_0_frontend_0_c::ERR_NONE) {

      p->errParam = assm->getErrorsParam();

      p->action = assemblerThread::ACT_ERROR;
      return 0;
    }

    if (!p->stopPressed)
      p->action = assemblerThread::ACT_REDUCE;

    assm->reduce();

  } else
    assm = (assm_0_frontend_0_c*)p->puzzle->probGetAssembler(p->prob);

  if (!p->stopPressed) {

    p->action = assemblerThread::ACT_ASSEMBLING;
    assm->assemble(p);
  
    if (assm->getFinished() >= 1)
      p->action = assemblerThread::ACT_FINISHED;
    else
      p->action = assemblerThread::ACT_PAUSING;

  } else
    p->action = assemblerThread::ACT_PAUSING;

  return 0;
}

assemblerThread::assemblerThread(puzzle_c * puz, unsigned int problemNum, unsigned int solAction) :
assemblies(0),
action(ACT_PREPARATION),
_solutionAction(solAction),
puzzle(puz),
prob(problemNum)
{
}

assemblerThread::~assemblerThread(void) {

  puzzle->probGetAssembler(prob)->stop();

  while (!puzzle->probGetAssembler(prob)->stopped())
#ifdef WIN32
    Sleep(1);
#else
    usleep(10000);
#endif

}

bool assemblerThread::assembly(assembly_c * a) {

  assemblies++;

  switch(_solutionAction) {
  case SOL_SAVE_ASM:
    puzzle->probAddSolution(prob, a);
    break;

  case SOL_DISASM:
    {
      action = ACT_DISASSEMBLING;

      assemblyVoxel_c *as = a->getVoxelSpace(puzzle, prob);
      disassembler_3_c d(as, puzzle->probPieceNumber(prob));
      separation_c * s = d.disassemble();
      action = ACT_ASSEMBLING;
  
      if (s)
        puzzle->probAddSolution(prob, a, s);

      delete as;
    }
    break;
  }

  return true;
}

void assemblerThread::stop(void) {
  action = ACT_WAIT_TO_STOP;
  puzzle->probGetAssembler(prob)->stop();
  stopPressed = true;
}

void assemblerThread::start(void) {

  stopPressed = false;

#ifdef WIN32
  CreateThread(NULL, 0, start_th, this, 0, 0);
#else
  pthread_t th; pthread_create(&th, 0, start_th, this);
#endif

}

