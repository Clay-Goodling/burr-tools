/* Burr Solver
 * Copyright (C) 2000-2005  Andreas R�ver
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

#include "configuration.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <FL/Fl_Pack.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl.H>

// retugns either the home directory or the current
// directory on system that don't know home dirs
static char * homedir() {

#if (SYSTEM == SYS_LINUX || SYSTEM == SYS_MACOSX)

  return getenv("HOME");

#else

  return "./";

#endif

}



static FILE *create_local_config_file(void) {

#if (SYSTEM == SYS_LINUX || SYSTEM == SYS_MACOSX)

  char n[200];
  snprintf(n, 199, "%s/.burrtools.rc", homedir());
  return fopen(n, "w");

#else

  return fopen("burrtools.conf", "w");

#endif

}

static FILE *open_local_config_file(void) {

#if (SYSTEM == SYS_LINUX || SYSTEM == SYS_MACOSX)

  char n[200];
  snprintf(n, 199, "%s/.burrtools.rc", homedir());
  return fopen(n, "r");

#else

  return fopen("burrtools.conf", "r");

#endif

}




static bool str2bool(char *s) {
  if (s) {
    if (!strcmp("yes", s) || !strcmp("true", s)) return true;
    else return (atoi(s) != 0);
  }
  return false;
}

void configuration::parse(FILE * in) {
  char line[201], param[201];

  while (!feof(in)) {
    if (fscanf(in, "%200s%*[\t ]%200s\n", line, param) == 2) {

      config_data *t = first_data;

      while (t) {
        if (strstr(line, t->cnf_name)) {
          switch (t->cnf_typ) {
          case CT_BOOL:
            *(bool *)t->cnf_var = str2bool(param);
            break;
          case CT_STRING:
            if (strlen(param) > 1) {
              param[strlen(param)-1] = '\0';
              strncpy((char *)t->cnf_var, param+1, t->maxlen);
            }
            break;
          case CT_INT:
            *(int *)t->cnf_var = atoi(param);
            break;
          default: assert(0);
          }
          break;
        }
        t = t->next;
      }
    }
  }
}

void configuration::register_entry(char *cnf_name, cnf_type  cnf_typ, void *cnf_var, long maxlen) {
  config_data *t = new config_data;

  t->next = first_data;
  first_data = t;

  t->cnf_name = cnf_name;
  t->cnf_typ = cnf_typ;
  t->cnf_var = cnf_var;
  t->maxlen = maxlen;
}

#define CNF_BOOL(a,b) register_entry(a, CT_BOOL, b, 0)
#define CNF_CHAR(a,b,c) register_entry(a, CT_STRING, b, c)
#define CNF_INT(a,b) register_entry(a, CT_INT, b, 0)

configuration::configuration(void) {

  i_use_tooltips = true;
  i_use_lightning = true;
  first_data = 0;

  CNF_BOOL("tooltips",            &i_use_tooltips);
  CNF_BOOL("lightning",           &i_use_lightning);

  FILE * f = open_local_config_file();
  need_save = (f == 0);
  if (f) {
    parse(f);
    fclose(f);
  }
}


configuration::~configuration(void) {

  if (need_save) {

    FILE * f = create_local_config_file();

    fseek(f, 0, SEEK_SET);

    config_data *t = first_data;

    while(t) {
      fprintf(f, "%s: ", t->cnf_name);

      switch (t->cnf_typ) {
      case CT_BOOL:
        fprintf(f, "%s", (*(bool *)t->cnf_var)?("yes"):("no"));
        break;
      case CT_STRING:
        fprintf(f, "\"%s\"", (char *)(t->cnf_var));
        break;
      case CT_INT:
        fprintf(f, "%i", *(int *)t->cnf_var);
        break;
      default: assert(0);
      }

      fprintf(f, "\n");

      t = t->next;
    }

    fclose(f);
  }

  config_data *t = first_data;

  while (t) {
    t = t->next;
    delete first_data;
    first_data = t;
  }
}


void configuration::dialog(void) {

  Fl_Double_Window * win = new Fl_Double_Window(200, 500);
  Fl_Pack * pack = new Fl_Pack(0, 0, 200, 500);


  config_data *t = first_data;
  int y = 0;

  while(t) {
    switch (t->cnf_typ) {
    case CT_BOOL:
      {
          Fl_Check_Button *w = new Fl_Check_Button(0, y, 200, 20, t->cnf_name);
        t->widget = w;
        if (*((bool*)t->cnf_var))
          w->value(1);
        else
          w->value(0);
      }
      break;
    case CT_STRING:
      break;
    case CT_INT:
      break;
    default: assert(0);
    }

    y += 20;

    t = t->next;
  }

  pack->end();
  win->end();

  win->show();

  while (win->visible())
    Fl::wait();

  t = first_data;

  while(t) {
    switch (t->cnf_typ) {
    case CT_BOOL:

      if (((Fl_Check_Button*)t->widget)->value())
        *((bool*)t->cnf_var) = true;
      else
        *((bool*)t->cnf_var) = false;;
      break;
    case CT_STRING:
      break;
    case CT_INT:
      break;
    default: assert(0);
    }
    t = t->next;
  }

  need_save = true;

  delete win;
}
configuration config;