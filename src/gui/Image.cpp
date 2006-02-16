#include "Image.h"

#include "../lib/bt_assert.h"

#include <png.h>
#include <stdio.h>

Image::Image(unsigned int w, unsigned int h, unsigned char r, unsigned char g, unsigned char b, unsigned char a) : width(w), height(h) {
  bitmap = new unsigned char[w*h*4];
  for (unsigned int x = 0; x < w*h; x++) {
    bitmap[4*x+0] = r;
    bitmap[4*x+1] = g;
    bitmap[4*x+2] = b;
    bitmap[4*x+3] = a;
  }
}


Image::Image(unsigned int w, unsigned int h, unsigned char *b) : width(w), height(h), bitmap(b) {
}

Image::Image(unsigned int w, unsigned int h, VoxelDrawer * dr, const VoxelView *v3d) : width(w), height(h) {
  bitmap = new unsigned char[w*h*4];

  TRcontext *tr = trNew();

  trTileSize(tr, v3d->w(), v3d->h(), 0);
  trPerspective(tr, 5 + v3d->getSize(), 1.0, 10, 1100);

  trImageSize(tr, w, h);
  trImageBuffer(tr, GL_RGBA, GL_UNSIGNED_BYTE, bitmap);

  GLfloat LightAmbient[]= { 0.01f, 0.01f, 0.01f, 1.0f };
  GLfloat LightDiffuse[]= { 1.5f, 1.5f, 1.5f, 1.0f };
  GLfloat LightPosition[]= { 700.0f, 200.0f, -90.0f, 1.0f };

  GLfloat AmbientParams[] = {0.1, 0.1, 0.1, 1};
  GLfloat DiffuseParams[] = {0.7, 0.7, 0.7, 0.1};
  GLfloat SpecularParams[] = {0.4, 0.4, 0.4, 0.5};

  glEnable(GL_COLOR_MATERIAL);
  glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
  glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);
  glEnable(GL_LIGHT1);

  glMaterialf(GL_FRONT, GL_SHININESS, 0.5);
  glMaterialfv(GL_FRONT, GL_AMBIENT, AmbientParams);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, DiffuseParams);
  glMaterialfv(GL_FRONT, GL_SPECULAR, SpecularParams);

  glClearColor(1, 1, 1, 0);

  do {
    trBeginTile(tr);
    dr->drawData();

  } while (trEndTile(tr));

  // now flip vertically
  for (unsigned int y = 0; y < h/2; y++)
    for (unsigned int x = 0; x < w*4; x++) {
      unsigned char tmp = bitmap[y*4*w+x];
      bitmap[y*4*w+x] = bitmap[(h-y-1)*4*w+x];
      bitmap[(h-y-1)*4*w+x] = tmp;
    }

  trDelete(tr);
}


Image::~Image(void) {
  delete [] bitmap;
}

int Image::saveToPNG(const char * fname) const {

  int sx = width;
  int sy = height;
  unsigned char * buffer = bitmap;

  png_structp png_ptr;
  png_infop info_ptr;
  unsigned char ** png_rows;
  int x, y;

  FILE *fi = fopen(fname, "wb");
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (png_ptr == NULL)
  {
    fclose(fi);
    png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
    fprintf(stderr, "\nError: Couldn't save the image!\n%s\n\n", fname);
    return 0;
  }

  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL)
  {
    fclose(fi);
    png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
    fprintf(stderr, "\nError: Couldn't save the image!\n%s\n\n", fname);
    return 0;
  }

  if (setjmp(png_jmpbuf(png_ptr)))
  {
    fclose(fi);
    png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
    fprintf(stderr, "\nError: Couldn't save the image!\n%s\n\n", fname);
    return 0;
  }

  png_init_io(png_ptr, fi);
  info_ptr->width = sx;
  info_ptr->height = sy;
  info_ptr->bit_depth = 8;
  info_ptr->color_type = PNG_COLOR_TYPE_RGB_ALPHA;
  info_ptr->interlace_type = PNG_INTERLACE_NONE;
  info_ptr->valid = 0;

  png_write_info(png_ptr, info_ptr);

  /* Save the picture: */

  png_rows = new unsigned char*[sy]; //     (unsigned char **)(malloc(sizeof(char *) * sy));
  for (y = 0; y < sy; y++)
  {
    png_rows[y] = new unsigned char[4*sx];// (unsigned char*)malloc(sizeof(char) * 3 * sx);

    for (x = 0; x < sx; x++)
    {
      png_rows[y][x * 4 + 0] = buffer[(y*sx+x)*4+0];
      png_rows[y][x * 4 + 1] = buffer[(y*sx+x)*4+1];
      png_rows[y][x * 4 + 2] = buffer[(y*sx+x)*4+2];
      png_rows[y][x * 4 + 3] = buffer[(y*sx+x)*4+3];
    }
  }

  png_write_image(png_ptr, png_rows);

  for (y = 0; y < sy; y++)
    delete [] png_rows[y]; //free(png_rows[y]);

  delete [] png_rows; //free(png_rows);

  png_write_end(png_ptr, NULL);

  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fi);

  return 1;
}

void Image::blit(const Image * i, int xpos, int ypos) {
  for (unsigned int x = 0; x < i->width; x++)
    for (unsigned int y = 0; y < i->height; y++)
      if ((x+xpos >= 0) && (x+xpos < width) && (y+ypos >= 0) && (y+ypos < height)) {
        unsigned char r1,r2, g1, g2, b1, b2, a1, a2;

        r1 = bitmap[((y+ypos)*width + (x+xpos)) * 4 + 0];
        r2 = i->bitmap[(y*i->width+x) * 4 + 0];

        g1 = bitmap[((y+ypos)*width + (x+xpos)) * 4 + 1];
        g2 = i->bitmap[(y*i->width+x) * 4 + 1];

        b1 = bitmap[((y+ypos)*width + (x+xpos)) * 4 + 2];
        b2 = i->bitmap[(y*i->width+x) * 4 + 2];

        a1 = bitmap[((y+ypos)*width + (x+xpos)) * 4 + 3];
        a2 = i->bitmap[(y*i->width+x) * 4 + 3];

        bitmap[((y+ypos)*width + (x+xpos)) * 4 + 0] = (r1 * (255-a2) + r2 * a2) / 255;
        bitmap[((y+ypos)*width + (x+xpos)) * 4 + 1] = (g1 * (255-a2) + g2 * a2) / 255;
        bitmap[((y+ypos)*width + (x+xpos)) * 4 + 2] = (b1 * (255-a2) + b2 * a2) / 255;

        bitmap[((y+ypos)*width + (x+xpos)) * 4 + 3] = (unsigned char)((1 - ((1-a1/255.0) * (1-a2/255.0))) * 255 + 0.5);
      }
}

void Image::deTransparentize(unsigned char r, unsigned char g, unsigned char b) {
  for (unsigned int x = 0; x < width; x++)
    for (unsigned int y = 0; y < height; y++)
      if (bitmap[4*(y*width + x) + 3] == 0) {
        bitmap[4*(y*width + x) + 0] = r;
        bitmap[4*(y*width + x) + 1] = g;
        bitmap[4*(y*width + x) + 2] = b;
        bitmap[4*(y*width + x) + 3] = 255;
      }
}

void Image::scaleDown(unsigned char by) {
  bt_assert(width % by == 0);
  bt_assert(height % by == 0);

  unsigned int nw = width / by;
  unsigned int nh = height / by;

  unsigned char * nb = new unsigned char [nw*nh*4];

  for (unsigned int x = 0; x < nw; x++)
    for (unsigned int y = 0; y < nh; y++) {
      unsigned int a = 0;
      unsigned int r = 0;
      unsigned int g = 0;
      unsigned int b = 0;

      for (unsigned int ax = 0; ax < by; ax++)
        for (unsigned int ay = 0; ay < by; ay++) {
          unsigned int pos = 4*((y*by+ay)*width + x*by+ax);
          if (bitmap[pos + 3]) {
            r += bitmap[pos + 0];
            g += bitmap[pos + 1];
            b += bitmap[pos + 2];
            a += bitmap[pos + 3];
          }
        }

      if (a) {
        nb[4*(y*nw+x) + 0] = 255 * r / a;
        nb[4*(y*nw+x) + 1] = 255 * g / a;
        nb[4*(y*nw+x) + 2] = 255 * b / a;
      } else {
        nb[4*(y*nw+x) + 0] = 0;
        nb[4*(y*nw+x) + 1] = 0;
        nb[4*(y*nw+x) + 2] = 0;
      }
      nb[4*(y*nw+x) + 3] = a / (by*by);
    }

  delete [] bitmap;
  bitmap = nb;

  height = nh;
  width = nw;
}

void Image::minimizeWidth(unsigned int border, unsigned int multiple) {

  unsigned int xmin = 0;

  // find the smallest used column
  do {

    bool allTran = true;
    // check column for a pixed that is not completely transparent
    for (unsigned int y = 0; y < height; y++)
      if (bitmap[y*width*4+4*xmin+3] != 0) {
        allTran = false;
        break;
      }

    if (!allTran)
      break;

    xmin++;

  } while (xmin < width);

  unsigned int xmax = width-1;

  // find the smallest used column
  do {

    bool allTran = true;
    // check column for a pixed that is not completely transparent
    for (unsigned int y = 0; y < height; y++)
      if (bitmap[y*width*4+4*xmax+3] != 0) {
        allTran = false;
        break;
      }

    if (!allTran)
      break;

    xmax--;

  } while (xmax > 0);

  // include the border
  if (xmin > border)
    xmin -= border;
  else
    xmin = 0;

  if (xmax + border + 1 < width)
    xmax += border;
  else
    xmax = width-1;

  if (xmin > xmax) {
    xmin = xmax = 0;
  }

  unsigned int nw = xmax-xmin+1;

  if ((nw % multiple) != 0) {
    unsigned int add = multiple - (nw % multiple);

    if (xmin >= add/2) {
      xmin -= add/2;
      add -= add/2;
    } else {
      add -= xmin;
      xmin = 0;
    }

    xmax += add;

    nw = xmax-xmin+1;
  }

  // create new bitmap
  unsigned char * nb = new unsigned char[nw*height*4];

  // copy image information
  for (unsigned int y = 0; y < height; y++)
    for (unsigned int x = xmin; x <= xmax; x++) {
      nb[(y*nw+x-xmin)*4+0] = bitmap[(y*width+x)*4+0];
      nb[(y*nw+x-xmin)*4+1] = bitmap[(y*width+x)*4+1];
      nb[(y*nw+x-xmin)*4+2] = bitmap[(y*width+x)*4+2];
      nb[(y*nw+x-xmin)*4+3] = bitmap[(y*width+x)*4+3];
    }

  delete [] bitmap;
  bitmap = nb;
  width = nw;
}

