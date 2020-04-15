#include "ift.h"

float max(float a, float b) {
  return a > b? a : b;
}

float absolute(float a) {
  return a > 0? a: -1 * a;
}

iftImage *GetSliceXY(iftImage *img, int z)
{
  iftImage *slc = iftCreateImage(img->xsize,img->ysize,1);
  iftVoxel  u;
  int       q=0;
  
  u.z = z;
  for (u.y = 0; u.y < img->ysize; u.y++)
    for (u.x = 0; u.x < img->xsize; u.x++){
      int p = iftGetVoxelIndex(img,u);      
      slc->val[q] = img->val[p];
      q++;
    }

  return(slc);
}

iftImage *GetSliceSagital(iftImage *img, int x, short vp)
{
  iftImage *slc = iftCreateImage(img->ysize,img->zsize,1);
  iftVoxel  u;
  int       q=0;
  
  u.x = x;
  if (vp == 0) {
    for (u.z = img->zsize-1; u.z >= 0; u.z--)
      for (u.y = 0; u.y < img->ysize; u.y++){
        int p = iftGetVoxelIndex(img,u);      
        slc->val[q] = img->val[p];
        q++;
      }
  } else {
    for (u.z = img->zsize-1; u.z >= 0; u.z--)
      for (u.y = img->ysize-1; u.y >= 0; u.y--){
        int p = iftGetVoxelIndex(img,u);      
        slc->val[q] = img->val[p];
        q++;
      }
  }

  return(slc);
}

iftImage *GetSliceCoronal(iftImage *img, int y, short vp)
{
  iftImage *slc = iftCreateImage(img->xsize,img->zsize,1);
  iftVoxel  u;
  int       q=0;
  
  u.y = y;
  if (vp == 0) {
    for (u.z = img->zsize-1; u.z >= 0; u.z--)
      for (u.x = 0; u.x < img->xsize; u.x++){
        int p = iftGetVoxelIndex(img,u);      
        slc->val[q] = img->val[p];
        q++;
      }
    
  } else {
    for (u.z = img->zsize-1; u.z >= 0; u.z--)
      for (u.x = img->xsize-1; u.x >= 0; u.x--){
        int p = iftGetVoxelIndex(img,u);      
        slc->val[q] = img->val[p];
        q++;
      }
  }

  return(slc);
}

iftImage *GetSliceAxial(iftImage *img, int z, short vp)
{
  iftImage *slc = iftCreateImage(img->xsize,img->ysize,1);
  iftVoxel  u;
  int       q=0;
  
  u.z = z;
  if (vp == 0) {
    for (u.y = 0; u.y < img->ysize; u.y++)
    for (u.x = 0; u.x < img->xsize; u.x++){
      int p = iftGetVoxelIndex(img,u);      
      slc->val[q] = img->val[p];
      q++;
    }
  } else {
    for (u.y = img->ysize-1; u.y >= 0; u.y--)
      for (u.x = img->xsize-1; u.x >= 0; u.x--){
        int p = iftGetVoxelIndex(img,u);      
        slc->val[q] = img->val[p];
        q++;
      }
  }

  return(slc);
}


iftImage *applyRainBowColorTable(iftImage *img, int h){
  iftImage *colored = iftCreateColorImage(img->xsize,img->ysize, img->zsize, 16);

  iftVoxel u;
 
  for(u.x = 0; u.x < img->xsize; u.x++) {
    for(u.y = 0; u.y < img->ysize; u.y++) {
      for(u.z = 0; u.z < img->zsize; u.z++) {
        int p = iftGetVoxelIndex(img, u);
        float v = (float) img->val[p];
        
        v = v/h;
        v = 4*v + 1;

        int r = iftRound(h * max(0.0, (3 - absolute(v -4) - absolute(v -5))/2));
        int g = iftRound(h * max(0.0, (4 - absolute(v -2) - absolute(v -4))/2));
        int b = iftRound(h * max(0.0, (3 - absolute(v -1) - absolute(v -2))/2));

        iftColor rgb_color;
        rgb_color.val[0] = r;
        rgb_color.val[1] = g;
        rgb_color.val[2] = b;
        //rgb_color.alpha = 1;

        iftColor YCbCr_color = iftRGBtoYCbCr(rgb_color, h);
      
        colored->val[p] = YCbCr_color.val[0];
        colored->Cb[p] = YCbCr_color.val[1];
        colored->Cr[p] = YCbCr_color.val[2];

      }
    }
  }

  return colored;
}

iftImage *window_level(iftImage *img, int level, int midpoint, int h){
  int l1 = iftRound(midpoint - level/2);
  int l2 = iftRound(level/2 + midpoint);

  iftImage *streached = iftCreateImage(img->xsize,img->ysize, img->zsize);

  iftVoxel u;

  for(u.x = 0; u.x < img->xsize; u.x++) {
    for(u.y = 0; u.y < img->ysize; u.y++) {
      for(u.z = 0; u.z < img->zsize; u.z++) {
        int p = iftGetVoxelIndex(img, u);
        int l = img->val[p];
        int k = 0;
        if (l > l2) {
          k = h;
        } else if (l < l1) {
          k = 0;
        } else {
          k = iftRound(h/(l2 - l1)*(l - l1));
        }
        streached->val[p] = k;
      }
    }
  }

  return streached;

}

int main(int argc, char *argv[]) 
{
  timer *tstart = NULL;
  int    MemDinInicial, MemDinFinal;
  
  MemDinInicial = iftMemoryUsed(1);

  if (argc != 8){
    iftError("Usage: getslice <...>\n"
	     "[1] input image .scn \n"
	     "[2] x coordinate from which the code will extract a sagital slice \n"
       "[3] y coordinate from which the code will extract a coronal slice \n"
       "[4] z coordinate from which the code will extract a axial slice \n"
       "[5] view n 0/1 for the point of view of the radiologists/neuroradiologists \n"
       "[6] windows size for linear stretching \n"
       "[7] level for linear stretching \n",
	     "main");
  }

  tstart = iftTic();

  /* ----------------------- Coding Area -------------------------- */

  iftImage *img    = iftReadImageByExt(argv[1]);
  int x            = atoi(argv[2]);
  int y            = atoi(argv[3]);
  int z            = atoi(argv[4]);
  short vp         = atoi(argv[5]);
  double window_p  = atof(argv[6]);
  double level_p   = atof(argv[7]);

  if ((x < 0)||(x >= img->xsize))
    iftError("x-coordinate must be in [0,%d]", "main", img->xsize-1);
  
  if ((y < 0)||(y >= img->ysize))
    iftError("y-coordinate must be in [0,%d]", "main", img->ysize-1);

  if ((z < 0)||(z >= img->zsize))
    iftError("z-coordinate must be in [0,%d]", "main", img->zsize-1);

  if (vp != 0 && vp != 1)
    iftError("view point must be 0 or 1", "main");

  if (window_p < 0 || window_p > 1)
    iftError("window size must be in [0, 1]", "main");

  if (level_p < 0 || level_p > 1)
    iftError("window size must be in [0, 1]", "main");

  int img_min = iftMinimumValue(img);
  int img_max = iftMaximumValue(img);

  int window = iftRound((img_max - img_min)*window_p);
  int level = iftRound(img_max * level_p);
  
  int h = 256*256-1;

  iftImage *norm = window_level(img, window, level, h);
  
  iftImage *slc_sagital = GetSliceSagital(norm, x, vp);
  iftImage *colored_slc_sagital = applyRainBowColorTable(slc_sagital, h);

  if (vp == 0) {
    iftWriteImageByExt(slc_sagital, "radiologist-sagital-gray.pgm");
    iftWriteImageByExt(colored_slc_sagital, "radiologist-sagital.png");
  } else {
    iftWriteImageByExt(slc_sagital, "_neuroradiologist-sagital-gray.pgm");
    iftWriteImageByExt(colored_slc_sagital, "neuroradiologist-sagital.png");
  }

  iftImage *slc_coronal = GetSliceCoronal(norm, y, vp); 
  iftImage *colored_slc_coronal = applyRainBowColorTable(slc_coronal, h);

  if (vp == 0) {
    iftWriteImageByExt(slc_coronal, "radiologist-coronal-gray.pgm");
    iftWriteImageByExt(colored_slc_coronal, "radiologist-coronal.png");
  } else {
    iftWriteImageByExt(slc_coronal, "neuroradiologist-coronal-gray.pgm");
    iftWriteImageByExt(colored_slc_coronal, "neuroradiologist-coronal.png");
  }

  iftImage *slc_axial = GetSliceAxial(norm, z, vp);
  iftImage *colored_slc_axial = applyRainBowColorTable(slc_axial, h);

  if (vp == 0) {
    iftWriteImageByExt(slc_axial, "radiologist-axial-gray.pgm");
    iftWriteImageByExt(colored_slc_axial, "radiologist-axial.png");
  } else {
    iftWriteImageByExt(slc_axial, "neuroradiologist-axial-gray.pgm");
    iftWriteImageByExt(colored_slc_axial, "neuroradiologist-axial.png");
  }

  iftDestroyImage(&slc_sagital);
  iftDestroyImage(&slc_coronal);
  iftDestroyImage(&slc_axial);
  iftDestroyImage(&colored_slc_sagital);
  iftDestroyImage(&colored_slc_coronal);
  iftDestroyImage(&colored_slc_axial);
  iftDestroyImage(&img);
  iftDestroyImage(&norm);
  
  /* -------------------- End of the coding area ----------------- */
    
  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));
    
  MemDinFinal = iftMemoryUsed();
  iftVerifyMemory(MemDinInicial, MemDinFinal);
  
  return(0);
}
