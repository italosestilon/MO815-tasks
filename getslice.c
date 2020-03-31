#include "ift.h"

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
int main(int argc, char *argv[]) 
{
  timer *tstart = NULL;
  int    MemDinInicial, MemDinFinal;
  
  MemDinInicial = iftMemoryUsed(1);

  if (argc != 4){
    iftError("Usage: getslice <...>\n"
	     "[1] input image .scn \n"
	     "[2] z-coordinate of a xy slice \n"
	     "[3] output image .pgm \n",
	     "main");
  }

  tstart = iftTic();

  /* ----------------------- Coding Area -------------------------- */

  iftImage *img    = iftReadImageByExt(argv[1]);
  int z            = atoi(argv[2]);
  if ((z < 0)||(z >= img->zsize))
    iftError("z-coordinate must be in [0,%d]", "main", img->zsize-1);
  
  iftImage *slc    = GetSliceXY(img, atoi(argv[2])); 
  iftWriteImageByExt(slc,argv[3]);

  iftDestroyImage(&slc);
  iftDestroyImage(&img);
  
  /* -------------------- End of the coding area ----------------- */
    
  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));
    
  MemDinFinal = iftMemoryUsed();
  iftVerifyMemory(MemDinInicial, MemDinFinal);
  
  return(0);
}
