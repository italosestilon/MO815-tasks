#include "ift.h"

typedef struct _face {
  iftVector normal;
  iftPoint center;
} Face;

typedef struct _graphcontext{
  Face face[6];
  float alpha, beta;
  iftMatrix *Rx, *Ry, *Txyz, *Tuv, *Phi, *Phi_r;
  float diag;
} GraphicalContext;

iftPoint addPoints(iftPoint a, iftPoint b) {
  iftPoint c = {a.x+b.x, a.y + b.y, a.z + b.z};

  return c;
}
float absolute(float a) {
  return a > 0? a: -1 * a;
}

iftPoint mulByScalar(iftPoint a, float s) {
  iftPoint c = {a.x * s, a.y * s, a.z * s};
  return c;
}

int validPoint(iftImage *img, iftPoint p){
  return p.z >= 0 && p.z < img->zsize;
}

float PointsDotProd(iftPoint a, iftPoint b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

int isPointsEqual(iftPoint a, iftPoint b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

GraphicalContext *create_graphical_context(iftImage *img, float alpha, float beta) {
    GraphicalContext *gc = calloc(1, sizeof(GraphicalContext));

    int xsize = img->xsize;
    int ysize = img->ysize;
    int zsize = img->zsize;

    gc->face[0].normal.x =  1; gc->face[0].normal.y =  0; gc->face[0].normal.z =  0; 
    gc->face[1].normal.x = -1; gc->face[1].normal.y =  0; gc->face[1].normal.z =  0; 
    gc->face[2].normal.x =  0; gc->face[2].normal.y =  1; gc->face[2].normal.z =  0; 
    gc->face[3].normal.x =  0; gc->face[3].normal.y = -1; gc->face[3].normal.z =  0; 
    gc->face[4].normal.x =  0; gc->face[4].normal.y =  0; gc->face[4].normal.z =  1; 
    gc->face[5].normal.x =  0; gc->face[5].normal.y =  0; gc->face[5].normal.z = -1; 
    gc->face[0].center.x =  xsize-1; gc->face[0].center.y = ysize/2; gc->face[0].center.z = zsize/2; 
    gc->face[1].center.x =        0; gc->face[1].center.y = ysize/2; gc->face[1].center.z = zsize/2; 
    gc->face[2].center.x =  xsize/2; gc->face[2].center.y = ysize-1; gc->face[2].center.z = zsize/2; 
    gc->face[3].center.x =  xsize/2; gc->face[3].center.y =       0; gc->face[3].center.z = zsize/2; 
    gc->face[4].center.x =  xsize/2; gc->face[4].center.y = ysize/2; gc->face[4].center.z = zsize-1; 
    gc->face[5].center.x =  xsize/2; gc->face[5].center.y = ysize/2; gc->face[5].center.z =       0; 

    gc->alpha = -alpha;
    gc->beta = -beta;

    float diag = sqrtf(xsize*xsize+ysize*ysize+zsize*zsize);
    gc->diag = diag;
    iftVector d = {.x=-diag/2.0, .y=-diag/2.0, .z=-diag/2.0};
    iftVector c = {.x=xsize/2.0, .y=ysize/2.0, .z=zsize/2.0};
    
    gc->Tuv = iftTranslationMatrix(d);

    gc->Rx = iftRotationMatrix(IFT_AXIS_X, -alpha);
    gc->Ry = iftRotationMatrix(IFT_AXIS_Y, -beta);
    gc->Phi_r = iftMultMatrices(gc->Rx,gc->Ry);

    gc->Txyz = iftTranslationMatrix(c);

    iftMatrix *aux = iftMultMatrices(gc->Phi_r, gc->Tuv);

    gc->Phi = iftMultMatrices(gc->Txyz, aux);

    return gc;
}

int find_max_intensity(iftImage *img, iftPoint p1, iftPoint pn) {
  float n;
  float dx, dy, dz;
  int max_intensity;

  if(isPointsEqual(p1, pn)){
    n = 1;
  } else {
    float D_x = pn.x - p1.x;
    float D_y = pn.y - p1.y;
    float D_z = pn.z - p1.z;

    if (absolute(D_x) >= absolute(D_y) && absolute(D_x) >= absolute(D_z)) {
      n = absolute(D_x) + 1;
      dx = iftSign(D_x);
      dy = (dx*D_y)/D_x;
      dz = (dx*D_z)/D_x;
    } else if (absolute(D_y) >= absolute(D_x) && absolute(D_y) >= absolute(D_z)) {
      n = absolute(D_y) + 1;
      dy = iftSign(D_y);
      dx = (dy*D_x)/D_y;
      dz = (dy*D_z)/D_y;
    } else {
      n = absolute(D_z) + 1;
      dz = iftSign(D_z);
      dx = (dz*D_x)/D_z;
      dy = (dz*D_y)/D_z;
    }
  }

  iftPoint pp = p1;

  int val = iftImageValueAtPoint(img, pp);
  max_intensity = val;

  iftPoint d = {dx, dy, dz};
  /*if (max_intensity > 0)
    return max_intensity;*/
  for (int k = 1; k < n-1; k++) {
      pp = addPoints(pp, d);
      val = iftImageValueAtPoint(img, pp);
      max_intensity = val > max_intensity ? val : max_intensity; 
      /*if(max_intensity > 0)
        return max_intensity;*/
  }

  pp = pn;
  
  val = iftImageValueAtPoint(img, pp);
  max_intensity = val > max_intensity ? val : max_intensity; 
  
  return max_intensity;
}

iftImage *maximum_intesity_projection(iftImage *img, GraphicalContext *gc) {
  iftImage *projection = iftCreateImage(gc->diag, gc->diag, 1);

  float diag = gc->diag;

  iftPoint p = {.x=0, .y=0, .z=-diag/2};

  iftPoint n = {.x=0, .y=0, .z=1};

  iftPoint np = iftTransformPoint(gc->Phi_r, n); 

  for(p.x = 0; p.x < diag; p.x++){
    for(p.y = 0; p.y < diag; p.y++){
      iftPoint p0 = iftTransformPoint(gc->Phi, p);

      iftPoint p1;
      iftPoint pn;

      float min_lambda = IFT_INFINITY_FLT;
      float max_lambda = IFT_INFINITY_FLT_NEG;

      for(int f = 0; f < 6; f++) {
        Face face = gc->face[f];
        iftPoint normal = face.normal;
        iftPoint center = face.center;
        
        //dot product between p0 and face vector
        double dp_p0_fn = PointsDotProd(p0, normal);
        //dot product between vector face center vector 
        double dp_fc_fn = PointsDotProd(center, normal);
        //dot product between face vector and n'
        double dp_fn_np = PointsDotProd(normal, np);

        if(iftAlmostZero(dp_fn_np)){
          continue;
        }
        

        double lambda = (dp_fc_fn - dp_p0_fn)/dp_fn_np;
        double dp = dp_fn_np * lambda - dp_fc_fn + dp_p0_fn;


        if(iftAlmostZero(dp)){
          if(lambda < min_lambda){
            min_lambda = lambda;
          }

          if(lambda > max_lambda){
            max_lambda = lambda;
          }
        }
      }
      p1 = addPoints(p0, mulByScalar(np, min_lambda));
      pn = addPoints(p0, mulByScalar(np, max_lambda));

      iftVoxel p_ = {p.x, p.y, 0};
      int index = iftGetVoxelIndex(projection, p_);
      projection->val[index] = find_max_intensity(img, p1, pn);
    }
  }

  return projection;
}

int main(int argc, char *argv[]) 
{
  timer *tstart = NULL;
  int    MemDinInicial, MemDinFinal;
  
  MemDinInicial = iftMemoryUsed(1);

  if (argc != 6){
    iftError("Usage: getslice <...>\n"
	    "[1] input image .scn \n"
        "[2] the tilt angle alpha \n"
        "[3] is the spin angle beta \n"
        "[4] is an optional parameter — an .scn object mask \n"
        "[5] is the output .png image of the maximum intensity projection \n",
	     "main");
  }

  tstart = iftTic();

  /* ----------------------- Coding Area -------------------------- */

  iftImage *img    = iftReadImageByExt(argv[1]);

  float alpha = atof(argv[2]);
  float beta = atof(argv[3]);

  GraphicalContext *gc = create_graphical_context(img, alpha, beta);

  iftImage *projection = maximum_intesity_projection(img, gc);

  iftWriteImageByExt(projection, "teste.pgm");

  iftDestroyImage(&img);


  
  /* -------------------- End of the coding area ----------------- */
    
  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));
    
  MemDinFinal = iftMemoryUsed();
  iftVerifyMemory(MemDinInicial, MemDinFinal);
  
  return(0);
}