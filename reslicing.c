#include "ift.h"

typedef struct _graphcontext{
  iftMatrix *Rx, *Ry, *Tuv, *Psi, *Psi_r;
  float diag;
} GraphicalContext;

iftPoint addPoints(iftPoint a, iftPoint b) {
  iftPoint c = {a.x+b.x, a.y + b.y, a.z + b.z};

  return c;
}

float absolute(float a) {
  return a > 0? a: -1 * a;
}

float max(float a, float b) {
  return a > b? a : b;
}

float pointNorm(iftPoint p) {
  return sqrtf(p.x*p.x + p.y * p.y + p.z * p.z);
}

iftPoint mulByScalar(iftPoint a, float s) {
  iftPoint c = {a.x * s, a.y * s, a.z * s};
  return c;
}

int validPoint(iftImage *img, iftPoint p){
  return ceil(p.x) >= 0 && ceil(p.x) < img->xsize && ceil(p.y) >= 0 && ceil(p.y) < img->ysize && ceil(p.z) >= 0 && ceil(p.z) < img->zsize;
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

    float diag = sqrtf(xsize*xsize+ysize*ysize+zsize*zsize);
    gc->diag = diag;
    iftVector d = {.x=-diag/2.0, .y=-diag/2.0, .z=diag/2.0};

    gc->Tuv = iftTranslationMatrix(d);

    gc->Rx = iftRotationMatrix(IFT_AXIS_X, -alpha);
    gc->Ry = iftRotationMatrix(IFT_AXIS_Y, beta);
    gc->Psi_r = iftMultMatrices(gc->Rx,gc->Ry);

    //iftPrintMatrix(gc->Psi_r);

    gc->Psi = iftMultMatrices(gc->Psi_r, gc->Tuv);


    return gc;
}

void get_slice(GraphicalContext *gc,  iftImage *original_img, iftImage *resliced_img, iftPoint p_0, int slice){
    
    iftPoint p = {0, 0, -gc->diag/2};
    
    //p_0 = mulByScalar(p_0, -1);
    iftVector v_p0 = {.x=p_0.x, .y=p_0.y, .z=p_0.z};
    iftMatrix *Txyz = iftTranslationMatrix(v_p0);

    iftMatrix *transform = iftMultMatrices(Txyz, gc->Psi);

    for(p.x = 0; p.x < resliced_img->xsize; p.x++){
      for(p.y = 0; p.y < resliced_img->ysize; p.y++){
        iftPoint pp = iftTransformPoint(transform, p);

        if(validPoint(original_img, pp)){
            iftVoxel p_ = {p.x, p.y, slice};
            int i = iftGetVoxelIndex(resliced_img, p_);
            int val = iftImageValueAtPoint(original_img, pp);
            resliced_img->val[i] = val;
        }
      }
    }
}

void find_alpha_and_betha(iftPoint p_np, float *alpha, float *betha) {

  float x_magnitude = (sqrtf(p_np.x * p_np.x + p_np.z*p_np.z));
  float cos_alpha = iftAlmostZero(x_magnitude) ? 1 : (p_np.z)/x_magnitude;
  float _alpha = acosf(cos_alpha);

  float z_ = p_np.x *-sinf(-_alpha) + p_np.z * cosf(-_alpha);
  float y_magnitude = (sqrtf(p_np.y * p_np.y + z_*z_));
  float cos_betha = iftAlmostZero(y_magnitude) ? 1 : (z_)/y_magnitude;
  float _betha = acosf(cos_betha);

  //convert alpha and betha to degrees
  _alpha = _alpha*180/PI;
  _betha = _betha*180/PI;

  *alpha = _alpha;
  *betha = _alpha;
}

iftImage *reslice_image(GraphicalContext *gc, iftImage *img, int num_slices, iftPoint p0, iftPoint np, float lambda){
  iftImage *resliced_image = iftCreateImage(gc->diag, gc->diag, num_slices);

  iftPoint pk = p0;

  get_slice(gc, img, resliced_image, pk, 0);

  for(int i=1; i < num_slices; i++){
    pk = addPoints(p0, mulByScalar(np, lambda*i));
    get_slice(gc, img, resliced_image, pk, i);
  }

  return resliced_image;
}

int main(int argc, char *argv[]) 
{
  timer *tstart = NULL;
  int    MemDinInicial, MemDinFinal;
  
  MemDinInicial = iftMemoryUsed(1);

  if (argc != 10){
    iftError("Usage: getslice <...>\n"
	    "[1] input image .scn \n"
        "[2] P_0 ex.: 139 139 0 \n"
        "[3] P_n-1 ex: 139 270 200\n"
        "[4] The number of n of axial slices of the new scene \n"
        "[5] the output .scn scene \n",
	     "main");
  }

  tstart = iftTic();

  /* ----------------------- Coding Area -------------------------- */

  iftImage *img = iftReadImageByExt(argv[1]);

  iftPoint p_0 = {atoi(argv[2]), atoi(argv[3]), atoi(argv[4])};
  iftPoint p_n = {atoi(argv[5]), atoi(argv[6]), atoi(argv[7])};

  iftPoint p_np = addPoints(p_n, mulByScalar(p_0, -1));
  float norm = pointNorm(p_np);
  p_np.x = p_np.x/norm;
  p_np.y = p_np.y/norm;
  p_np.z = p_np.z/norm;

  int n = atoi(argv[8]);
  float lambda = norm/n;

  float alpha, betha;

  find_alpha_and_betha(p_np, &alpha, &betha);

  GraphicalContext *gc = create_graphical_context(img, alpha, betha);

  iftImage *resliced_image = reslice_image(gc, img, n, p_0, p_np, lambda);

  iftWriteImageByExt(resliced_image, argv[9]);

  iftDestroyImage(&resliced_image);
  iftDestroyImage(&img);
  /* -------------------- End of the coding area ----------------- */
    
  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));
    
  MemDinFinal = iftMemoryUsed();
  iftVerifyMemory(MemDinInicial, MemDinFinal);
  
  return(0);
}