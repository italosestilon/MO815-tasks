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

iftPoint toPoint(float *a) {
  iftPoint p = {.x=a[0], .y=a[1], .z=a[2]};

  return p;
}

iftPoint addPoints(iftPoint a, iftPoint b) {
  iftPoint c = {a.x+b.x, a.y + b.y, a.z + b.z};

  return c;
}

iftPoint subPoints(iftPoint a, iftPoint b) {
  iftPoint c = {a.x -b.x, a.y - b.y, a.z - b.z};

  return c;
}

float pointNorm(iftPoint a){
  return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
}

float absolute(float a) {
  return a > 0? a: -1 * a;
}

float max(float a, float b) {
  return a > b? a : b;
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
  return fabs(a.x - b.x) <= IFT_EPSILON && fabs(a.y - b.y) <= IFT_EPSILON && fabs(a.z - b.z) <= IFT_EPSILON;
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

iftPoint find_surface_point(iftImage *img, iftPoint p1, iftPoint pn) {
  float n;
  float dx, dy, dz;

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

  iftPoint d = {dx, dy, dz};
  
  if (val > 0) {
    return pp;
  }

  for (int k = 1; k < n-1; k++) {
      pp = addPoints(pp, d);
      val = iftImageValueAtPoint(img, pp);
      if (val > 0) {
        return pp;
      }
  }

  pp.x = -1; pp.y = -1; pp.z = -1;

  return pp;
}

iftMImage *compute_gradient(iftImage *img, float adjacency_radius) {
  iftMImage *grad = iftCreateMImage(img->xsize, img->ysize, img->zsize, 3);
  
  iftAdjRel *A = iftSpheric(adjacency_radius);

  #pragma omp parallel for 
  for(int i = 0; i < img->n; i++) {
    iftVoxel u = iftGetVoxelCoord(img, i);
    iftPoint grad_vector;

    for( int j = 0; j < A->n; j++) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, j);
      int q = iftGetVoxelIndex(img, v);
      
      float intensity_diff = img->val[i] - img->val[q];
      iftPoint diff = subPoints(iftVoxelToPoint(v), iftVoxelToPoint(u));
      grad_vector = mulByScalar(diff, intensity_diff/pointNorm(diff));

    }

    grad->val[i][0] = grad_vector.x;
    grad->val[i][1] = grad_vector.y;
    grad->val[i][2] = grad_vector.z;

  }

  return grad;

}

iftPoint interpolate_grad(iftMImage *grad, iftVoxel pp) {
    iftVoxel u[8];
    int p[8], i;
    float dx,dy,dz;
    iftPoint val[6], value;
    
    u[0].x = (int)pp.x;      u[0].y = (int)pp.y;       u[0].z = (int)pp.z;
    u[1].x = u[0].x+1;      u[1].y = u[0].y;         u[1].z = u[0].z;
    u[2].x = u[0].x;        u[2].y = u[0].y + 1;     u[2].z = u[0].z;
    u[3].x = u[0].x+1;      u[3].y = u[0].y + 1;     u[3].z = u[0].z;
    u[4].x = u[0].x;        u[4].y = u[0].y;         u[4].z = u[0].z + 1;
    u[5].x = u[0].x+1;      u[5].y = u[0].y;         u[5].z = u[0].z + 1;
    u[6].x = u[0].x;        u[6].y = u[0].y + 1;     u[6].z = u[0].z + 1;
    u[7].x = u[0].x+1;      u[7].y = u[0].y + 1;     u[7].z = u[0].z + 1;

    for (i=0; i < 8; i++) {
        if (iftValidVoxel(grad,u[i])){
            p[i] = iftGetVoxelIndex(grad,u[i]);
        }else{
            p[0] = iftGetVoxelIndex(grad,u[0]);
            iftPoint grad_p = {
              .x=grad->val[p[0]][0], 
              .y=grad->val[p[0]][1], 
              .z=grad->val[p[0]][2]};
            
            return grad_p;
        }
    }
    
    val[0] = addPoints(mulByScalar(toPoint(grad->val[p[1]]), dx), mulByScalar(toPoint(grad->val[p[0]]),(1.0-dx)));
    val[1] = addPoints(mulByScalar(toPoint(grad->val[p[3]]), dx), mulByScalar(toPoint(grad->val[p[2]]),(1.0-dx)));
    val[2] = addPoints(mulByScalar(toPoint(grad->val[p[5]]), dx), mulByScalar(toPoint(grad->val[p[4]]),(1.0-dx)));
    val[3] = addPoints(mulByScalar(toPoint(grad->val[p[7]]), dx), mulByScalar(toPoint(grad->val[p[8]]),(1.0-dx)));
    val[4] = addPoints(mulByScalar(val[1],dy), mulByScalar(val[0],(1.0-dy)));
    val[5] = addPoints(mulByScalar(val[3],dy), mulByScalar(val[2],(1.0-dy)));
    value  = addPoints(mulByScalar(val[5],dz), mulByScalar(val[4],(1.0-dz)));

    return value;
}

float phones_illumination(iftImage *img, iftImage *label_image, iftPoint grad_pp, iftPoint np, iftPoint pp, float k_a, float k_d, float n_s, float r_a) {

    float alpha = 1.0;
    iftPoint normal_pp;
    iftPoint p_extended = addPoints(pp, mulByScalar(grad_pp, alpha));


    int i = iftGetVoxelIndex(label_image, iftPointToVoxel(pp));
    int j = iftGetVoxelIndex(label_image, iftPointToVoxel(p_extended));


    if(label_image->val[i] == label_image->val[j]) {
      normal_pp = mulByScalar(grad_pp, -1/pointNorm(grad_pp));
    } else {
      normal_pp = mulByScalar(grad_pp, 1/pointNorm(grad_pp));
    }

    float tetha = acosf(PointsDotProd(normal_pp, np)/(pointNorm(normal_pp)*pointNorm(pp)));

    float r = k_a*r_a;

}

iftImage *surface_rendering(iftImage *img, GraphicalContext *gc) {
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

      if(min_lambda < max_lambda) {
        p1 = addPoints(p0, mulByScalar(np, min_lambda));
        pn = addPoints(p0, mulByScalar(np, max_lambda));

        iftVoxel p_ = {iftRound(p.x), iftRound(p.y), 0};
        int index = iftGetVoxelIndex(projection, p_);
        iftPoint pp = find_surface_point(img, p1, pn);
  
      }
    }
  }

  return projection;
}

void change_intesity_interval(iftImage *img, int h) {
  int min = iftMinimumValue(img);
  int max = iftMaximumValue(img);

  iftVoxel u = {0, 0, 0};

  for(u.x = 0; u.x < img->xsize; u.x++) {
    for(u.y = 0; u.y < img->ysize; u.y++) {
      int p = iftGetVoxelIndex(img,u);
      int val = img->val[p];
      float new_val = h*(val - min)/((float)(max - min));
      img->val[p] = iftRound(new_val);
    }
  }

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

iftImage *window_level(iftImage *img, float window_p, float level_p, int h){
  int img_min = iftMinimumValue(img);
  int img_max = iftMaximumValue(img);

  int window = iftRound((img_max - img_min)*window_p);
  int level = iftRound(img_max * level_p);
  int l1 = iftRound(level - level/2);
  int l2 = iftRound(level/2 + window);

  iftImage *streached = iftCreateImage(img->xsize,img->ysize, img->zsize);

  iftVoxel u = {.x = 0, .y = 0, .z = 0};

  for(u.x = 0; u.x < img->xsize; u.x++) {
    for(u.y = 0; u.y < img->ysize; u.y++) {
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

  return streached;

}

int main(int argc, char *argv[]) 
{
  timer *tstart = NULL;
  int    MemDinInicial, MemDinFinal;
  
  MemDinInicial = iftMemoryUsed(1);

  if (argc != 6){
    iftError("Usage: Surface rendering <...>\n"
	    "[1] the name of the original scene (.scn).\n"
      "[2] the name of the label scene (.scn).\n"
      "[3] is the viewing tilt angle alhpa.\n"
      "[4] is the viewing spin angle beta. \n",
      "[5] is the resulting rendition. \n",
      "main");
  }

  tstart = iftTic();

  /* ----------------------- Coding Area -------------------------- */

  iftImage *img = iftReadImageByExt(argv[1]);
  iftImage *label_image = iftReadImageByExt(argv[2]);

  float alpha = atof(argv[3]);
  float beta = atof(argv[4]);
  int h = 256*256-1;

  GraphicalContext *gc = create_graphical_context(img, alpha, beta);

  iftImage *rendering = surface_rendering(img, gc);

  iftWriteImageByExt(rendering, argv[5]);

  iftDestroyImage(&img);
  
  /* -------------------- End of the coding area ----------------- */
    
  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));
    
  MemDinFinal = iftMemoryUsed();
  iftVerifyMemory(MemDinInicial, MemDinFinal);
  
  return(0);
}