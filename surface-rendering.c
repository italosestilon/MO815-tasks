// Made by Italos
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
  float min_dist;
  float max_dist;
  float k_a, k_s, k_d, n_s, r_a;
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

GraphicalContext *create_graphical_context(iftImage *img, float alpha, float beta, int H) {
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

    gc->k_a = 0.1;
    gc->k_s = 0.7;
    gc->k_d = 0.2;
    gc->n_s = 5;
    gc->r_a = H;

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

iftImage *GetSliceSagital(iftImage *img, iftImage *label_image, int x, GraphicalContext *gc) {
  iftImage *slc = iftCreateColorImage(img->ysize,img->zsize,1, 8);
  iftColorTable *colorTabel = iftCategoricalColorTable(iftMaximumValue(label_image) + 1);
  
  iftVoxel  u;
  int q=0;
  
  u.x = x;
  
  for (u.z = img->zsize-1; u.z >= 0; u.z--) {
    for (u.y = 0; u.y < img->ysize; u.y++){
      int p = iftGetVoxelIndex(img,u);      

      int label = label_image->val[p];
      iftColor rgb_color = colorTabel->color[label];
  
      if(label == 0) {
        rgb_color.val[0] = img->val[p];
        rgb_color.val[1] = img->val[p];
        rgb_color.val[2] = img->val[p];
      }
      
      iftColor YCbCr_color = iftRGBtoYCbCr(rgb_color, gc->r_a);
      slc->val[q] = img->val[p];
      
      slc->Cb[q] = YCbCr_color.val[1];
      slc->Cr[q] = YCbCr_color.val[2];
      q++;
    }
  }
  return(slc);
}

iftImage *GetSliceCoronal(iftImage *img, iftImage *label_image, int y, GraphicalContext *gc) {
  iftImage *slc = iftCreateColorImage(img->xsize,img->zsize,1, 8);
  iftColorTable *colorTabel = iftCategoricalColorTable(iftMaximumValue(label_image) + 1);
  iftVoxel  u;
  int q=0;
  
  u.y = y;

  for (u.z = img->zsize-1; u.z >= 0; u.z--){
    for (u.x = 0; u.x < img->xsize; u.x++){
      int p = iftGetVoxelIndex(img,u);      
      slc->val[q] = img->val[p];
      int label = label_image->val[p];
      iftColor rgb_color = colorTabel->color[label];
  
      if(label == 0) {
        rgb_color.val[0] = img->val[p];
        rgb_color.val[1] = img->val[p];
        rgb_color.val[2] = img->val[p];
      }
      
      iftColor YCbCr_color = iftRGBtoYCbCr(rgb_color, gc->r_a);
      slc->val[q] = img->val[p];
      
      slc->Cb[q] = YCbCr_color.val[1];
      slc->Cr[q] = YCbCr_color.val[2];
      q++;
    }
  }
    

  return(slc);
}

iftImage *GetSliceAxial(iftImage *img, iftImage *label_image, int z, GraphicalContext *gc)
{
  iftImage *slc = iftCreateColorImage(img->xsize,img->ysize,1, 8);
  iftColorTable *colorTabel = iftCategoricalColorTable(iftMaximumValue(label_image) + 1);
  iftVoxel  u;
  int       q=0;
  
  u.z = z;

  for (u.y = 0; u.y < img->ysize; u.y++) {
    for (u.x = 0; u.x < img->xsize; u.x++){
      int p = iftGetVoxelIndex(img,u);      
      slc->val[q] = img->val[p];
      int label = label_image->val[p];
      iftColor rgb_color = colorTabel->color[label];
  
      if(label == 0) {
        rgb_color.val[0] = img->val[p];
        rgb_color.val[1] = img->val[p];
        rgb_color.val[2] = img->val[p];
      }
      
      iftColor YCbCr_color = iftRGBtoYCbCr(rgb_color, gc->r_a);
      slc->val[q] = img->val[p];
      
      slc->Cb[q] = YCbCr_color.val[1];
      slc->Cr[q] = YCbCr_color.val[2];

      q++;
    }
  }

  return(slc);
}

iftPoint find_surface_point(iftImage *label_image, iftPoint p1, iftPoint pn) {
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

  iftPoint d = {dx, dy, dz};
  int index;
  iftVoxel rounded = {iftRound(pp.x), iftRound(pp.y), iftRound(pp.z)};
  if(iftValidVoxel(label_image, rounded)) {
    index = iftGetVoxelIndex(label_image, rounded);

    if (label_image->val[index] != 0) {
      return pp;
    }
  }

  for (int k = 1; k < n-1; k++) {
      pp = addPoints(pp, d);
      rounded.x = iftRound(pp.x);
      rounded.y = iftRound(pp.y);
      rounded.z = iftRound(pp.z);
      if(iftValidVoxel(label_image, rounded)) {
        index = iftGetVoxelIndex(label_image, rounded);
        if (label_image->val[index] != 0) {
          return pp;
        }
      }
  }

  pp.x = -1; pp.y = -1; pp.z = -1;

  return pp;
}

iftMImage *compute_gradient(iftImage *img, iftImage *label_image, float adjacency_radius) {
  iftMImage *grads = iftCreateMImage(img->xsize, img->ysize, img->zsize, 3);
  
  iftAdjRel *A = iftSpheric(adjacency_radius);
  //iftImage *border = iftBorderImage(label_image, 1);
  float alpha = 2.0;

  #pragma omp parallel for 
  for(int i = 0; i < img->n; i++) {
    iftVoxel u = iftGetVoxelCoord(img, i);

    //if(border->val[i] == 0) continue;

    iftPoint grad_vector = {.x=0, .y=0, .z=0};

    for( int j = 1; j < A->n; j++) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, j);
      if(iftValidVoxel(img, v)) {
        int q = iftGetVoxelIndex(img, v);
        
        float intensity_diff = img->val[i] - img->val[q];
        iftPoint diff = subPoints(iftVoxelToPoint(v), iftVoxelToPoint(u));
        grad_vector = addPoints(grad_vector, mulByScalar(diff, intensity_diff/pointNorm(diff)));
      }

    }

    iftPoint p_near = addPoints(iftVoxelToPoint(u), mulByScalar(grad_vector, alpha));
    
    if(iftValidVoxel(label_image, iftPointToVoxel(p_near))) {
        int j = iftGetVoxelIndex(label_image, iftPointToVoxel(p_near));
        if((label_image->val[i] == label_image->val[j])) {
            grad_vector = mulByScalar(grad_vector, -1/pointNorm(grad_vector));
        } else {
            grad_vector = mulByScalar(grad_vector, 1/pointNorm(grad_vector));
        }
    } else {
        grad_vector = mulByScalar(grad_vector, 1/pointNorm(grad_vector));
    }

    grads->val[i][0] = grad_vector.x;
    grads->val[i][1] = grad_vector.y;
    grads->val[i][2] = grad_vector.z;

  }

  return grads;

}

iftPoint interpolate_grad(iftMImage *grad, iftVoxel pp) {
    iftVoxel u[8];
    int p[8], i;
    float dx,dy,dz;
    iftPoint val[6], value;
    
    u[0].x = (int)pp.x;     u[0].y = (int)pp.y;      u[0].z = (int)pp.z;
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
            pp.x = iftRound(pp.x);
            pp.y = iftRound(pp.y);
            pp.z = iftRound(pp.z);
            int index = iftGetVoxelIndex(grad, pp);
            iftPoint grad_p;
            grad_p.x = grad->val[index][0];
            grad_p.y = grad->val[index][1];
            grad_p.z = grad->val[index][2];
            return grad_p;
        }
    }
    dx = 1.0, dy = 1.0, dz = 1.0;
    val[0] = addPoints(mulByScalar(toPoint(grad->val[p[1]]), dx), mulByScalar(toPoint(grad->val[p[0]]),(1.0-dx)));
    val[1] = addPoints(mulByScalar(toPoint(grad->val[p[3]]), dx), mulByScalar(toPoint(grad->val[p[2]]),(1.0-dx)));
    val[2] = addPoints(mulByScalar(toPoint(grad->val[p[5]]), dx), mulByScalar(toPoint(grad->val[p[4]]),(1.0-dx)));
    val[3] = addPoints(mulByScalar(toPoint(grad->val[p[7]]), dx), mulByScalar(toPoint(grad->val[p[6]]),(1.0-dx)));
    val[4] = addPoints(mulByScalar(val[1],dy), mulByScalar(val[0],(1.0-dy)));
    val[5] = addPoints(mulByScalar(val[3],dy), mulByScalar(val[2],(1.0-dy)));
    value  = addPoints(mulByScalar(val[5],dz), mulByScalar(val[4],(1.0-dz)));
    
    return value;
}

double phones_model(iftImage *img, iftPoint grad_pp, iftPoint np, iftPoint pp, 
                    iftPoint p0, GraphicalContext *gc) {
    double r = 0;
    if (grad_pp.x == 0 && grad_pp.y == 0 && grad_pp.z == 0) {
      printf("Zero gradient...\n");
    }

    double magnitudes = (pointNorm(grad_pp)*pointNorm(np));
    
    double tetha = acos(PointsDotProd(grad_pp, mulByScalar(np, -1))/magnitudes);

    if(tetha >=  0 && tetha < PI/2.0) {
      double r_d = gc->r_a*(pointNorm(subPoints(pp, p0)) - gc->min_dist)/(gc->max_dist - gc->min_dist);
      float specular = tetha >= 0 && tetha < PI/4.0 ? gc->k_s * pow(cos(2*tetha), gc->n_s) : 0.0;
      r = gc->k_a*gc->r_a + r_d * (gc->k_d * cos(tetha) + specular);
    }
    if(r < 1) {
      printf("Oopds");
    }
    return r;
}

void distance_map(iftImage *img, iftPoint p, GraphicalContext *gc) {
  iftVoxel u = {.x=0, .y=0, .z=0};
  float max_dist = 1;
  float min_dist = FLT_MAX;
  for(u.x = 0; u.x < img->xsize; u.x++) {
    for(u.y = 0; u.y < img->ysize; u.y++) {
      for(u.z = 0; u.z < img->zsize; u.z++) {
        float dist = pointNorm(subPoints(iftVoxelToPoint(u), p));

        if (dist > max_dist) {
          max_dist = dist;
        }

        if (dist > 0 && dist < min_dist) {
          min_dist = dist;
        }
      }
    }
  }

  gc->max_dist = max_dist;
  gc->min_dist = min_dist;
}

iftImage *surface_rendering(iftImage *img, iftImage *label_image, GraphicalContext *gc) {
  iftImage *projection =  iftCreateColorImage(gc->diag, gc->diag, 1, 8);

  iftMImage *grads = compute_gradient(img, label_image, 3.0);
  
  float diag = gc->diag;

  iftPoint p = {.x=0, .y=0, .z=-diag/2};

  iftPoint c= {.x=diag/2, .y=diag/2, .z=diag/2};
  iftPoint c_p = iftTransformPoint(gc->Phi, c);

  distance_map(img, c_p, gc);

  iftPoint n = {.x=0, .y=0, .z=1};

  iftPoint np = iftTransformPoint(gc->Phi_r, n);

  iftColorTable *colorTabel = iftCategoricalColorTable(iftMaximumValue(label_image) + 1);

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

        iftPoint pp = find_surface_point(label_image, p1, pn);

        if (pp.x >= 0 && pp.y >= 0 && pp.z >= 0) {
          iftPoint grad_point = interpolate_grad(grads, iftPointToVoxel(pp));
          float r = phones_model(img, grad_point, np, pp, c_p, gc);

          iftVoxel p_label = {iftRound(pp.x), iftRound(pp.y), iftRound(pp.z)};
          int index_label = iftGetVoxelIndex(label_image, p_label);

          int label = label_image->val[index_label];
          iftColor rgb_color = colorTabel->color[label];

          iftColor YCbCr_color = iftRGBtoYCbCr(rgb_color, gc->r_a);
      
          projection->val[index] = iftRound(r);
          projection->Cb[index] = YCbCr_color.val[1];
          projection->Cr[index] = YCbCr_color.val[2];

        }
      }
    }
  }
  iftDestroyMImage(&grads);
  return projection;
}

void change_intesity_interval(iftImage *img, int h) {
  int min = iftMinimumValue(img);
  int max_v = iftMaximumValue(img);

  iftVoxel u = {0, 0, 0};

  for(u.x = 0; u.x < img->xsize; u.x++) {
    for(u.y = 0; u.y < img->ysize; u.y++) {
      int p = iftGetVoxelIndex(img,u);
      int val = img->val[p];
      float new_val = h*(val - min)/((float)(max_v - min));
      img->val[p] = iftRound(new_val);
    }
  }

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
  int h = 256-1;

  GraphicalContext *gc = create_graphical_context(img, alpha, beta, h);

  iftImage *rendering = surface_rendering(img, label_image, gc);
  iftImage *slc_sagital = GetSliceSagital(img, label_image, 120, gc);
  iftImage *slc_coronal = GetSliceCoronal(img, label_image, 120, gc);
  iftImage *slc_axial = GetSliceAxial(img, label_image, 120, gc);

  change_intesity_interval(slc_sagital, h);
  change_intesity_interval(slc_coronal, h);
  change_intesity_interval(slc_axial, h);

  iftWriteImageByExt(rendering, argv[5]);
  iftWriteImageByExt(slc_sagital, "sagital.png");
  iftWriteImageByExt(slc_coronal, "coronal.png");
  iftWriteImageByExt(slc_axial, "axial.png");  

  iftDestroyImage(&img);
  iftDestroyImage(&rendering);
  iftDestroyImage(&slc_sagital);
  iftDestroyImage(&slc_coronal);
  iftDestroyImage(&slc_axial);
   
  /* -------------------- End of the coding area ----------------- */
    
  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));
    
  MemDinFinal = iftMemoryUsed();
  iftVerifyMemory(MemDinInicial, MemDinFinal);
  
  return(0);
}