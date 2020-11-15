#include "kdtree.hpp"
#include "pc_class.hpp"
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

using namespace std;

PointCloud::PointCloud() {
  char x[] = "None";
  memcpy(this->path,x,5);
  this->points = NULL;
  this->normals = NULL;
  this->colors = NULL;
  this->mask = NULL;
  this->num_points = 0;
  this->is_kdtree_built = false;
}

void PointCloud::read_ply(char* path) {
  char line[101];
  float pt[3];
  float n_pt[3];
  uint8 color[3];
  FILE* ply_file = fopen(path,"rb");
  if (ply_file != NULL) {
    for (int i = 0; i < 13; i++) {
      fgets(line,100,ply_file);
      // cout << line;
    }
    size_t num_points = 157006;
    this->num_points = num_points;
    this->points = new float[3*num_points];
    this->normals = new float[3*num_points];
    this->colors = new uint8[3*num_points];
    float* coords = this->points;
    float* normals = this->normals;
    uint8* colors = this->colors;
    while (!feof(ply_file)) {
      fread(pt,sizeof(float),3,ply_file);
      fread(n_pt,sizeof(float),3,ply_file);
      fread(color,sizeof(uint8),3,ply_file);

      if (!feof(ply_file)) {
        uint8 bytes[4];
        for (int i = 0; i < 3; i++) {
          memcpy(bytes,&pt[i],4);
          int big_endian = bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24);
          memcpy(&pt[i],&big_endian,4);
          memcpy(bytes,&(n_pt[i]),4);
          big_endian = bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24);
          memcpy(&n_pt[i],&big_endian,4);
        }
        memcpy(coords,pt,3*sizeof(float));
        memcpy(normals,n_pt,3*sizeof(float));
        memcpy(colors,color,3*sizeof(uint8));
        coords = coords + 3;
        normals = normals + 3;
        colors = colors + 3;
      }
    }
    fclose(ply_file);
  } else {
    printf("Failed to read file\n");
  }
  this->num_points = 157006;
}

void PointCloud::read_xyzm(char* path) {
  FILE* xyzm_file = fopen(path,"rb");
  if (xyzm_file == NULL) {
    printf("failed to open file\n");
    exit(1);
  }
  int height;
  int width;
  int ret = fscanf(xyzm_file, "image size width x height = %d x %d",&width,&height);
  if (ret!=2) {
    printf("unsupported file format\n");
    fclose(xyzm_file);
    exit(1);
  }
  this->width= width;
  this->height = height;
  this->num_points = width * height;
  while(fgetc(xyzm_file)==0);
  fseek(xyzm_file,-1,SEEK_CUR);
  float* xyz = new float[3*width*height];
  uint8* rgb = new uint8[3*width*height];
  uint8* mask = new uint8[width*height];
  fread(xyz,sizeof(float),3*width*height,xyzm_file);
  fread(rgb,sizeof(uint8),3*width*height,xyzm_file);
  fread(mask,sizeof(uint8),width*height,xyzm_file);
  fclose(xyzm_file);
  this->points = xyz;
  this->colors = rgb;
  this->mask = mask;
  printf("successfully read xyzm file\n");
}

PointCloud::PointCloud(char* path) {
  this->points = NULL;
  this->normals = NULL;
  this->colors = NULL;
  this->mask = NULL;
  this->path = new char[100];
  memcpy(this->path,path,100);
  int path_len = strlen(path);
  if (!strcmp(path+path_len-4,".ply")) {
    printf(".ply detected\n");
    this->read_ply(path);
  } else if (!strcmp(path+path_len-5,".xyzm")) {
    printf(".xyzm detected\n");
    this->read_xyzm(path);
  } else {
    fprintf(stderr,"couldn't read file\n");
    exit(1);
  }
  this->is_kdtree_built = false;
}

PointCloud::~PointCloud() {
  if (this->path != NULL) { delete[] this->path; }
  if (this->points != NULL) { delete[] this->points; }
  if (this->colors != NULL) { delete[] this->colors; }
  if (this->normals != NULL) { delete[] this->normals; }
  if (this->mask != NULL) { delete[] this->mask; }
}

void PointCloud::build_KDTree() {
    this->kdtree.set(this->points,3,this->num_points);
    this->is_kdtree_built = true;
    return;
}

int PointCloud::find_k_nearest(float* query, int k, size_t* k_nearest, float* k_distances) {
  int a =  this->kdtree.find_k_nearest(query,k,k_nearest,k_distances);
  return a;
}
