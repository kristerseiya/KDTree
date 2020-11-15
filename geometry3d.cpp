#include "kdtree.hpp"
#include "geometry3d.hpp"
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

using namespace std;

Geometry3D::Geometry3D() {
  this->path = new char[100];
  this->path[0] = '\0';
  this->points = NULL;
  this->normals = NULL;
  this->colors = NULL;
  this->n_points = 0;
  this->is_kdtree_built = false;
  this->is_empty = true;
}

Geometry3D& Geometry3D::Clear() {
  this->path[0] = '\0';
  if (this->points != NULL) { delete[] this->points; }
  if (this->colors != NULL) { delete[] this->colors; }
  if (this->normals != NULL) { delete[] this->normals; }
  this->is_empty = true;
  return *this;
}

void Geometry3D::read_stl(char* path) {
  if (!this->is_empty) {
    this->Clear();
  }
  FILE* fp = fopen(path,"rb");
  if (fp == NULL) {
    printf("failed to open file\n");
    exit(1);
  }
  char header[80];
  fread(header,sizeof(char),80,fp);
  unsigned int n_face;
  fread(&n_face, sizeof(unsigned int), 1, fp);
  this->n_points = ((size_t)n_face) * 3;

  float* normals = new float[((size_t)n_face) * 3];
  float* vertices = new float[((size_t)n_face) * 9];
  short attr;

  for (size_t i = 0; i < n_face; i++) {
    fread(normals+i*3, sizeof(float), 3, fp);
    fread(vertices+i*9, sizeof(float), 9, fp);
    fread(&attr, sizeof(short), 1, fp);
  }

  fclose(fp);

  this->points = vertices;
  this->normals = normals;
  this->is_empty = false;
  printf("successfully read stl file\n");
}

void Geometry3D::read_ply(char* path) {
  if (!this->is_empty) {
    this->Clear();
  }
  char line[101];
  float pt[3];
  float n_pt[3];
  uint8 color[3];
  FILE* ply_file = fopen(path,"rb");
  if (ply_file != NULL) {
    for (int i = 0; i < 13; i++) {
      fgets(line,100,ply_file);
    }
    size_t n_points = 157006;
    this->n_points = n_points;
    this->points = new float[3*n_points];
    this->normals = new float[3*n_points];
    this->colors = new uint8[3*n_points];
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
  this->n_points = 157006;
  this->is_empty = false;
}

void Geometry3D::read_xyzm(char* path) {
  if (!this->is_empty) {
    this->Clear();
  }
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
  this->n_points = width * height;
  while(fgetc(xyzm_file)==0);
  fseek(xyzm_file,-1,SEEK_CUR);
  float* xyz = new float[3*width*height];
  uint8* rgb = new uint8[3*width*height];
  // uint8* mask = new uint8[width*height];
  fread(xyz,sizeof(float),3*width*height,xyzm_file);
  fread(rgb,sizeof(uint8),3*width*height,xyzm_file);
  // fread(mask,sizeof(uint8),width*height,xyzm_file);
  fclose(xyzm_file);
  this->points = xyz;
  this->colors = rgb;
  this->is_empty = false;
  // this->mask = mask;
  printf("successfully read xyzm file\n");
}

Geometry3D::Geometry3D(char* path) {
  this->points = NULL;
  this->normals = NULL;
  this->colors = NULL;
  this->is_empty = true;
  this->path = new char[100];
  memcpy(this->path,path,100);
  int path_len = strlen(path);
  if (!strcmp(path+path_len-4,".ply")) {
    printf(".ply detected\n");
    this->read_ply(path);
  } else if (!strcmp(path+path_len-5,".xyzm")) {
    printf(".xyzm detected\n");
    this->read_xyzm(path);
  } else if (!strcmp(path+path_len-4,".stl")) {
    printf(".stl detected\n");
    this->read_stl(path);
  } else {
    fprintf(stderr,"couldn't read file\n");
    exit(1);
  }
  this->is_kdtree_built = false;
}

Geometry3D::~Geometry3D() {
  if (this->path != NULL) { delete[] this->path; }
  if (this->points != NULL) { delete[] this->points; }
  if (this->colors != NULL) { delete[] this->colors; }
  if (this->normals != NULL) { delete[] this->normals; }
}

void Geometry3D::build_KDTree() {
    this->kdtree.set(this->points,3,this->n_points);
    this->is_kdtree_built = true;
    return;
}

int Geometry3D::find_k_nearest(float* query, int k, size_t* k_nearest, float* k_distances) {
  int a =  this->kdtree.find_k_nearest(query,k,k_nearest,k_distances);
  return a;
}
