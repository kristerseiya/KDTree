//
// Created by Krister Ulvog on 1/1/21.
//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
//#include <Eigen/Core>
#include <Eigen/StdVector>

static void swap_bytes(char* x, int size) {
    char* tmp = (char*)malloc(sizeof(char)*size);
    for (int i = 0; i < size; i++) {
        tmp[i] = x[size-1-i];
    }
    for (int i = 0; i < size; i++) {
        x[i] = tmp[i];
    }
    free(tmp);
}

static bool is_little_endian() {
    int x = 1;
    char* y = (char*)malloc(sizeof(int));
    memcpy(y,&x,sizeof(int));
    char z = y[0];
    free(y);
    return z;
}

void ReadPLY(const std::string& filename, std::vector<Eigen::Vector3d>& points,
             std::vector<Eigen::Vector3d>& normals,
             std::vector<Eigen::Vector3d>& colors) {

    points.clear();
    normals.clear();
    colors.clear();

    FILE* fp = fopen(filename.c_str(), "rb");
    if (fp == nullptr) {
        // fprintf(stderr, "could not open file\n");
        // exit(1);
        throw std::runtime_error("could not open file");
    }

    char buffer[100];
    char format[100];
    int size;
    char property_type[100];
    char property_name[100];
    int n_property = 0;
    int vertex_byte_offset[3];
    int byte_size = 0;
    bool is_double = false;
    bool normal_exists = false;
    int normal_byte_offset[3];
    bool color_exists = false;
    int color_byte_offset[3];

    fgets(buffer, 99, fp);
    // char str_ply[] = "ply\n";
    if (strcmp(buffer,"ply\n")!=0) {
        fprintf(stderr,"corrupted .ply file\n");
    }

    fgets(buffer, 99, fp);
    sscanf(buffer, "format %s", format);

    fgets(buffer, 99, fp);
    // char str_element_vertex[] = "element vertex";
    while (strncmp(buffer, "element vertex", 14) != 0) {
        fgets(buffer, 99, fp);
    }

    sscanf(buffer, "element vertex %d", &size);

    // char str_float[] = "float";
    // char str_double[] = "double";
    // char str_uchar[] = "uchar";
    fgets(buffer, 99, fp);
    while (sscanf(buffer, "property %s %s",
                  property_type, property_name)==2) {
        n_property++;
        if (strcmp(property_name, "x")==0) {
            vertex_byte_offset[0] = byte_size;
            if (strcmp(property_type, "double")==0) {
                is_double = true;
            }
        }
        if (strcmp(property_name, "y")==0) {
            vertex_byte_offset[1] = byte_size;
        }
        if (strcmp(property_name, "z")==0) {
            vertex_byte_offset[2] = byte_size;
        }
        if (strcmp(property_name, "nx")==0) {
            normal_exists = true;
            normal_byte_offset[0] = byte_size;
        }
        if (strcmp(property_name, "ny")==0) {
            normal_byte_offset[1] = byte_size;
        }
        if (strcmp(property_name, "nz")==0) {
            normal_byte_offset[2] = byte_size;
        }
        if (strcmp(property_name, "red")==0) {
            color_exists = true;
            color_byte_offset[0] = byte_size;
        }
        if (strcmp(property_name, "green")==0) {
            color_byte_offset[1] = byte_size;
        }
        if (strcmp(property_name, "blue")==0) {
            color_byte_offset[2] = byte_size;
        }

        if (strcmp(property_type,"float")==0) {
            byte_size += sizeof(float);
        } else if (strcmp(property_type,"double")==0) {
            byte_size += sizeof(double);
        } else if (strcmp(property_type,"uchar")==0) {
            byte_size += sizeof(unsigned char);
        } else {
            // fprintf(stderr,"unknown property type\n");
            // exit(1);
            throw std::runtime_error("unknown property type");
        }
        fgets(buffer, 99, fp);
    }

    // char str_end_header[] = "end_header\n";
    while(strcmp(buffer, "end_header\n")!=0) {
        fgets(buffer, 99, fp);
    }

    // char str_ascii[] = "ascii";
    // char str_binary[] = "binary_little_endian";
    if (strcmp(format,"binary_little_endian")==0) {
        points.resize(size);
        if (normal_exists) {
            normals.resize(size);
        }
        if (color_exists) {
            colors.resize(size);
        }
        if (is_double) {
            double x, y, z;
            double nx, ny, nz;
            unsigned char r, g, b;
            bool swap_endian = !is_little_endian();
            for (int i = 0; i < size; i++) {
                fread(buffer, 1, byte_size, fp);
                if (swap_endian) {
                    swap_bytes(buffer+vertex_byte_offset[0],sizeof(double));
                    swap_bytes(buffer+vertex_byte_offset[1],sizeof(double));
                    swap_bytes(buffer+vertex_byte_offset[2],sizeof(double));
                }
                memcpy(&x,buffer+vertex_byte_offset[0],sizeof(double));
                memcpy(&y,buffer+vertex_byte_offset[1],sizeof(double));
                memcpy(&z,buffer+vertex_byte_offset[2],sizeof(double));
                points[i][0] = x;
                points[i][1] = y;
                points[i][2] = z;
                if (normal_exists) {
                    if (swap_endian) {
                        swap_bytes(buffer+normal_byte_offset[0],sizeof(double));
                        swap_bytes(buffer+normal_byte_offset[1],sizeof(double));
                        swap_bytes(buffer+normal_byte_offset[2],sizeof(double));
                    }
                    memcpy(&nx,buffer+normal_byte_offset[0],sizeof(double));
                    memcpy(&ny,buffer+normal_byte_offset[1],sizeof(double));
                    memcpy(&nz,buffer+normal_byte_offset[2],sizeof(double));
                    normals[i][0] = nx;
                    normals[i][1] = ny;
                    normals[i][2] = nz;
                }
                if (color_exists) {
                    memcpy(&r,buffer+color_byte_offset[0],sizeof(unsigned char));
                    memcpy(&g,buffer+color_byte_offset[1],sizeof(unsigned char));
                    memcpy(&b,buffer+color_byte_offset[2],sizeof(unsigned char));
                    colors[i][0] = r / 255.;
                    colors[i][1] = g / 255.;
                    colors[i][2] = b / 255.;
                }
            }
        } else {
            float x, y, z;
            float nx, ny, nz;
            unsigned char r, g, b;
            bool swap_endian = !is_little_endian();
            for (int i = 0; i < size; i++) {
                fread(buffer, 1, byte_size, fp);
                if (swap_endian) {
                    swap_bytes(buffer+vertex_byte_offset[0],sizeof(float));
                    swap_bytes(buffer+vertex_byte_offset[1],sizeof(float));
                    swap_bytes(buffer+vertex_byte_offset[2],sizeof(float));
                }
                memcpy(&x,buffer+vertex_byte_offset[0],sizeof(float));
                memcpy(&y,buffer+vertex_byte_offset[1],sizeof(float));
                memcpy(&z,buffer+vertex_byte_offset[2],sizeof(float));
                points[i][0] = x;
                points[i][1] = y;
                points[i][2] = z;
                if (normal_exists) {
                    if (swap_endian) {
                        swap_bytes(buffer+normal_byte_offset[0],sizeof(float));
                        swap_bytes(buffer+normal_byte_offset[1],sizeof(float));
                        swap_bytes(buffer+normal_byte_offset[2],sizeof(float));
                    }
                    memcpy(&nx,buffer+normal_byte_offset[0],sizeof(float));
                    memcpy(&ny,buffer+normal_byte_offset[1],sizeof(float));
                    memcpy(&nz,buffer+normal_byte_offset[2],sizeof(float));
                    normals[i][0] = nx;
                    normals[i][1] = ny;
                    normals[i][2] = nz;
                }
                if (color_exists) {
                    memcpy(&r,buffer+color_byte_offset[0],sizeof(unsigned char));
                    memcpy(&g,buffer+color_byte_offset[1],sizeof(unsigned char));
                    memcpy(&b,buffer+color_byte_offset[2],sizeof(unsigned char));
                    colors[i][0] = r / 255.;
                    colors[i][1] = g / 255.;
                    colors[i][2] = b / 255.;
                }
            }
        }
    } else if (strcmp(format,"ascii")==0) {
        points.resize(3*size);
        for (int i = 0; i < size; i++) {
            fgets(buffer, 99, fp);
            sscanf(buffer, "%lf %lf %lf", &points[i][0], &points[i][1], &points[i][2]);
        }
    } else {
        // fprintf(stderr, "unknown format\n");
        // exit(1);
        throw std::runtime_error("unknown format");
    }

    fclose(fp);
    return;
}
