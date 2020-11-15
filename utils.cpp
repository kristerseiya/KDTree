
float distance2(float* x, float* y, int dim) {
  float distance = 0;
  for (int i = 0; i < dim; i++) {
    float a = x[i] - y[i];
    distance += a * a;
  }
  return distance;
}
