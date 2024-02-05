double ker(double x, double y) {
  double out, pi;
  pi = 3.14159265;
  out = 0.0;
  if ((x*x + y*y) <= 1.0){
    out = 2.0/pi * (1.0 - x*x - y*y);
  }
  return out;
}
