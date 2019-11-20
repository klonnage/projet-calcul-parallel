Mat fsta(double Lx, double Ly, int Nx, int Ny) {
  double dx = Lx / (double)Nx;
  double dy = Ly / (double)Ny;
  Mat result(Nx, Ny, 0.);

  for (size_t i = 0; i < Nx; ++i) {
    for (size_t j = 0; j < Ny; ++j) {
      double x = (double) i * dx;
      double y = (double) j * dy;
      result(i, j) = (x - x*x + y - y*y);
    }
  }
  
  return result;
}

Mat gsta(double Lx, double Ly, int Nx, int Ny) {
  return Mat(Nx, Ny, 0.);
}

Mat hsta(double Lx, double Ly, int Nx, int Ny) {
  return Mat(Nx, Ny, 0.);
}

Mat fper(double Lx, double Ly, int Nx, int Ny) {

}

Mat gper(double Lx, double Ly, int Nx, int Ny) {

}

Mat hper(double Lx, double Ly, int Nx, int Ny) {

}
