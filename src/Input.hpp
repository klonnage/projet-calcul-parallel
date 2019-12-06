struct InputData {
    unsigned int Nx, Ny;/* problem dimensions */
    float Lx, Ly;
    float D;
    float dt;
    float tMax;
    float beta;
    int kmax;
};

InputData ReadInput(const char *file);