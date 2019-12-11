struct InputData {
    unsigned int rowCount, colCount;/* problem dimensions */
    float Lrow, Lcol;
    float D;
    float dt;
    float tMax;
    float beta;
    int kmax;
    float eps;
    int mode;
};

InputData ReadInput(const char *file);