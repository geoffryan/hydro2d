#ifndef HYDRO2D_IO
#define HYDRO2D_IO

struct parList;
struct grid;

struct io
{
    int nTot;
    int current;
    double tnext;
    double tmin;
    double tmax;
};

const static struct io IO_DEFAULT = {
    .nTot = 0,
    .current = 0,
    .tnext = 0.0,
    .tmin = 0.0,
    .tmax = 0.0
};

void (*io_out)(struct io *, struct grid *, struct parList *);

int set_io(struct parList *par);
void io_init(struct io *iom, struct parList *par);
void io_out_ascii(struct io *iom, struct grid *g, struct parList *par);
void io_out_hdf5(struct io *iom, struct grid *g, struct parList *par);

#endif
