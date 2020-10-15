// C/C++ headers
#include <cstring>
#include <cassert>

// Athena++ headers
#include "utils.hpp"

void ReadTabular(char const *fname, double** data, int *rows, int *cols)
{
  FILE *fp = fopen(fname, "r");

  // read dimension of the table
  fscanf(fp, "%d%d", rows, cols);

  NewCArray(data, *rows, *cols);

  char buf[256];
  int irow = 0, icol;
  char *p;
  while (fgets(buf, 256, fp) != NULL) {
    p = std::strtok(buf," ,");
    int icol = 0;
    while (p != NULL) {
      sscanf(p, "%lf", &data[irow][icol]);
      p = std::strtok(NULL," ,");
      assert(icol++ < *cols);
    }
    assert(irow++ < *rows);
  }

  fclose(fp);
}
