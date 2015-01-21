#ifndef __MUTINFO_H
#define __MUTINFO_H

#include <vector>
#include <string>
using namespace std;

class partitions {
private:
  vector<vector<float> > inData;
  vector<vector<float> > bounds;
  void quickSort(vector<float>&, int, int);
public:
  partitions(vector<vector<float> >&);
  ~partitions();
  void part_run(vector<vector<float> >&, int);
};

class mutinfo {
private:
  vector<vector<float> > inData;
  vector<string> inNames;
  vector<vector<float> > miGrid;
  vector<vector<vector<float> > > cmiGrid;
  void quickSort(vector<float>&, int, int);
  void randomize(vector<float>&);
public:
  mutinfo(vector<vector<float> >&, vector<string>&);
  ~mutinfo();
  void mi_batch_combs(vector<vector<float> >&);
  void mi_batch_pairs(vector<vector<float> >&);
  float mi_batch_set_threshold(int);
  void cmi_batch(vector<vector<float> >&, vector<vector<vector<float> > >&);
  void writeMiFile(const string& outPath);
};

class xy_data {
private:
  vector<float> x, y, x_b, y_b;
  vector<vector<int> > bin_XY; // <-- Joint distribution
  vector<int> bin_X; // <-- Marginal distributions
  vector<int> bin_Y;
  int n, data_length;
public:
  xy_data(vector<float>&, vector<float>&, vector<float>&, vector<float>&);
  ~xy_data();
  float mi_calc();
  bool independence_test(float);
};

class xyz_data {
private:
  vector<float> x, y, z, x_b, y_b, z_b;
  vector<vector<vector<int> > > bin_XYZ;
  vector<vector<int> > bin_XZ, bin_YZ;
  vector<int> bin_Z;
  int n, data_length;
public:
  xyz_data(vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&, vector<float>&);
  ~xyz_data();
  float cmi_calc();
};

class undirected_graph {
private:
  vector<vector<float> > mi;
  vector<vector<short> > adj, adj_d, adj_c, adj_t;
  vector<string> inNames;
  int n;
public:
  undirected_graph(vector<vector<float> >&, vector<string>&);
  ~undirected_graph();
  void miGraphInference(vector<vector<short> >&, float);
  void dpi(vector<vector<short> >&, float);
  void cmiGraphInference(vector<vector<float> >&, vector<vector<short> >&, float);
  void writeEdgesFile(const string&);
};

#endif /* __MUTINFO_H */
