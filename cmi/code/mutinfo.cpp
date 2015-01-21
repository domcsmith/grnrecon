#include <mutinfo.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <ctime>
using namespace std;

/********************************************************
partitions class functions
 ********************************************************/

// Constructor for partitions
partitions::partitions(vector< vector<float> >& varData) {
  inData=varData; // copy varData to class member inData
}

partitions::~partitions() {
}

/* 
   Function to find the boundaries which partition each data
   variable into equiprobable bins subject to Cochran's criterion
*/
void partitions::part_run(vector< vector<float> >& bounds, int dim) {

  int N_p = inData.size();
  int N_d = inData[0].size();
  int N_e, pos;
  vector<float> row;

  N_e = (int)floor(pow(N_d/5, 1/(double)dim));

  for (int i=0; i<N_p; i++) {

    // Destroy the contents of bounds for the next probe
    row.clear();

    // Sort each probe's data
    partitions::quickSort(inData[i], 0, N_d);
    
    // get the partition boundaries and output them to "bounds"
    for (int j=0; j<=N_e; j++) {
      pos = (int)(j * round(N_d/N_e));
      row.push_back(inData[i][pos]);
    }

    bounds.push_back(row);
  }
  
}

// Quick sort to aid partitioning of each vector
void partitions::quickSort(vector<float>& vec, int left, int right) {
  
  int i = left, j = right;
  float tmp;
  float pivot = vec[(left + right) / 2];

  /* partition */

  while (i <= j) {
    while (vec[i] < pivot)
      i++;
    while (vec[j] > pivot)
      j--;
    if (i <= j) {
      tmp = vec[i];
      vec[i] = vec[j];
      vec[j] = tmp;
      i++;
      j--;
    }
  }
      
  /* recursion */
  if (left < j) quickSort(vec, left, j);
  if (i < right) quickSort(vec, i, right);

}

void mutinfo::quickSort(vector<float>& vec, int left, int right) {
  
  int i = left, j = right;
  float tmp;
  float pivot = vec[(left + right) / 2];

  /* partition */

  while (i <= j) {
    while (vec[i] < pivot)
      i++;
    while (vec[j] > pivot)
      j--;
    if (i <= j) {
      tmp = vec[i];
      vec[i] = vec[j];
      vec[j] = tmp;
      i++;
      j--;
    }
  }
      
  /* recursion */
  if (left < j) quickSort(vec, left, j);
  if (i < right) quickSort(vec, i, right);

}

/********************************************************
mi_calc class functions
 ********************************************************/

// Constructor for xy_data
xy_data::xy_data (vector<float>& var1, vector<float>& var2, vector<float>& x_bounds, vector<float>& y_bounds) {
  x = var1;
  y = var2;
  x_b = x_bounds;
  y_b = y_bounds;
  n = x_b.size()-1;
  bin_XY.assign(n, vector<int>(n, 0)); // allocate and initialize joint frequency grid
  bin_X.assign(n, 0); // allocate marginal frequency
  bin_Y.assign(n, 0);
}

// Destructor for xy_data
xy_data::~xy_data () {
}

// Function to calculate the mutual information
float xy_data::mi_calc() {
  int x_it, y_it;
  float mi=0, t1=0, t2=0, t3=0;

  data_length = x.size();

  for (int i=0; i<x.size(); i++) {
    x_it = distance(x_b.begin(), upper_bound(x_b.begin(), x_b.end(), x[i])); // get the bin number for the data point x[i]
    y_it = distance(y_b.begin(), upper_bound(y_b.begin(), y_b.end(), y[i])); // get the bin number for the data point y[i]
    x_it = min(x_it,n) -1;
    y_it = min(y_it,n) -1;
    bin_XY[x_it][y_it]++; // increment joint frequency
    bin_X[x_it]++; // increment marginal frequency
    bin_Y[y_it]++; // increment marginal frequency
  }


  for (int i=0; i<n; i++) {
    //cout << bin_X[i] << ", ";
  }

  //cout <<endl;

  for (int i=0; i<n; i++) {
    //cout << bin_Y[i] << ", ";
  }

  //cout << endl;

  for (int i=0; i < n; i++) {
    for (int j=0; j< n; j++) {

      //cout << bin_XY[i][j] << ", ";

      t1 = (float)bin_XY[i][j] / data_length;
      if(t1 > 0) {
	t2 = (float)bin_XY[i][j] / ((float)bin_X[i]);
	t3 = (float)bin_Y[j] / data_length;
	mi += t1 * log2 (t2/t3);	
      }
    }
    //cout << endl;
  }

  return mi;
}

bool xy_data::independence_test (float sig_level) {

  float nu, chisq=0, E_XY, p_null;
  bool retval=false;

  // calculate chi^2
  for (int i=0; i<bin_X.size(); i++) {
    for (int j=0; j<bin_Y.size(); j++) {
      E_XY = (float)bin_X[i] * (float)bin_Y[j] / data_length;
      chisq += pow((float)bin_XY[i][j] - E_XY, 2) / E_XY;
    }
  }
  
  //calculate degrees of freedom
  nu = (bin_X.size() - 1) * (bin_Y.size() - 1);
  
  // integrate upper incomplete gamma function Q(nu/2, chisq/2)
  // to yield probability of the null hypothesis
  p_null = gsl_sf_gamma_inc_P(nu/2,chisq/2);

  if (p_null > 1-sig_level) {
    retval = true;
  }

  cout << chisq << "\t" << nu << endl;
  cout << "P_null = " << p_null << endl;

  return retval;  
}


/********************************************************
mutinfo class functions
 ********************************************************/

// Constructor for mutinfo
mutinfo::mutinfo (vector<vector<float> >& varData, vector<string>& varNames) {
  inData = varData;
  inNames = varNames;
}

// Destructor for mutinfo
mutinfo::~mutinfo () {
}

// Performs MI calculation for all combinations of variables in the
// input data
void mutinfo::mi_batch_combs(vector<vector<float> >& outData) {

  vector<vector<float> > bounds;
  partitions parts(inData);
  parts.part_run(bounds, 2);
  float f=0;

  for(int i=0; i<inNames.size(); i++) {
    vector<float> row;
    for(int j=i+1; j<inNames.size(); j++) {
      xy_data xy(inData[i], inData[j], bounds[i], bounds[j]);
      f = xy.mi_calc();
      row.push_back(f);

      //cout << "MI(" << inNames[i] << ", " << inNames[j] << ") = " << f << endl;
    }
    miGrid.push_back(row);
  }
  
  outData = miGrid;

  /*
  for (int i=0; i< miGrid.size(); i++) {
    for (int j=0; j<miGrid[i].size(); j++) {
      cout << miGrid[i][j] << ", ";
    }
    cout << endl;
  }
  */
  
}

void mutinfo::writeMiFile(const string& outPath) {
  ofstream outFile;
  outFile.open(outPath.c_str());

  for (int i=0; i<inNames.size(); i++) {
    for (int j=i; j<inNames.size()-1; j++) {

      //cout << inNames[i] << "\t" << inNames[j] << "\t" << miGrid[i][j] << "\n";
      outFile << inNames[i] << "\t" << inNames[j] << "\t" << miGrid[i][j] << "\n";
    }
  }

  cout << "Done";
  outFile.close();
}


/*
// Runs a batch calculation of MI for a sequence of jointly Gaussian pairs of 
// variables in the input data. Runs from correlation k=0 to k=0.99 in 0.01
// increments. Same variables run 100 times each.
void mutinfo::mi_batch_pairs(vector<vector<float> >& outData) {

  vector<vector<float> > bounds;
  vector<float> row;
  partitions parts(inData);
  parts.part_run(bounds, 2);
  float f, a=0, mi_ana, e_n=0, e_d=0;
  ofstream myfile;
  myfile.open("mi_out.dat");

  for(int k=0; k<100; k++) {
    mi_ana = (-0.5)*log(1-(float)k/100);
    a=0;
    e_n=0;
    e_d=0;
    for(int i=k*200; i<(k+1)*200;i+=2) {
      xy_data xy(inData[i], inData[i+1], bounds[i], bounds[i+1]);
      f = xy.mi_calc();
      row.push_back(f);

      a+=f;
      e_n+=pow((f-mi_ana),2);
      e_d+=pow(mi_ana,2);

      //cout << "MI(" << inNames[i] << ", " << inNames[i+1] << ") = " << f << endl;
    }

    a=a/100;

    cout << "Analytic MI = " << mi_ana << endl;
    cout << "Average MI = " << a << endl;
    cout << "Error vs Analytic MI = " << e_n/e_d << endl;

    myfile << a << "\t";
  
    outData = miGrid;
  }
  myfile.close();
}
*/

float mutinfo::mi_batch_set_threshold(int numSamples) {

  float f;
  vector<float> out;

  for (int k=0; k<numSamples; k++) {
    vector<vector<float> > bounds;
    for (int i=0; i<inNames.size(); i+=2) {
      randomize(inData[i+1]);
    }

    partitions parts(inData);
    parts.part_run(bounds, 2);

    for(int i=0; i<inNames.size()-2;i+=2) {
      xy_data xy(inData[i], inData[i+1], bounds[i], bounds[i+1]);
      f = xy.mi_calc();
      out.push_back(f);
    }
  }

  quickSort(out, 0, out.size());
  cout << out[0] << ", " << out[out.size()-1] << "," << out.size() << endl;
  
  return out[out.size()-1];
}


void mutinfo::mi_batch_pairs(vector<vector<float> >& outData) {

  vector<vector<float> > bounds;
  vector<float> row;
  partitions parts(inData);
  parts.part_run(bounds, 2);
  float f;
  ofstream myfile;
  myfile.open("mi_pairs_out.dat");

  for(int i=0; i<inNames.size();i+=2) {
    xy_data xy(inData[i], inData[i+1], bounds[i], bounds[i+1]);
    for (int j=0; j<10; j++) {
      f = xy.mi_calc();
      myfile << inNames[i] << "\t" << inNames[i+1] << "\t" << f << endl;
    }
  }
  
  outData = miGrid;
  myfile.close();
}

// Constructor for xyz_data
xyz_data::xyz_data (vector<float>& var1, vector<float>& var2, vector<float>& var3, vector<float>& x_bounds, vector<float>& y_bounds, vector<float>& z_bounds) {
  x = var1;
  y = var2;
  z = var3;
  x_b = x_bounds;
  y_b = y_bounds;
  z_b = z_bounds;
  n = x_b.size()-1;
  bin_XYZ.assign(n, vector<vector<int> >(n, vector<int>(n, 0))); // allocate and initialise three-way joint frequency grid
  bin_XZ.assign(n, vector<int>(n, 0)); // allocate and initialize joint frequency grids
  bin_YZ.assign(n, vector<int>(n, 0));
  bin_Z.assign(n, 0); // allocate marginal frequency for Z
}

// Destructor for xyz_data
xyz_data::~xyz_data () {
}

// Function to calculate the conditional mutual information
float xyz_data::cmi_calc() {
  int x_it, y_it, z_it;
  float cmi=0, t1=0, t2=0, t3=0;

  data_length = x.size();

  for (int i=0; i<x.size(); i++) {
    x_it = distance(x_b.begin(), upper_bound(x_b.begin(), x_b.end(), x[i])); // get the bin number for the data point x[i]
    y_it = distance(y_b.begin(), upper_bound(y_b.begin(), y_b.end(), y[i])); // get the bin number for the data point y[i]
    z_it = distance(z_b.begin(), upper_bound(z_b.begin(), z_b.end(), z[i])); // get the bin number for data point z[i]
    x_it = min(x_it,n) -1;
    y_it = min(y_it,n) -1;
    z_it = min(z_it,n) -1;
    bin_XYZ[x_it][y_it][z_it]++; // increment joint frequency
    bin_XZ[x_it][z_it]++; // increment marginal frequency
    bin_YZ[y_it][z_it]++; // increment marginal frequency
    bin_Z[z_it]++;
  }

  /*
  for (int i=0; i<n; i++) {
    //cout << bin_X[i] << ", ";
  }

  //cout <<endl;

  for (int i=0; i<n; i++) {
    //cout << bin_Y[i] << ", ";
  }

  //cout << endl;
  */
  
  for (int i=0; i < n; i++) {
    for (int j=0; j < n; j++) {
      for (int k=0; k < n; k++) {

	//cout << bin_XY[i][j] << ", ";

	t1 = (float)bin_XYZ[i][j][k] / data_length;
	if(t1 > 0) {
	  t2 = t1 * ((float)bin_Z[k] / data_length);
	  t3 = ((float)bin_XZ[i][k]) * ((float)bin_YZ[j][k]) / pow(data_length, 2);
	  cmi += t1 * log2 (t2/t3);	
	}
      }
    }
    //cout << endl;
  }

  return cmi;
}

/********************************************************************
BATCH CMI CALCULATIONS
 *******************************************************************/

  /*
    1) Loop through a grid of previously calculated MI values I(X,Y)
    2) Only those values with an MI greater than the threshold level included
    3) If significance thresholding is done to the MI grid in an earlier step,
    that can be included by setting them equal to zero in the MI grid
    4) Calculate CMI against all other nodes
    5) Return a 3D vector of all CMIs where position X=i, Y = i+j, Z=k
   
void mutinfo::cmi_batch(vector<vector<float> >& miIn, vector<vector<vector<float> > >& cmiOut) {
  
  vector<vector<float> > bounds;
  partitions parts(inData);
  parts.part_run(bounds, 3); // 3D so use cube root of N_d/5

  vector<vector<float> > plane;
  vector<float> row;
  float cmi;

  for (int i=0; i < miIn.size(); i++) {
    for (int j=0; j < miIn[i].size(); j++) {
      if (miIn[i][j] > 0 && miIn[i][j] > thresholdMI) {
	for (int k=0; k < inData.size(); k++) {
	  xyz_data xyz(inData[i], inData[i+j], inData[k], bounds[i], bounds[i+j], bounds[k]);

	  cmi = xyz.cmi_calc();
	  row.push_back(cmi);
	}
	else {
	  row.push_back(0);
	}
      }
      plane.push_back(row);
    }
    cmiGrid.push_back(plane);
  }

  cmiOut = cmiBatch;
}
  */
/*******************************************************************
NETWORK INFERENCE BASED ON MI
 ******************************************************************/

undirected_graph::undirected_graph(vector<vector<float> >& miGrid, vector<string>& varNames) {
  mi = miGrid;
  inNames = varNames;
  n = varNames.size();
  adj.assign(n, vector<short>(n, 0)); // initialise adjacency matrix
  adj_d.assign(n, vector<short>(n, 0));
}

undirected_graph::~undirected_graph() {
}

void undirected_graph::miGraphInference(vector<vector<short> >& outAdj, float threshold) {

  cout << n << endl;

  // populate adj with edges
  for (int i=0; i < n; i++) {
    vector<short> row;
    for (int j=0; j < mi[i].size(); j++) {
      if (mi[i][j] > threshold) {
	adj[i][i+j+1] = 1;
	adj[1+i+j][i] = 1; // because adj is symmetric
	row.push_back(1);
      }
      else {
	row.push_back(0);
      }
    }
    adj_t.push_back(row);
  }  
  outAdj = adj; // return outAdj

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      cout << adj[i][j] << ",";
    }
    cout << endl;
  }

  cout<<endl;

  // adj_t is a triangular adjancency matrix
  for (int i=0; i<adj_t.size(); i++) {
    for (int j=0; j<adj_t[i].size(); j++) {
      cout << adj_t[i][j] << ",";
    }
    cout << endl;
  }

  cout<<endl;
}

// implement the data processing inequality
void undirected_graph::dpi(vector<vector<short> >& outAdj, float threshold) {
  cout << n << endl;

  // loop through all A(X,Y). If A(X,Y)==1, loop through all other nodes to see if X is connected to another node, e.g. Z. If it is, check if Y is also connected to Z. If it is, remove the edge with the lowest MI provided that MI*(1+threshold) is still lower than the minimum of the other two.

  float i_xy=0, i_xz=0, i_yz=0;
  vector<float> comp;

  // initialise the DPI-adjusted adjacency matrix
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      adj_d[i][j]=adj[i][j];
    }

  for (int i=0; i<adj.size()-2; i++) {
    for (int j=i+1; j<adj[i].size()-1; j++) {
      if (j<n && adj[i][j]==1) {
	for (int k=j+1; k<adj[i].size(); k++) {
	  // j < k < i so i should always be first index, then k, then j
	  if (k<n && adj[i][k]==1 && adj[j][k]==1) {
	    i_xy=mi[i][j-i-1];
	    i_xz=mi[i][k-i-1];
	    i_yz=mi[j][k-j-1];

	    if (i_xy < i_xz) {
	      if (i_xy < i_yz) {
		// i_xy is smallest
		if (i_xy*(1+threshold)<fminf(i_xy,i_yz)) {
		  adj_d[i][j]=0;
		  adj_d[j][i]=0;
		}
	      }
	      else {
		if (i_xz < i_yz) {
		  // i_xz is smallest
		  if (i_xz*(1+threshold)<fminf(i_yz,i_xy)) {
		    adj_d[i][k]=0;
		    adj_d[k][i]=0;
		  }
		}
		else {
		  // i_yz is smallest
		  if (i_yz*(1+threshold)<fminf(i_xy,i_xz)) {
		    adj_d[j][k]=0;
		    adj_d[k][j]=0;
		  }
		}
	      }
	    }
	    else {
	      if (i_xz < i_yz) {
		// i_xz is smallest
		if (i_xz*(1+threshold)<fminf(i_xy,i_yz)) {
		    adj_d[i][k]=0;
		    adj_d[k][i]=0;
		}
	      }
	      else {
		// i_yz is smallest
		if (i_yz*(1+threshold)<fminf(i_xy,i_xz)) {
		    adj_d[j][k]=0;
		    adj_d[k][j]=0;
		}
	      }
	    }
	    
	  }
	}
      }
    }
  }


  outAdj = adj_d; // return outAdj

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      cout << adj_d[i][j] << ",";
    }
    cout << endl;
  }

  cout<<endl;
}


void undirected_graph::cmiGraphInference(vector<vector<float> >& varData, vector<vector<short> >& outAdj, float zero_threshold) {
  vector<vector<float> > bounds;
  partitions parts(varData);
  parts.part_run(bounds, 3);
  float cmi_xyz, cmi_xzy, cmi_yzx;

  // initialise the cmi-adjusted adjacency matrix as the mi adjacency matrix
  adj_c.assign(n, vector<short>(n,0));

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      adj_c[i][j]=adj_d[i][j];
    }


  for (int i=0; i<n; i++) {
    for (int j=i; j<n; j++) {
      // rule out any I(X,X) situations
      if (i!=j) {
	if (adj_d[i][j]==1) {
	  for (int k=0; k<n; k++) {
	    // rule out any I(X,Y|X) and I(X,Y|Y) situations
	    if (k!=i && k!=j) {
	      // calculate I(X,Y|Z)
	      xyz_data xyz(varData[i], varData[j], varData[k], bounds[i], bounds[j], bounds[k]); // instantiate class for I(X,Y|Z)
	      cmi_xyz = xyz.cmi_calc();

	      //cout << "I(" << inNames[i] << "," << inNames[j] << "|" << inNames[k] << ") = " << cmi_xyz << endl;

	      if (cmi_xyz < zero_threshold) {
		adj_c[i][j]=0;
		adj_c[j][i]=0;
	      }
	    }
	  }
	}
      }
    }
  }

  /*

  // algorithm to remove
  for (int i=0; i<n; i++) {
    for (int j=i; j<n; j++) {

      if (i==j) {
	// do nothing on the diagonal of the adjacency matrix
	// move to the next row of adj if i has reached the final element
	// in that row
      }
      else {      
	// find an edge between nodes i and j
	if (adj[i][j]==1) {
	  // if the node exists, look for an edge between nodes i and k
	  // if there are none, move to the next i node
	  for (int k=j+1; k<n-i; k++) {
	    if (k<n) {
	      if (adj[i][k]==1) {
		if (adj[j][k]==1) {
		  xyz_data xyz(varData[i], varData[j], varData[k], bounds[i], bounds[j], bounds[k]); // instantiate class for I(X,Y|Z)
		  xyz_data xzy(varData[i], varData[k], varData[j], bounds[i], bounds[k], bounds[j]); // instantiate class for I(X,Z|Y)
		  xyz_data yzx(varData[j], varData[k], varData[i], bounds[j], bounds[k], bounds[i]); // instantiate class for I(Y,Z|X)

		  // compute the conditional mutual information
		  cmi_xyz = xyz.cmi_calc();
		  cmi_xzy = xzy.cmi_calc();
		  cmi_yzx = yzx.cmi_calc();
	      
		  // case 1: coregulation of Y and Z by X
		  if (cmi_xyz < zero_threshold) {
		    if (cmi_xzy < zero_threshold) {
		      if (cmi_yzx < lt_threshold*mi[j][k]) {
			// might need to make this mi[k][j] (lower/upper triangular)
			adj_c[j][k]=0;
			adj_c[k][j]=0;
		      }
		    }
		  }

		  // case 2: coregulation of X and Y by Z
		  if (cmi_yzx < zero_threshold) {
		    if (cmi_xzy < zero_threshold) {
		      if (cmi_xyz < lt_threshold*mi[i][j]) {
			// might need to make this mi[k][j] (lower/upper triangular)
			adj_c[i][j]=0;
			adj_c[j][i]=0;
		      }
		    }
		  }

		  // case 3: coregulation of X and Z by Y
		  if (cmi_yzx < zero_threshold) {
		    if (cmi_xyz < zero_threshold) {
		      if (cmi_xzy < lt_threshold*mi[i][k]) {
			// might need to make this mi[k][j] (lower/upper triangular)
			adj_c[i][k]=0;
			adj_c[k][i]=0;
		      }
		    }
		  }
		
		}
	      }
	    }
	  }
	}
	else {
	  for (int k=i+1; k<n-i; k++) {
	    if (adj[i][k]==0) {
	      if (adj[j][k]==0) {
		// calculate the CMI
	      }
	    }
	  }
	}
      }
    }
  }

  */

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      cout << adj_c[i][j] << ",";
    }
    cout << endl;
  }

  outAdj=adj_c;
  
}


void undirected_graph::writeEdgesFile(const string& outPath) {
  ofstream outFile;
  outFile.open(outPath.c_str());
  
 
  for (int i=0; i < n; i++) {
    for (int j=i; j < n; j++) {
      if (adj_c[i][j]==1) {
	outFile << inNames[i] << "\t" << inNames[j] << "\n";
      }
      
    }
  }
  outFile.close();
}

void mutinfo::randomize(vector<float>& rv) {
  std::srand (unsigned(std::time(0)));
  std::random_shuffle (rv.begin(), rv.end());

}
