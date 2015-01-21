/************************************************************************

cmi.cpp
Dominic Smith
Theoretical Systems Biology Group
Division of Molecular Biosciences
Imperial College London

16/07/13

Conditional mutual information based on Cellucci's adaptive
partitioning algorithm from Phys Rev E 71 (2005) pp066208



 ************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <mutinfo.h>
using namespace std;

int main(int argc, char **argv)
{
  string inFile, outPath, miPath;
  // check if we have enough parameters
  if(argc < 5) {
    cout << "Usage is -in <infile> -out <outdir>\n";
    cin.get();
    exit(0);
  }
  else {
    for (int i=1; i<argc; i++) {
      if (i<argc) {
	if (std::string(argv[i])=="-in") {
	  inFile = argv[++i];
	} else if (std::string(argv[i])=="-out") {
	  outPath = argv[++i];
	} else if (std::string(argv[i])=="-mi") {
	  miPath = argv[++i];
	}
      }
      else {
	cout << "Too few or invalid arguments (" << argc << "). Try again.\n" << endl;
	exit(0);
      }
   }
    cout << "Input file: "<< inFile << endl << "Output path: " << outPath << endl;
  }

  ifstream dataFile;
  string s;
  int i=0, j=0;
  vector< vector<float> > varData; // dynamic storage for input data
  vector<string> varNames;

  dataFile.open(inFile.c_str());
  
  if (dataFile==NULL) {
    cout << "Error opening file " << inFile;
  }

  // Read in data from a file 
  while(getline(dataFile,s))
    {
      stringstream linestream(s);
      string varName;
      float data;
      vector<float> row;
      char next;

      getline(linestream, varName, '\t');
      //cout << varName << endl;

      
      
      while (linestream.eof()==false)
	{
	  linestream >> data; //grab the current element from linestream
	  row.push_back(data); //feed that into the row vector
	  j++;
	}

      varNames.push_back(varName);
      varData.push_back(row);
      i++;
    }

  dataFile.close();

  /*
  vector<vector<float> > bounds;

  // Declare partitions class
  partitions parts(varData);
  parts.part_run(bounds, 2);

  vector<float> x = varData[199];
  vector<float> y = varData[200];
  vector<float> x_bounds = bounds[199];
  vector<float> y_bounds = bounds[200];

  xy_data xy(x,y,x_bounds,y_bounds);
  float f = xy.mi_calc();
  bool b = xy.independence_test(0.05);

  cout << f << endl;
  
  if(b==false) {
    cout << "Dependence is not statistically significant at alpha = 0.05" << endl;
  }

  vector<vector<float> > miGrid;
  mutinfo mi(varData, varNames);

  */

  int m = varNames.size();
  vector<vector<float> > outGrid;
  vector<vector<short> > adj;
  float threshold;

  // instantiate the mi class
  mutinfo mi(varData, varNames);
  //threshold = mi.mi_batch_set_threshold(1000); // set the threshold
  // (threshold=0.266051 for p<10^-7 for Basso's data)
  // INSERT method for doing accurate hypothesis testing on the data
  //threshold=0.26;

  // calculate the mi for all nodes
  mi.mi_batch_combs(outGrid);
  mi.writeMiFile(miPath);

  undirected_graph G(outGrid, varNames);
  G.miGraphInference(adj, threshold);

  G.dpi(adj, 0.2); // apply the DPI. A larger threshold means fewer edges will be removed
  G.cmiGraphInference(varData, adj, 7e-2); // apply CMI


  cout << outPath << endl;
  G.writeEdgesFile(outPath);
  

  return 0;
}


  /*

  Auto testing of method
  float adj_werhli[11][11] = {{0,1,0,0,0,0,0,1,1,0,0},
			      {1,0,0,0,0,1,0,1,1,0,0},
			      {0,0,0,1,1,0,0,0,1,0,0},
			      {0,0,1,0,1,0,0,0,1,0,0},
			      {0,0,1,1,0,0,1,0,0,0,0},
			      {0,1,0,0,0,0,1,1,0,0,0},
			      {0,0,0,0,1,1,0,1,0,0,0},
			      {1,1,0,0,0,1,1,0,1,1,1},
			      {1,1,1,1,0,0,0,1,0,1,1},
			      {0,0,0,0,0,0,0,1,1,0,0},
			      {0,0,0,0,0,0,0,1,1,0,0}};

  int tp=0, fp=0, p;

  for (int i=0; i<11; i++) {
    for (int j=0; j<11; j++) {
      if(adj_werhli[i][j]==1) {
	p+=1;
	if (adj[i][j]==1) {
	  tp+=1;
	}
      }
      else {
	if (adj[i][j]==1) {
	  fp+=1;
	}
      }
    }
  }

  cout << endl << "TP = " << tp << ", FP = " << fp << ", P = " << p << endl;

  */
