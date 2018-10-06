#ifndef MAPPER_H_
#define MAPPER_H_

#include "Pair.h"

#include <string>
#include <list>

using namespace std;

void InitNum(int num_of_var, int &num_of_num, int &num_of_cat, char* type);

void InitFactor(int num_of_var, int num_of_num, int num_of_cat, double * factor, double * &num_factor, double * &cat_factor, char * type);
void DelFactor(double * &num_factor, double * &cat_factor);

void InitDataSet(int num_of_data, int &num_of_num, int &num_of_cat, char* type, double ** &num_data, string ** &cat_data);
void DelDataSet(int num_of_data, double** &num_data, string ** &cat_data);

void ReadData(int num_of_data, int num_of_var, char* type, double ** &num_data, string ** &cat_data, char * inputFilePath);

void InitDistMat(int num_of_data, double ** &dist_mat);
void DelDistMat(int num_of_data, double ** &dist_mat);
void Distance(int num_of_data, int num_of_num, int num_of_cat, double* num_factor, double* cat_factor, double ** num_data, string ** cat_data, double ** &dist_mat);
double DistanceHelper(int num_of_num, int num_of_cat, double* num_factor, double* cat_factor, double * num_row1, double * num_row2, string * cat_row1, string * cat_row2);

void InitFilter(int num_of_data, double * &filter);
void DelFilter(double * &filter);
void Filter(int num_of_data, double ** dist_mat, double * &filter);
void Filter(int num_of_data, double ** dist_mat, double * &filter, double exponent);
double Filter_eccentricity(int num_of_data, double * dist_row, double exponent);

void FindExtreme(int num_of_data, double &min_filt, double &max_filt, double * filter);
void CreateBin(int num_of_data, int num_of_bin, double ratio_of_overlap, double min_filt, double max_filt, double * filter, list<Pair> * &bin);
void DelBin(list<Pair> * &bin);

void CreateCluster(int num_of_bin, double threshold, double ** dist_mat, list<Pair> * &bin, list<list<Pair>> * &cluster, int * &num_of_cluster);
void DelCluster(list<list<Pair>> * &cluster);

void InitClusterInfo(int num_of_bin, int * num_of_cluster, int ** &cluster_info, list<list<Pair>> * &cluster);
void DelClusterInfo(int num_of_bin, int ** &cluster_info);

void DivideCluster(int num_of_bin, list<list<Pair>> * cluster, list<int> ** &cluster_U, list<int> ** &cluster_D);
void DelTypedCluster(int num_of_bin, list<int>** &cluster_T);

void InitGraph(int num_of_bin, int * num_of_cluster, list<int> ** cluster_U, list<int> ** cluster_D, list<int> ** &graph);
bool CompareCluster(list<int> cluster_U, list<int> cluster_D);
void DelGraph(int num_of_bin, list<int> ** &graph);

void MakeDotFile(int num_of_bin, int * num_of_cluster, int ** cluster_info, list<int> ** graph, char * outputFilePath);

#endif
