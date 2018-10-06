#include "Pair.h"
#include "Mapper.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <list>

using namespace std;

void InitNum(int num_of_var, int &num_of_num, int &num_of_cat, char* type) {
	for(int i=0; i<num_of_var; i++) {
		if(type[i] == 'N') num_of_num++;
		else num_of_cat++;
	}
}

void InitFactor(int num_of_var, int num_of_num, int num_of_cat, double * factor, double * &num_factor, double * &cat_factor, char * type) {
	num_factor = new double[num_of_num];
	cat_factor = new double[num_of_cat];
	int num_cnt = 0;
	int cat_cnt = 0;

	for(int i=0; i<num_of_var; i++) {
		if(type[i] == 'N') {
				num_factor[num_cnt] = factor[i];
				num_cnt++;
			}
			else {
				cat_factor[cat_cnt] = factor[i];
				cat_cnt++;
			}
	}
}

void DelFactor(double * &num_factor, double * &cat_factor) {
	delete[] num_factor;
	delete[] cat_factor;
}

void InitDataSet(int num_of_data, int &num_of_num, int &num_of_cat, char* type, double ** &num_data, string ** &cat_data) {
	num_data = new double*[num_of_data];
	cat_data = new string*[num_of_data];

	for(int i=0; i<num_of_data; i++){
		num_data[i] = new double[num_of_num];
		cat_data[i] = new string[num_of_cat];

		for(int j=0; j<num_of_num; j++) num_data[i][j] = 0.0;
		for(int j=0; j<num_of_cat; j++) cat_data[i][j] = "";
	}
}

void DelDataSet(int num_of_data, double ** &num_data, string ** &cat_data) {
	for(int i=0; i<num_of_data; i++) {
		delete[] num_data[i];
		delete[] cat_data[i];
	}
	delete[] num_data;
	delete[] cat_data;
}

void ReadData(int num_of_data, int num_of_var, char* type, double ** &num_data, string ** &cat_data, char * inputFilePath) {
	ifstream ip;
	ip.open(inputFilePath);
	if(!ip.is_open()) cout << "ERROR : File Open" << endl;

	for(int tot_cnt = 0; tot_cnt < num_of_data; tot_cnt++) {
		int num_cnt = 0;
		int cat_cnt = 0;
		string raw;

		for(int i=0; i<num_of_var; i++) {
			if(i!=num_of_var-1)
				getline(ip, raw ,',');
			else
				getline(ip, raw ,'\n');

			if(type[i] == 'N') {
				num_data[tot_cnt][num_cnt] = atof((char*)raw.c_str());
				num_cnt++;
			}
			else {
				cat_data[tot_cnt][cat_cnt] += raw;
				cat_cnt++;
			}
		}
	}
	ip.close();
}

void InitDistMat(int num_of_data, double ** &dist_mat) {
	dist_mat = new double*[num_of_data];

	for(int i=0; i<num_of_data; i++) {
		dist_mat[i] = new double[num_of_data];
		for(int j=0; j<num_of_data; j++) dist_mat[i][j] = 0.0;
	}
}

void DelDistMat(int num_of_data, double ** &dist_mat) {
	for(int i=0; i<num_of_data; i++)
		delete[] dist_mat[i];

	delete[] dist_mat;
}

void Distance(int num_of_data, int num_of_num, int num_of_cat, double * num_factor, double * cat_factor, double ** num_data, string ** cat_data, double ** &dist_mat) {
	for(int i=0; i<num_of_data; i++) {
		for(int j=0; j<i; j++) {
			dist_mat[i][j] = DistanceHelper(num_of_num, num_of_cat, num_factor, cat_factor, num_data[i], num_data[j], cat_data[i], cat_data[j]);
			dist_mat[j][i] = dist_mat[i][j];
		}
	}
}

double DistanceHelper(int num_of_num, int num_of_cat, double* num_factor, double* cat_factor, double * num_row1, double * num_row2, string * cat_row1, string * cat_row2) {
	double dist = 0.0;
	for(int i=0; i<num_of_num ; i++)
		dist += (num_factor[i]*(num_row1[i]-num_row2[i])*(num_row1[i]-num_row2[i]));
	for(int i=0; i<num_of_cat ; i++) {
		if(cat_row1[i].compare(cat_row2[i]) != 0)
			dist += cat_factor[i]*1;
	}
	return dist;
}

void InitFilter(int num_of_data, double * &filter) {
	filter = new double[num_of_data];
}

void DelFilter(double * &filter) {
	delete[] filter;
}

void Filter(int num_of_data, double ** dist_mat, double * &filter) {
	Filter(num_of_data, dist_mat, filter, 1.0);
}

void Filter(int num_of_data, double ** dist_mat, double * &filter, double exponent) {
	for(int i=0; i<num_of_data; i++)
		filter[i] = Filter_eccentricity(num_of_data, dist_mat[i], exponent);
}

double Filter_eccentricity(int num_of_data, double * dist_row, double exponent){
	double filter = 0;
	for(int i=0; i<num_of_data; i++) {
		filter += pow(dist_row[i], exponent);
	}
	filter = pow((filter/((double)num_of_data)), 1/exponent);
	return filter;
}

void FindExtreme(int num_of_data, double &min_filt, double &max_filt, double * filter) {
	min_filt = filter[0];
	max_filt = filter[0];
	for(int i=1; i<num_of_data; i++) {
		if(filter[i] < min_filt) min_filt = filter[i];
		else if(filter[i] > max_filt) max_filt = filter[i];
	}
}

void CreateBin(int num_of_data, int num_of_bin, double ratio_of_overlap, double min_filt, double max_filt, double * filter, list<Pair> * &bin) {
	bin = new list<Pair>[num_of_bin];

	for(int i=0; i<num_of_bin; i++)
		bin[i] = list<Pair>();

	double x = (max_filt-min_filt) / ((num_of_bin-1)*(1-ratio_of_overlap)+1);
	int first_bin;
	int last_bin;

	for(int i=0; i<num_of_data; i++) {
		double temp;
		temp = (filter[i]-min_filt-x)/((1-ratio_of_overlap)*x);
		if(temp > num_of_bin-1) {
			first_bin = num_of_bin - 1;
		} else if(temp <= -1) {
			first_bin = 0;
		} else {
			first_bin = (int) ceil(temp);
		}

		temp = (filter[i]-min_filt)/((1-ratio_of_overlap)*x);
		if(temp >= num_of_bin) {
			last_bin = num_of_bin -1;
		} else if(temp < 0) {
			last_bin = 0;
		} else {
			last_bin = (int) floor(temp);
		}

		if(first_bin == last_bin) {
			bin[first_bin].push_back(Pair(i,'N'));
		} else {
			bin[first_bin].push_back(Pair(i,'U'));
			for(int cnt=first_bin+1; cnt<last_bin; cnt++)
				bin[cnt].push_back(Pair(i,'B'));
			bin[last_bin].push_back(Pair(i,'D'));
		}
	}
}

void DelBin(list<Pair> * &bin) {
	delete[] bin;
}

void CreateCluster(int num_of_bin, double threshold, double ** dist_mat, list<Pair> * &bin, list<list<Pair>> * &cluster, int * &num_of_cluster) {
	cluster = new list<list<Pair>>[num_of_bin];
	num_of_cluster = new int[num_of_bin];

	for(int i=0; i<num_of_bin; i++){
		cluster[i] = list<list<Pair>>();
		cluster[i].push_back(list<Pair>());
		list<list<Pair>>::iterator iter_tot_cluster = cluster[i].begin();
		num_of_cluster[i] = 0;

		while(!bin[i].empty()) {
			//(*cluster[i].rbegin()).push_back(bin[i].front());
			(*iter_tot_cluster).push_back(bin[i].front());
			bin[i].pop_front();
			num_of_cluster[i]++;

			for(list<Pair>::iterator iter_cluster = (*iter_tot_cluster).begin(); iter_cluster != (*iter_tot_cluster).end(); iter_cluster++) {
				list<Pair>::iterator iter_bin = bin[i].begin();
				while(iter_bin != bin[i].end()) {
					if(dist_mat[(*iter_bin).data][(*iter_cluster).data] < threshold) {
						(*iter_tot_cluster).push_back((*iter_bin));
						iter_bin = bin[i].erase(iter_bin);
					} else {
						iter_bin++;
					}
				}

				if(bin[i].empty())	goto out;
			}

			cluster[i].push_back(list<Pair>());
			iter_tot_cluster++;
		}
		out:;
	}
}

void DelCluster(list<list<Pair>> * &cluster) {
	delete[] cluster;
}

void InitClusterInfo(int num_of_bin, int * num_of_cluster, int ** &cluster_info, list<list<Pair>> * &cluster) {
	cluster_info = new int * [num_of_bin];
	for(int i=0; i<num_of_bin; i++) {
		cluster_info[i] = new int[num_of_cluster[i]];

		int cnt = 0;
		for(list<list<Pair>>::iterator iter = cluster[i].begin(); iter != cluster[i].end(); iter++, cnt++) {
			cluster_info[i][cnt] = (*iter).size();
		}
	}
}

void DelClusterInfo(int num_of_bin, int ** &cluster_info){
	for(int i=0; i<num_of_bin; i++)
		delete[] cluster_info[i];

	delete[] cluster_info;
}

void DivideCluster(int num_of_bin, list<list<Pair>> * cluster, list<int> ** &cluster_U, list<int> ** &cluster_D) {
	cluster_U = new list<int>*[num_of_bin];
	cluster_D = new list<int>*[num_of_bin];

	for(int i=0; i<num_of_bin; i++) {
		int num_of_cluster = cluster[i].size();
		cluster_U[i] = new list<int>[num_of_cluster];
		cluster_D[i] = new list<int>[num_of_cluster];

		int cnt = 0;
		for(list<list<Pair>>::iterator iter_tot_cluster = cluster[i].begin(); iter_tot_cluster != cluster[i].end(); iter_tot_cluster++, cnt++) {
			for(list<Pair>::iterator iter_cluster = (*iter_tot_cluster).begin(); iter_cluster != (*iter_tot_cluster).end(); iter_cluster++) {
				cluster_U[i][cnt] = list<int>();
				cluster_D[i][cnt] = list<int>();

				switch((*iter_cluster).type) {
				case 'X': break;
				case 'N': break;
				case 'U':
					cluster_U[i][cnt].push_back((*iter_cluster).data);
					break;
				case 'D':
					cluster_D[i][cnt].push_back((*iter_cluster).data);
					break;
				case 'B':
					cluster_U[i][cnt].push_back((*iter_cluster).data);
					cluster_D[i][cnt].push_back((*iter_cluster).data);
					break;
				}
			}

			cluster_U[i][cnt].sort();
			cluster_D[i][cnt].sort();
			//Sort for later purpose
		}
	}
}

void DelTypedCluster(int num_of_bin, list<int>** &cluster_T) {
	for(int i=0; i<num_of_bin; i++)
		delete[] cluster_T[i];

	delete[] cluster_T;
}

void InitGraph(int num_of_bin, int * num_of_cluster, list<int> ** cluster_U, list<int> ** cluster_D, list<int> ** &graph) {
	graph = new list<int> * [num_of_bin-1];

	for(int i=0; i<num_of_bin-1; i++) {
		graph[i] = new list<int> [num_of_cluster[i]];

		for(int j=0; j<num_of_cluster[i]; j++) {
			for(int k=0; k<num_of_cluster[i+1]; k++) {
				if(CompareCluster(cluster_U[i][j], cluster_D[i+1][k]))
					graph[i][j].push_back(k);
			}
		}
	}
}

bool CompareCluster(list<int> cluster_U, list<int> cluster_D) {
	list<int>::iterator iter_U = cluster_U.begin();
	list<int>::iterator iter_D = cluster_D.begin();

	while((iter_U != cluster_U.end()) && (iter_D != cluster_D.end())) {
		if((*iter_U) < (*iter_D)) {
			iter_U++;
		} else if((*iter_U) > (*iter_D)) {
			iter_D++;
		} else {
			return true;
		}
	}
	return false;
}

void DelGraph(int num_of_bin, list<int> ** &graph){
	for(int i=0; i<num_of_bin-1; i++)
		delete[] graph[i];

	delete[] graph;
}

void MakeDotFile(int num_of_bin, int * num_of_cluster, int ** cluster_info, list<int> ** graph, char * outputFilePath) {
	ofstream op;
	op.open(outputFilePath);
	if(!op.is_open()) cout << "ERROR : File Open" << endl;

	op << "digraph G{" << endl;
	op << " edge[style=solid, len=2]" << endl;
	op << " model = mds { " << endl;
	op << "  node[shape=circle, style=filled]" << endl;

	for(int i=0; i<num_of_bin; i++) {
		for(int j=0; j<num_of_cluster[i]; j++) {
			op << "  A" << i << "_" << j << "[fillcolor=" << "\"" << 0.5*i/num_of_bin << " 0.914 " << "0.90" << "\"";
			op << " ,width=" << 0.5*sqrt((double)cluster_info[i][j]) << "]" << endl;
		}
	}

	op << " }" << endl;

	for(int i=0; i<num_of_bin-1; i++) {
		for(int j=0; j<num_of_cluster[i]; j++) {
			for(list<int>::iterator iter = graph[i][j].begin(); iter != graph[i][j].end(); iter++)
				op << " A" << i << "_" << j << " -> " <<  "A" << i+1 << "_" << (*iter) << endl;
		}
	}

	op << "}";
	op.close();
}


int main() {
//-----------------------------------------------------------------------------
//-- Need to change -------------------------------------------------

	int num_of_data = 4000;
	int num_of_var = 9;
	char type[num_of_var];
	double factor[num_of_var];
	double exponent = 2.0;
	int num_of_bin = 10;
	double ratio_of_overlap = 0.8;
	double threshold = 2;
	char inputFilePath[] = "./Data/check.txt";
	char outputFilePath[] = "./Result/result.txt";

	for(int i=0; i<num_of_var; i++) {
		type[i] = 'N';
		factor[i] = 1;
	}
	type[0] = 'C';
	/* N: Numerical, C: Categorical */

//------------------------------------------------------------------------------
//-- Do not change below ---------------------------------------------

	int num_of_num = 0;
	int num_of_cat = 0;
	double * num_factor;
	double * cat_factor;
	double ** num_data;
	//char *** cat_data;
	string ** cat_data;
	double ** dist_mat;
	double * filter;
	double min_filt;
	double max_filt;
	list<Pair>* bin;
	list<list<Pair>>* cluster;
	int * num_of_cluster;
	int ** cluster_info; // For visualization
	list<int> ** cluster_U;
	list<int> ** cluster_D;
	list<int> ** graph;

	InitNum(num_of_var, num_of_num, num_of_cat, type);
	InitFactor(num_of_var, num_of_num, num_of_cat, factor, num_factor, cat_factor, type);
	InitDataSet(num_of_data, num_of_num, num_of_cat, type, num_data, cat_data);
	ReadData(num_of_data, num_of_var, type, num_data, cat_data, inputFilePath);
	InitDistMat(num_of_data, dist_mat);
	Distance(num_of_data, num_of_num, num_of_cat, num_factor, cat_factor, num_data, cat_data, dist_mat);
	DelDataSet(num_of_data, num_data, cat_data);
	InitFilter(num_of_data, filter);
	Filter(num_of_data, dist_mat, filter, exponent);
	FindExtreme(num_of_data, min_filt, max_filt, filter);
	CreateBin(num_of_data, num_of_bin, ratio_of_overlap, min_filt, max_filt, filter, bin);
	DelFilter(filter);
	CreateCluster(num_of_bin, threshold, dist_mat, bin, cluster, num_of_cluster);
	InitClusterInfo(num_of_bin, num_of_cluster, cluster_info, cluster);
	DelDistMat(num_of_data, dist_mat);
	DelBin(bin);
	DivideCluster(num_of_bin, cluster, cluster_U, cluster_D);
	DelCluster(cluster);
	InitGraph(num_of_bin, num_of_cluster, cluster_U, cluster_D, graph);
	DelTypedCluster(num_of_bin, cluster_U);
	DelTypedCluster(num_of_bin, cluster_D);
	MakeDotFile(num_of_bin, num_of_cluster, cluster_info, graph, outputFilePath);
	DelGraph(num_of_bin, graph);
	DelClusterInfo(num_of_bin, cluster_info);

	cout << "END" << endl;
	return 0;
}
