#include "../SparseMatrix/SparseMatrix.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <ctime>
#include <stdlib.h>
using namespace std;

//--------------First functional block: Do the node selection----------------
double matrixSum(SparseMatrix* matrix, int k, int omitIndex, Edge* deleteEdge);
vector<pair<int, double>>* katzGreedy(SparseMatrix* matrix, int k, int num, int M, bool ifDelta, bool submodule = true);
double deltaMatrixSum(SparseMatrix* matrix, int k, int omitIndex, vector<vector<double>*>* base);
vector<vector<double>*>* matrixSumAllValue(SparseMatrix* matrix, int k);
vector<pair<int, double>>* katzGreedyWithPriorSet(SparseMatrix* matrix,int k, int num, int priorNum, bool ifDelta, bool submodule = true);
void nodeSelect();// Find those top nodes and save them into the local.

//--------------Second functional block: Generate a graph with top or not label-------------------
void labelNetworkEdge(string address, SparseMatrix* edgeSet, SparseMatrix* matrix);
void labelNetworkNode(string address, set<int>* nodeSet, set<int>* rawNodeSet, vector<string>* dictionary);

//-------------Thrid functional block: map a node list with their name-------------
void nodesToEdges(set<int>* nodeSet, SparseMatrix* matrix, string saveAddress);
void edgesToNodes(SparseMatrix* edgeSet, string saveAddress);

//-------------Forth functional block: input and output function-------------------
SparseMatrix* inputNetwork(string filename, bool ifSelf = false); //This file should have the nodeNum in the first line, edgeNum in the second line,and edges with their weight. Input is an actual.ifSelf is a bool value to indicate whether this network has a self ring.
void saveToLocal(vector<pair<int, double>>* result, string address); // Input is a vector of pair, int is the index of node and double is the remaining part in the network, address should be the actual address.
vector<string>* inputDictionary(string dictionary);
set<int>* inputNodeSet(string address);
void saveToLocal(pair<vector<int>*, vector<pair<int, double>>*>* result, string address);

//-------------Fifth functional block: Run a few baseline------------
vector<pair<int, double>>* pageRank(SparseMatrix* matrix, int resultNum, double alpha, double threshold);
vector<int>* diversedRank(SparseMatrix* matrix, int resultNum, double alpha, double threshold, vector<double>* p = NULL);
vector<pair< int, double>>* connecivity(SparseMatrix* matrix, int N);//greedy
vector<int>* maximalIndepSet(SparseMatrix* matrix);

//-------------Sixth function block: Run a heuristics baseline------------
vector<pair<int, double>>* katzGreedyWithAppoximateCriteria(SparseMatrix* matrix, int k, int size);

//--------------Tool function-------------
template <class Type>
string numToString(Type value);


void backbone();
void baseline_diversifiedRank(){};
void baseline_pageRank(){};
void baseline_connectivity();
void self(){
	clock_t begin = clock();
	for(int i =0; i<10000; i++){
		cout << i << endl;
	}
	clock_t end = clock();
	cout<<"Running time: "<<(double)(end-begin)/CLOCKS_PER_SEC*1000<<"ms"<<endl;
}

int main(){
	cout << "enter one algorithm: 1.backbone; 2.baseline_diversifiedRank; 3.baseline_pageRank; 4.baseline_connectivity; 5.self" <<endl;

	int algo_flag;
	cin >> algo_flag;

	if(algo_flag == 1) backbone();
	else if(algo_flag == 2) baseline_diversifiedRank();
	else if(algo_flag == 3) baseline_pageRank();
	else if(algo_flag == 4) baseline_connectivity();
	else if(algo_flag == 5) self();
	else{
		cout << "input error" <<endl;
		exit(0);
	}
}

//-------------First block----------------
//
//
//-------------Do some prelimilaries for the setting and save them to the local-------------
void backbone(){
	cout << "Enter the Dataset name, name, k, alpha, N, M, ifCS, ifDelta, ifSub" << endl;
	string dataset, name;
	int k, N, M;
	double alpha;
	bool ifDelta, ifSub, ifCS;

	cin >> dataset >> name >> k >> alpha >> N >> M >> ifCS >> ifDelta >> ifSub;

	string inputAddress = "../../Dataset/" + dataset + "/indexNetwork.txt";
	string saveAddress = "../../Output/" + dataset + "/topNode/" + name + "_ifCS=" + numToString<bool>(ifCS) + "_k=" + numToString<int>(k) + "_alpha=" + numToString<double>(alpha) + "_N=" + numToString<int>(N) + "_M=" + numToString<int>(M) + "_ifSub=" + numToString<bool>(ifSub) + "_ifDelta=" + numToString<bool>(ifDelta);

	clock_t begin = clock();
	SparseMatrix* tempInputGraph = inputNetwork(inputAddress);
	tempInputGraph->rowNormalize();
	SparseMatrix* inputGraph = (*tempInputGraph) * alpha;

	delete tempInputGraph;

	if(ifDelta) cout << "delta input!" <<endl;

	if(!ifCS){
		vector<pair<int, double>>* result = katzGreedy(inputGraph, k, N, M, ifDelta, ifSub);
		saveToLocal(result, saveAddress);
	}
	else{
		vector<pair<int, double>>* result = katzGreedyWithPriorSet(inputGraph, k, N, M, ifDelta, ifSub);
		saveToLocal(result, saveAddress);
	}
	clock_t end = clock();
	clock_t runningTime = (double)(end-begin)/CLOCKS_PER_SEC;
	string runningTimeAddress = "../../Output/" + dataset+ "/runningTime.txt";
    ofstream outFile;
    outFile.open(runningTimeAddress.c_str(), ios::app);
    outFile << endl;
	string filename = name + "_ifCS=" + numToString<bool>(ifCS) + "_k=" + numToString<int>(k) + "_alpha=" + numToString<double>(alpha) + "_N=" + numToString<int>(N) + "_M=" + numToString<int>(M) + "_ifSub=" + numToString<bool>(ifSub) + "_ifDelta=" + numToString<bool>(ifDelta);
	outFile << filename << " " << runningTime << endl;
    outFile.close();
}



//-------------Compute the objective function \sum_{i,j}(\sum_{l=1}^k A^l)_{ij}--------------
double matrixSum(SparseMatrix* matrix, int k, int omitIndex, double currentMinSum, Edge* deleteEdge = NULL){
	double result = 0;
	int size = matrix->getSize();
	double* sum0 = new double[size];
	double* sum = new double[size];
	double gamma;
	double* powGamma = new double[k];
	
	for(int i=0; i<size; i++){
		sum0[i] = 1;
		sum[i] = 0;
	}

	//compute the gamma
	set<Edge>::iterator it;

	if(omitIndex != -1)sum0[omitIndex] = 0;

	for(int i=0; i<k; i++){
		//cout << "compute the power " << i <<endl;
		for(int j=0; j<size; j++){
			for(it = matrix->getIterFRowBegin(j); it != matrix->getIterFRowEnd(j); it++){
				sum[it->dst] += it->value * sum0[it->src];
			}
		}
		if(omitIndex != -1)sum[omitIndex] = 0;
		if(deleteEdge != NULL)sum[deleteEdge->dst] -= deleteEdge->value * sum0[deleteEdge->src];
		//Test if it can be pruned.	
		double tempSum = 0;
		for(int j=0; j<size; j++){
			tempSum += sum[j];
			sum0[j] = sum[j];
			sum[j] = 0;
		}
		//End test
		result = result + tempSum;
	}

	delete[] sum0;
	delete[] sum;
	delete[] powGamma;
	return result;
}


vector<vector<double>*>* matrixSumAllValue(SparseMatrix* matrix, int k){
	vector<vector<double>*>* result = new vector<vector<double>*>;
	int size = matrix->getSize();
	double* sum0 = new double[size];
	double* sum = new double[size];
	
	for(int i=0; i<size; i++){
		sum0[i] = 1;
		sum[i] = 0;
	}

	set<Edge>::iterator it;

	for(int i=0; i<k; i++){
		result->push_back(new vector<double>);
		//cout << "compute the power " << i <<endl;
		for(int j=0; j<size; j++){
			for(it = matrix->getIterFRowBegin(j); it != matrix->getIterFRowEnd(j); it++){
				sum[it->dst] += it->value * sum0[it->src];
			}
		}
		//Test if it can be pruned.	
		double tempSum = 0;
		for(int j=0; j<size; j++){
			tempSum += sum[j];
			result->at(i)->push_back(sum[j]);
			sum0[j] = sum[j];
			sum[j] = 0;
		}
		result->at(i)->push_back(tempSum);
	}

	delete[] sum0;
	delete[] sum;
	return result;
}

vector<vector<double>*>* matrixSumAllValueBegin(SparseMatrix* matrix, int k){
	vector<vector<double>*>* result = new vector<vector<double>*>;
	int size = matrix->getSize();
	double* sum0 = new double[size];
	double* sum = new double[size];
	
	for(int i=0; i<size; i++){
		sum0[i] = 1;
		sum[i] = 0;
	}

	set<Edge>::iterator it;

	for(int i=0; i<k; i++){
		result->push_back(new vector<double>);
		//cout << "compute the power " << i <<endl;
		for(int j=0; j<size; j++){
			for(it = matrix->getIterFRowBegin(j); it != matrix->getIterFRowEnd(j); it++){
				sum[it->src] += it->value * sum0[it->dst];
			}
		}
		//Test if it can be pruned.	
		double tempSum = 0;
		for(int j=0; j<size; j++){
			tempSum += sum[j];
			result->at(i)->push_back(sum[j]);
			sum0[j] = sum[j];
			sum[j] = 0;
		}
		result->at(i)->push_back(tempSum);
	}

	delete[] sum0;
	delete[] sum;
	return result;
}

double deltaMatrixSum2(SparseMatrix* matrix, int k, int omitIndex, vector<vector<double>*>* base){
	double result = 0;
	double baseSum = 0;
	double* delta0 = new double[matrix->getSize()];
	double* delta = new double[matrix->getSize()];
	int* oddActivated = new int[matrix->getSize()];
	int* evenActivated = new int[matrix->getSize()];
	for(int i = 0; i<matrix->getSize(); i++){
		delta[i] = 0;
		delta0[i] = 0;
		oddActivated[i] = -1;
		evenActivated[i] = -1;
	}

	set<Edge>::iterator itEdge;
	set<int>::iterator itInt;

	double tempResult = 0;

	for(itEdge = matrix->getIterFRowBegin(omitIndex); itEdge != matrix->getIterFRowEnd(omitIndex); itEdge ++){
		delta0[itEdge->dst] += itEdge->value;
		tempResult += itEdge->value;
		evenActivated[itEdge->dst] = 0;
	}

	for(itEdge = matrix->getIterFColBegin(omitIndex); itEdge != matrix->getIterFColEnd(omitIndex); itEdge ++){
		delta0[omitIndex] += itEdge->value;
		tempResult += itEdge->value;
	}
	evenActivated[omitIndex] = 0;

	baseSum += base->at(0)->at(matrix->getSize());

	result = result + tempResult;

	for(int i=1; i<k; i++){
		tempResult = 0;
		int* oldActivated;
		int* newActivated;
		if(i%2 == 0){
			oldActivated = oddActivated;
			newActivated = evenActivated;
		}
		else{
			oldActivated = evenActivated;
			newActivated = oddActivated;
		}

		for(int j=0; j<matrix->getSize(); j++){
			if(oldActivated[j] < i-1) continue;
			for(itEdge = matrix->getIterFRowBegin(j); itEdge != matrix->getIterFRowEnd(j); itEdge++){
				newActivated[itEdge->dst] = i;
				delta[itEdge->dst] += delta0[itEdge->src]*itEdge->value;
			}
		}
		for(itEdge = matrix->getIterFRowBegin(omitIndex); itEdge != matrix->getIterFRowEnd(omitIndex); itEdge ++){
			if(oldActivated[itEdge->src] == i-1) delta[itEdge->dst] += itEdge->value * (base->at(i-1)->at(itEdge->src)-delta0[itEdge->src]);
			else delta[itEdge->dst] += itEdge->value * base->at(i-1)->at(itEdge->src);
			newActivated[itEdge->dst] = i;
		}

		for(itEdge = matrix->getIterFColBegin(omitIndex); itEdge != matrix->getIterFColEnd(omitIndex); itEdge ++){
			if(oldActivated[itEdge->src] == i-1) delta[itEdge->dst] += itEdge->value * (base->at(i-1)->at(itEdge->src)-delta0[itEdge->src]);
			else delta[itEdge->dst] += itEdge->value * base->at(i-1)->at(itEdge->src);
			newActivated[itEdge->dst] = i;
		}

		for(int j=0; j<matrix->getSize(); j++){
			tempResult = tempResult + delta[j];
			delta0[j] = delta[j];
			delta[j] = 0;
		}

	
		result += tempResult;
		baseSum = baseSum + base->at(i)->at(matrix->getSize());
	}

	delete[] delta;
	delete[] delta0;
	delete[] oddActivated;
	delete[] evenActivated;

	result = baseSum - result;
	return result;
}

double deltaMatrixSum(SparseMatrix* matrix, int k, int omitIndex, vector<vector<double>*>* base){
	double result = 0;
	double baseSum = 0;
	double* delta0 = new double[matrix->getSize()]();
	double* delta = new double[matrix->getSize()]();

	vector<int>* oldActivated = new vector<int>;
	vector<int>* newActivated = new vector<int>;
	double* tempDelta;

	set<Edge>::iterator itEdge;
	vector<int>::iterator itVecInt;

	double tempResult = 0;

	for(itEdge = matrix->getIterFRowBegin(omitIndex); itEdge != matrix->getIterFRowEnd(omitIndex); itEdge ++){
		delta0[itEdge->dst] += itEdge->value;
		tempResult += itEdge->value;
		oldActivated->push_back(itEdge->dst);
	}

	for(itEdge = matrix->getIterFColBegin(omitIndex); itEdge != matrix->getIterFColEnd(omitIndex); itEdge ++){
		delta0[omitIndex] += itEdge->value;
		tempResult += itEdge->value;
	}
	oldActivated->push_back(omitIndex);

	baseSum += base->at(0)->at(matrix->getSize());

	result = result + tempResult;


	for(int i=1; i<k; i++){
		for(itVecInt = oldActivated->begin(); itVecInt != oldActivated->end(); itVecInt++){
			for(itEdge = matrix->getIterFRowBegin(*itVecInt); itEdge != matrix->getIterFRowEnd(*itVecInt); itEdge++){
				if(delta[itEdge->dst] == 0)newActivated->push_back(itEdge->dst);
				delta[itEdge->dst] = delta[itEdge->dst] + delta0[*itVecInt] * itEdge->value;
			}
		}

		for(itEdge = matrix->getIterFRowBegin(omitIndex); itEdge != matrix->getIterFRowEnd(omitIndex); itEdge++){
			delta[itEdge->dst] = delta[itEdge->dst] + (base->at(i-1)->at(omitIndex) - delta0[omitIndex]) * itEdge->value;
		}

		if(delta[omitIndex] == 0) newActivated->push_back(omitIndex);
		delta[omitIndex] = base->at(i)->at(omitIndex);

		for(itVecInt = newActivated->begin(); itVecInt != newActivated->end(); itVecInt++){
			result = result + delta[*itVecInt];
		}

		delete delta0;
		delta0 = delta;
		delta = new double[matrix->getSize()]();

		delete oldActivated;
		oldActivated = newActivated;
		newActivated = new vector<int>;

		baseSum = baseSum + base->at(i)->at(matrix->getSize());
	}

	delete[] delta;
	delete[] delta0;
	delete oldActivated;
	delete newActivated;

	result = baseSum - result;
	return result;
}

//---------------Find the top node set---------------------
vector<pair<int, double>>* katzGreedy(SparseMatrix* matrix,int k, int num, int M, bool ifDelta, bool submodule){
	vector<pair<int, double>>* result = new vector<pair<int, double>>;
	double base = matrixSum(matrix,k,-1,matrix->getSize()*k*2);
	cout << "base is " <<base <<endl;

	double iterBase = base;
	double* delta = new double[matrix->getSize()];
	if(submodule){
		for(int i=0; i<matrix->getSize(); i++){
			if(matrix->getIterFRowBegin(i) == matrix->getIterFRowEnd(i) && matrix->getIterFColBegin(i) == matrix->getIterFColEnd(i)){
				delta[i] = 0;
			}
			else{
				delta[i] = matrix->getSize()*k*2;
			}
		}
	}

	for(int i=0; i<num; i++){
		cout<<"begin to select "<<i << " node"<<endl;
		int size = matrix->getSize();
		int optimalId = -1;
		double minValue = iterBase;
		double temp = -1;
		
		vector<vector<double>*>* currentBase;
		if(ifDelta) currentBase = matrixSumAllValue(matrix, k);

		// make the minvalue as small as possible by using fewer time.
		int priorNum = M;
		set<pair<double, int>>* priorSet = new set<pair<double, int>>;
		set<int>* priorNodeSet = new set<int>;


		for(int j=0; j<size; j++){
			set<Edge>::iterator it;
			double tempSum = 0;
			for(it = matrix->getIterFColBegin(j); it != matrix->getIterFColEnd(j); it++){
				tempSum = tempSum + it->value;
			}
			for(it = matrix->getIterFRowBegin(j); it != matrix->getIterFRowEnd(j); it++){
				tempSum = tempSum + it->value;
			}
			if(priorSet->size() < priorNum) priorSet->insert(pair<double, int>(tempSum, j));
			else if(tempSum > (priorSet->begin())->first){
				priorSet->erase(*(priorSet->begin()));
				priorSet->insert(pair<double, int>(tempSum, j));
			}
		}
		set<pair<double, int>>::iterator itPair;

		for(itPair = priorSet->begin(); itPair != priorSet->end(); itPair++){
			priorNodeSet->insert(itPair->second);
			if(submodule && delta[itPair->second] <= iterBase - minValue) continue;

			if(ifDelta) temp = deltaMatrixSum(matrix,k,itPair->second, currentBase);
			else temp = matrixSum(matrix,k,itPair->second,minValue);

			if(submodule) delta[itPair->second] = iterBase - temp;
			if(temp < minValue){
				minValue = temp;
				optimalId = itPair->second;
			}
		}

		if(ifDelta) cout << "delta!" <<endl;
		if(submodule) cout << "sub!" << endl;

		for(int j=0; j<size; j++){
			if(priorNodeSet->find(j) != priorNodeSet->end()) continue;
			if(submodule && delta[j] <= iterBase - minValue) continue;

			if(ifDelta) temp = deltaMatrixSum(matrix,k,j, currentBase);
			else temp = matrixSum(matrix,k,j,minValue);

			if(submodule)delta[j] = iterBase - temp;
			if(temp < minValue){
				minValue = temp;
				optimalId = j;
			}
		}

		//cout<<endl;

		//update the network

		if(optimalId == -1) break;
		set<Edge>::iterator it;
		for(it = matrix->getIterFRowBegin(optimalId); it!= matrix->getIterFRowEnd(optimalId); ){
			Edge temp(it->src, it->dst, 0);
			it++;
			matrix->erase(temp);
		}
		for(it = matrix->getIterFColBegin(optimalId); it!= matrix->getIterFColEnd(optimalId); ){
			Edge temp(it->src, it->dst, 0);
			it++;
			matrix->erase(temp);
		}
		cout << "The " << i <<"'s node is "<< optimalId <<" and the remaining value is "<<minValue/base;
		cout << endl;


		iterBase = minValue;
		result->push_back(pair<int, double>(optimalId, minValue/base));
		//delete currentBase;
		delete priorSet;
		delete priorNodeSet;
		if(ifDelta){
			for(int l=0; l<currentBase->size(); l++){
				delete currentBase->at(l);
			}
		}
	}

	delete[] delta;
	return result;
}

//---------------Find the top node set with the prior set in each iteration---------------------
vector<pair<int, double>>* katzGreedyWithPriorSet(SparseMatrix* matrix,int k, int num, int priorNum,bool ifDelta, bool submodule){
	vector<pair<int, double>>* result = new vector<pair<int, double>>;
	double base = matrixSum(matrix,k,-1, matrix->getSize()*k*2);

	cout << "base is " << base <<endl;
	double iterBase = base;
	double* delta = new double[matrix->getSize()];
	for(int i=0; i<matrix->getSize(); i++){
		if(matrix->getIterFRowBegin(i) == matrix->getIterFRowEnd(i) && matrix->getIterFColBegin(i) == matrix->getIterFColEnd(i)){
			delta[i] = 0;
		}
		else{
			delta[i] = matrix->getSize()*k*2;
		}
	}

	for(int i=0; i<num; i++){
		cout<<"begin to select "<<i << " node" << endl;
		int size = matrix->getSize();
		int optimalId = -1;
		double minValue = iterBase;
		double temp = -1;

		int optimalRankInPriorSet;

		set<pair<double, int>>* priorSet = new set<pair<double, int>>;
		for(int j=0; j<size; j++){
			set<Edge>::iterator it;
			double tempSum = 0;
			for(it = matrix->getIterFColBegin(j); it != matrix->getIterFColEnd(j); it++){
				tempSum = tempSum + it->value;
			}
			for(it = matrix->getIterFRowBegin(j); it != matrix->getIterFRowEnd(j); it++){
				tempSum = tempSum + it->value;
			}
			if(priorSet->size() < priorNum) priorSet->insert(pair<double, int>(tempSum, j));
			else if(tempSum > (priorSet->begin())->first){
				priorSet->erase(*(priorSet->begin()));
				priorSet->insert(pair<double, int>(tempSum, j));
			}
		}
		//cout << "finish the selection of prior set " << endl;

		
		vector<vector<double>*>* currentBase; 
		if(ifDelta) currentBase = matrixSumAllValue(matrix, k);
		//cout << "finishe the computation of currentBase" << endl;

		set<pair<double, int>>::iterator itPair;
		for(itPair = priorSet->begin(); itPair != priorSet->end(); itPair++){
			if(submodule && delta[itPair->second] <= iterBase - minValue){
				continue;
			}
			if(ifDelta) temp = deltaMatrixSum(matrix, k, itPair->second, currentBase);
			else temp = matrixSum(matrix, k, itPair->second, minValue);
			if(submodule) delta[itPair->second] = iterBase - temp;
			if(temp < minValue){
				minValue = temp;
				optimalId = itPair->second;
			}
		}

		if(optimalId == -1) break;
		//update the iterBase
		iterBase = minValue;

		set<Edge>::iterator it;
		for(it = matrix->getIterFRowBegin(optimalId); it!= matrix->getIterFRowEnd(optimalId); ){
			Edge temp(it->src, it->dst, 0);
			it++;
			matrix->erase(temp);
		}
		for(it = matrix->getIterFColBegin(optimalId); it!= matrix->getIterFColEnd(optimalId); ){
			Edge temp(it->src, it->dst, 0);
			it++;
			matrix->erase(temp);
		}
		//cout << priorSet->size()<<endl;
		cout << "The " << i <<"'s node is "<< optimalId <<" and the remaining value is "<<minValue/base << endl;


		result->push_back(pair<int, double>(optimalId, minValue/base));
		delete priorSet;
		if(ifDelta){
			for(int l=0; l<currentBase->size(); l++){
				delete currentBase->at(l);
			}
		}
	}

	delete delta;
	return result;
}



//-----------second block----------------

void labelNetworkEdge(string address, SparseMatrix* edgeSet, SparseMatrix* matrix){
        ofstream file;
        file.open(address, ios::out);
        file<<"edgeNum: "<< matrix->getNonZeroNum() <<"\n";
        set<Edge>::iterator it, it1, it2;
        for(int i=0; i<matrix->getSize(); i++){
                it1 = matrix->getIterFRowBegin(i);
                it2 = matrix->getIterFRowEnd(i);
                for(it = it1; it!= it2; it++){
			if(edgeSet->ifExist(*it)){
				file << it->src << " " << it->dst << " " << it->value << " top" << endl;
			}
			else{
				file << it->src << " " << it->dst << " " << it->value << " not" << endl;
			}
		}
        }
	file.close();
}

//-----------give labels to each node and save----------------
void labelNetworkNode(string address, set<int>* nodeSet, set<int>* rawNodeSet, vector<string>* dictionary){
	ofstream file;
	file.open(address, ios::out);
	file << "nodeNum: " << dictionary->size() << endl;
	for(int i=0; i<dictionary->size(); i++){
		if(rawNodeSet->find(i) != rawNodeSet->end()){
			if(nodeSet->find(i) != nodeSet->end()){
				file << "top " << dictionary->at(i) << endl;
			}
			else{
				file << "not " << dictionary->at(i) << endl;
			}
		}
	}
	file.close();
}


//-----------Third block----------
//
//
//-----------transformation between nodes and edges----------


void nodesToEdges(set<int>* nodeSet, SparseMatrix* matrix, string saveAddress){
	SparseMatrix* result = new SparseMatrix(matrix->getSize());

	set<Edge>::iterator it;
	for(int i=0; i<matrix->getSize(); i++){
		for(it = matrix->getIterFRowBegin(i); it != matrix->getIterFRowEnd(i); it++){
			if( (nodeSet->find(it->src) != nodeSet->end()) && nodeSet->find(it->dst) != nodeSet->end()){
				result->insert(Edge(it->src, it->dst, it->value));
			}
		}
	}

	result->saveToLocal(saveAddress);

	delete result;
}

void edgesToNodes(SparseMatrix* edgeSet, string saveAddress){
	set<int>* result = new set<int>;

	set<Edge>::iterator it;
	for(int i=0; i<edgeSet->getSize(); i++){
		for(it = edgeSet->getIterFRowBegin(i); it != edgeSet->getIterFRowEnd(i); it++){
			result->insert(it->src);
			result->insert(it->dst);
		}
	}

        ofstream file;
        file.open(saveAddress, ios::out);

	file << result->size() << endl;
	set<int>::iterator itInt;
	for(itInt = result->begin(); itInt != result->end(); itInt ++){
		file << *itInt << endl;
	}

	file.close();
}




//-----------sixth block--------------
//
//
//-----------read in a network--------------
SparseMatrix* inputNetwork(string filename, bool ifSelf){
	SparseMatrix* result;

    	int nodeNum;
	int edgeNum;

	ifstream infile(filename.c_str());
	cout<<"Begin to read the file: "<<filename<<endl;

	string s,s2;
	infile >> s >> nodeNum;
	getline(infile,s);
	cout<<"NodeNum is "<<nodeNum<<endl;
	result = new SparseMatrix(nodeNum);


    	getline(infile,s);
    	istringstream line (s);
    	line >> s2 >> edgeNum;
    	cout << "EdgeNum is "<< edgeNum<<endl;

    	int tempNode1, tempNode2;
	double temp;

    	for(int i=0; i<edgeNum; i++){
        	getline(infile,s);
        	istringstream line (s);
        	line >> tempNode1 >> tempNode2 >> temp;
/*		cout << tempNode1 << " " << tempNode2 <<endl;
		if(tempNode1 > nodeNum || tempNode2 > nodeNum){
			cout << "out of index. wrong with indexNetwork.txt" <<endl;
			exit(0);
		}
*/		if(tempNode1 == tempNode2 && (!ifSelf)) continue;
        	result->insert(Edge(tempNode1, tempNode2, temp));
    	}	


        cout << "Finish the reading! Real edgeNum is" << result->getNonZeroNum()<<endl;
    	return result;
}

//-----------read in a dictionary--------------
vector<string>* inputDictionary(string dictionary){
        ifstream infile(dictionary.c_str());
        string s, s1;
        int num, temp;
        getline(infile, s);
        istringstream line (s);
        line >> num;
        //cout << num;
        vector<string>* indexToName = new vector<string>;
        for(int i=0; i<num; i++){
                getline(infile, s);
                istringstream line (s);
                line >> temp;
                //cout << temp << endl;
		if(temp != i){
			cout << "format error with indexToName.txt" << endl;
			exit(0);
		}
                indexToName->push_back(s);
        }
//        cout << indexToName->size();

        infile.close();

	return indexToName;
}


//-----------read in a node set---------------
set<int>* inputNodeSet(string address){
        set<int>* result = new set<int>;
        ifstream infile(address.c_str());
        int nodeNum;
        infile >> nodeNum;
        cout << nodeNum <<endl;
        string s;
        getline(infile,s);
        for(int i=0; i<nodeNum; i++){
                string s;
                getline(infile, s);
                istringstream line(s);
                int nodeId;
                line >> nodeId;
                result->insert(nodeId);
        }
        return result;
}

//-------------save the node list-------------
void saveToLocal(vector<pair<int, double>>* result, string address){
	ofstream file;
	cout << "begin to save " << address << endl;
	file.open(address, ios::out);
	int size = result->size();
	file << size << endl;
	for(int i=0; i<size; i++){
		file << result->at(i).first << " " << result->at(i).second << endl;
	}
	cout << "finish saving " << address <<endl;
	file.close();
}

void saveToLocal(vector<int>* result, string address){
	ofstream file;
	cout << "begin to save " << address << endl;
	file.open(address, ios::out);
	int size = result->size();
	file << size << endl;
	for(int i=0; i<size; i++){
		file << result->at(i) << endl;
	}
	cout << "finish saving " << address <<endl;
	file.close();
}
void saveToLocal(pair<vector<int>*, vector<pair<int, double>>*>* result, string address){
	ofstream file;
	cout << "begin to save " << address << endl;
	file.open(address, ios::out);
	int size = result->second->size();
	file << size << endl;
	for(int i=0; i<size; i++){
		file << result->second->at(i).first << " " << result->second->at(i).second << " " << result->first->at(i) << endl;
	}
	cout << "finish saving " << address <<endl;
	file.close();
}




template <class Type>
string numToString(Type value) {
    stringstream ss;
    ss << value;
    return ss.str();
}



vector<pair<int, double>>* pageRank(SparseMatrix* matrix, int resultNum, double alpha, double threshold){
	vector<pair<int, double>>* result = new vector<pair<int, double>>;
	double* oldRank = new double[matrix->getSize()];
	double* newRank = new double[matrix->getSize()];

	for(int i=0; i<matrix->getSize(); i++){
		oldRank[i] = 1.0/matrix->getSize();
		newRank[i] = 0;
	}

	bool ifConverge = false;

	int index = 0;
	while(!ifConverge){
		cout << index << " iteration" << endl;
		index ++;
		set<Edge>::iterator itEdge;
		ifConverge = true;
		for(int i=0; i<matrix->getSize(); i++){
			for(itEdge = matrix->getIterFColBegin(i); itEdge != matrix->getIterFColEnd(i); itEdge++){
				newRank[i] += oldRank[itEdge->src] * itEdge->value;
			}
			newRank[i] = newRank[i] * alpha;
			newRank[i] = newRank[i] + double(1-alpha)/matrix->getSize();
		}
		for(int i=0; i<matrix->getSize(); i++){
			if(newRank[i] - oldRank[i] > threshold || oldRank[i] - newRank[i] > threshold){
				ifConverge = false;
			}
			oldRank[i] = newRank[i];
			newRank[i] = 0;
		}
	}

	set<pair<double, int>>* topSet = new set<pair<double, int>>;
	for(int i=0; i<matrix->getSize(); i++){
		if(topSet->size() < resultNum){
			pair<double, int>* temp = new pair<double, int>(oldRank[i], i);
			topSet->insert(*temp);
		}
		else if(oldRank[i] > (topSet->begin())->first){
			topSet->erase(pair<double, int>(topSet->begin()->first, topSet->begin()->second));
			pair<double, int>* temp = new pair<double, int>(oldRank[i], i);
			topSet->insert(*temp);
		}
	}

	set<pair<double, int>>::reverse_iterator itPair;
	for(itPair = topSet->rbegin(); itPair != topSet->rend(); itPair ++){
		cout << itPair->second << " " << itPair->first <<endl;
		result->push_back(pair<int, double>(itPair->second, itPair->first));
	}

	
	delete topSet;
	delete[] oldRank;
	delete[] newRank;
	return result;
}

vector<int>* diversifiedRank(SparseMatrix* matrix, int resultNum, double alpha, double threshold, vector<double>* p){
	vector<int>* result = new vector<int>;
	
	if(p == NULL){
		cout << "success" << endl;
		p = new vector<double>;
		for(int i=0; i<matrix->getSize(); i++){
			p->push_back(1.0/matrix->getSize());
		}
	}

	double* oldRank = new double[matrix->getSize()];
	double* newRank = new double[matrix->getSize()];

	for(int i=0; i<matrix->getSize(); i++){
		oldRank[i] = 1.0/matrix->getSize();
		newRank[i] = 0;
	}

	bool ifConverge = false;

	int index= 0;
	while(!ifConverge){
		cout << index << " iteration" << endl;
		index++;
		set<Edge>::iterator itEdge;
		ifConverge = true;
		for(int i=0; i<matrix->getSize(); i++){
			for(itEdge = matrix->getIterFColBegin(i); itEdge != matrix->getIterFColEnd(i); itEdge++){
				newRank[i] += oldRank[itEdge->src] * itEdge->value;
			}
			newRank[i] = newRank[i] * alpha;
			newRank[i] = newRank[i] + double(1-alpha) * p->at(i);
		}
		for(int i=0; i<matrix->getSize(); i++){
			if(newRank[i] - oldRank[i] > threshold || oldRank[i] - newRank[i] > threshold){
				ifConverge = false;
			}
			oldRank[i] = newRank[i];
			newRank[i] = 0;
		}
	}

	double* s = new double[matrix->getSize()];
	double* u = new double[matrix->getSize()];
	double* v = new double[matrix->getSize()];

	for(int i=0; i<matrix->getSize(); i++){
		s[i] = (2 - (1 - alpha) * p->at(i)) * oldRank[i];
		u[i] = 0;
		v[i] = 0;
	}

	set<int>* resultSet = new set<int>;	
	for(int iter=0; iter<resultNum; iter++){
		cout << "selecting "<<iter << "'s node" << endl;
		//compute the score vector s = s - u\cdot r - v;
		for(int j=0; j<matrix->getSize(); j++){
			s[j] = s[j] - u[j] * oldRank[j] - v[j];
		}
		//Find i = argmax_j s_j(j = 1,2,cdots, n; j\in S)
		int optimalId;
		double maxValue = -matrix->getSize();
		for(int j = 0; j<matrix->getSize(); j++){
			if(s[j] > maxValue){
				if(resultSet->find(j) == resultSet->end()){
					optimalId = j;
					maxValue = s[j];
				}
			}
		}
		//Add node i into S
		result->push_back(optimalId);
		resultSet->insert(optimalId);
		//update u := u + c A(:,i) + (1-c) p(i)1_{n\times 1};
		for(int j = 0; j<matrix->getSize(); j++){
			u[j] = u[j] + (1 - alpha) * p->at(optimalId);
		}
		set<Edge>::iterator itEdge;
		for(itEdge = matrix->getIterFColBegin(optimalId); itEdge != matrix->getIterFColEnd(optimalId); itEdge++){
			u[itEdge->src] = u[itEdge->src] + itEdge->value * alpha;
		}
		//update v := v + c A'(:,i) r(i) + (1 - c) r(i) p;
		for(int j=0; j<matrix->getSize(); j++){
			v[j] = v[j] + oldRank[optimalId] * p->at(j) * (1 - alpha);
		}
		for(itEdge = matrix->getIterFRowBegin(optimalId); itEdge != matrix->getIterFRowEnd(optimalId); itEdge++){
			v[itEdge->dst] = v[itEdge->dst] + alpha * itEdge->value * oldRank[optimalId];
		}
	}

	delete[] s;
	delete[] u;
	delete[] v;
	delete[] oldRank;
	delete[] newRank;
	delete resultSet;
	return result;
}

void baseline_connectivity(){
	cout << "dataset name, name, N" <<endl;
	string dataset, name;
	int N;
	cin >> dataset >> name >> N;
	string inputAddress = "../../Dataset/" + dataset + "/indexNetwork.txt";
	string saveAddress = "../../Output/" + dataset + "/baseline/connectivity-2.txt";
	SparseMatrix* tempInputGraph = inputNetwork(inputAddress);
	vector<pair<int, double>>* result = connecivity(tempInputGraph, N);
	saveToLocal(result, saveAddress);
}

vector<pair<int, double>>* connecivity(SparseMatrix* matrix, int N){
	vector<pair<int, double>>* result = new vector<pair<int, double>>;

	vector<int>* maxIS = maximalIndepSet(matrix);

	cout << "maximal indepedent set: " << endl;
	for(vector<int>::iterator it = maxIS->begin(); it!= maxIS->end(); it++){
		cout << *it << endl;
	}

	int* connectedId = new int[matrix->getSize()];
	int* connectedNum = new int[matrix->getSize()];
	bool* idIfSelected = new bool[matrix->getSize()];

	for(int i=0; i<matrix->getSize(); i++){
		connectedId[i] = -1;
		connectedNum[i] = 0;
		idIfSelected[i] = false;
	}
	for(int i=0; i<maxIS->size(); i++){
		connectedId[maxIS->at(i)] = i;
		connectedNum[i] = 1;
	}

	//If the size of maximal indepedent set has already been larger than the goal, we are done in this situation.
	if(maxIS->size() > matrix->getSize() - N){
		for(int i=0; i<matrix->getSize(); i++){
			if(connectedId[i] != -1)result->push_back(pair<int, double>(i,0) );
		}
		return result;
	}
	for(int i=0; i<maxIS->size(); i++){
			result->push_back( pair<int, double>(maxIS->at(i), 0) );
	}

	cout << "the second step: " <<endl;

	set<Edge>::iterator itEdge;

	double minValue;
	double testValue = 0;

	for(int i=maxIS->size(); i<matrix->getSize() - N; i++){
		minValue= matrix->getSize() * matrix->getSize();
		int optId;
		for(int j=0; j<matrix->getSize(); j++){
			if(connectedId[j] != -1) continue;
			double tempValue = 0;

			vector<int>* selectedId = new vector<int>;
			for(itEdge = matrix->getIterFRowBegin(j); itEdge != matrix->getIterFRowEnd(j); itEdge++){
				if(!idIfSelected[ connectedId[itEdge->dst] ] && connectedId[itEdge->dst] != -1){
					idIfSelected[ connectedId[itEdge->dst] ] = true;
					selectedId->push_back(connectedId[itEdge->dst]);
					tempValue += connectedNum[connectedId[itEdge->dst]];
				}
			}

			for(itEdge = matrix->getIterFColBegin(j); itEdge != matrix->getIterFColEnd(j); itEdge++){
				if(!idIfSelected[ connectedId[itEdge->src] ] && connectedId[itEdge->src] != -1){
					idIfSelected[ connectedId[itEdge->src] ] = true;
					selectedId->push_back(connectedId[itEdge->src]);
					tempValue += connectedNum[connectedId[itEdge->src]];
				}
			}

			tempValue = (tempValue + 1) * (tempValue) / 2;

			for(int l=0; l<selectedId->size(); l++){
				tempValue = tempValue - connectedNum[selectedId->at(l)] * (connectedNum[selectedId->at(l)] - 1)/2;
				idIfSelected[selectedId->at(l)] = false;
			}

			if(tempValue < minValue){
				minValue = tempValue;
				optId = j;
			}

			delete selectedId;
		}


		cout << optId << " " << minValue << " ";

		testValue += minValue;


		vector<int>* selectedId = new vector<int>;
		for(itEdge = matrix->getIterFRowBegin(optId); itEdge != matrix->getIterFRowEnd(optId); itEdge++){
			if(!idIfSelected[ connectedId[itEdge->dst] ] && connectedId[itEdge->dst] != -1){
				idIfSelected[ connectedId[itEdge->dst] ] = true;
				selectedId->push_back(connectedId[itEdge->dst]);
			}
		}

		for(itEdge = matrix->getIterFColBegin(optId); itEdge != matrix->getIterFColEnd(optId); itEdge++){
			if(!idIfSelected[ connectedId[itEdge->src] ] && connectedId[itEdge->src] != -1){
				idIfSelected[ connectedId[itEdge->src] ] = true;
				selectedId->push_back(connectedId[itEdge->src]);
			}
		}


//		cout << selectedId->at(0) << endl;

		connectedId[optId] = selectedId->at(0);
		connectedNum[selectedId->at(0)] += 1;


		for(int l=1; l<selectedId->size(); l++){
			connectedNum[selectedId->at(0)] += connectedNum[selectedId->at(l)];
			connectedNum[selectedId->at(l)] = 0;
		}


		double optValue = 0;
		for(int l=0; l<matrix->getSize(); l++){
			optValue = optValue + double (connectedNum[l] * (connectedNum[l] - 1)) / 2;
		}

		if(optValue != testValue){
			cout << "error!" <<endl;
			cout << optValue << " " << testValue << endl;
			exit(0);
		}

		result->push_back(pair<int, double>(optId,optValue * 2 /(matrix->getSize() * (matrix->getSize()- 1)  )  ) );

		//cout << optValue * 2 /(i * (i + 1)  ) << endl;

		

		for(int l=0; l<matrix->getSize(); l++){
			if(idIfSelected[connectedId[l]]){
				connectedId[l] = selectedId->at(0);
			}
		}


		for(int l=0; l<selectedId->size(); l++){
			idIfSelected[selectedId->at(l)] = false;
		}

		delete selectedId;
	}



	if(result->size() != matrix->getSize() - N){
		cout << "There must be a bug here !" << endl;
		exit(0);
	}

	//NOTICE: Here we assume that the original graph is connected;
	delete[] connectedId;
	delete[] connectedNum;
	delete[] idIfSelected;
	return result;
}

vector<int>* maximalIndepSet(SparseMatrix* matrix){
	vector<int>* result = new vector<int>;
	bool* ifNeb = new bool[matrix->getSize()];
	for(int i=0; i<matrix->getSize(); i++){
		ifNeb[i] = false;
	}

	int pointer = 0;

	set<Edge>::iterator itEdge;

	while(pointer != matrix->getSize()){
		result->push_back(pointer);
		ifNeb[pointer] = true;
		for(itEdge = matrix->getIterFRowBegin(pointer); itEdge != matrix->getIterFRowEnd(pointer); itEdge++){
			ifNeb[itEdge->dst] = true;
		}
		for(itEdge = matrix->getIterFColBegin(pointer); itEdge != matrix->getIterFColEnd(pointer); itEdge++){
			ifNeb[itEdge->src] = true;
		}
		while(ifNeb[pointer]&&pointer < matrix->getSize()){
			pointer++;
		}
	}
	delete[] ifNeb;
	return result;
}


//-----------find top node set with a faster approximate method---------
vector<pair<int, double>>* katzGreedyWithAppoximateCriteria(SparseMatrix* matrix, int k, int size){
	vector<pair<int, double>>* result = new vector<pair<int, double>>;

	//*********************Accuracy**************
	double base = matrixSum(matrix,k,-1, matrix->getSize()*k*2);
//	double base = 1;
	cout << "base is "<< base << endl;



	for(int i=0; i<size; i++){
		cout << "begin to select " << i <<"'s node" << endl;
		vector<vector<double>*>* baseEndWith = matrixSumAllValue(matrix, k);
		vector<vector<double>*>* baseBeginWith = matrixSumAllValueBegin(matrix, k);

		double* delta = new double[matrix->getSize()];
		for(int j=0; j<matrix->getSize(); j++){
			delta[j] = 0;
		}
//		cout << "finish computing the current base" << endl;


		for(int j=0; j<matrix->getSize(); j++){
			baseBeginWith->at(0)->at(j) += 1;
		}
		for(int l=1; l<k; l++){
			for(int j=0; j<matrix->getSize(); j++){
				baseBeginWith->at(l)->at(j) += baseBeginWith->at(l-1)->at(j);
//				baseBeginWith->at(l)->at(j) += 1;
			}
		} 

		set<Edge>::iterator itEdge, itEdge2;

		for(int j=0; j<matrix->getSize(); j++){
			for(itEdge = matrix->getIterFRowBegin(j); itEdge != matrix->getIterFRowEnd(j); itEdge++){
				for(itEdge2 = matrix->getIterFColBegin(j); itEdge2 != matrix->getIterFColEnd(j); itEdge2++){
					double temp = itEdge->value * itEdge2->value;
					for(int l=0; l<k-3; l++){
						delta[j] += baseEndWith->at(l)->at(itEdge2->src) * temp * baseBeginWith->at(k-4-l)->at(itEdge->dst);
					}
					delta[j] += baseEndWith->at(k-3)->at(itEdge2->src) * temp;
					delta[j] += temp * baseBeginWith->at(k-3)->at(itEdge->dst);
				}
				delta[j] += itEdge->value * baseBeginWith->at(k-2)->at(itEdge->dst);
			}
			for(itEdge2 = matrix->getIterFColBegin(j); itEdge2 != matrix->getIterFColEnd(j); itEdge2++){
				delta[j] += itEdge2->value * baseEndWith->at(k-2)->at(itEdge2->src);
			}
		}

		double maxValue = -1;
		int optId;
		for(int j=0; j<matrix->getSize(); j++){
			if(delta[j] > maxValue){
				maxValue = delta[j];
				optId = j;
			}
		}

		double remaining = deltaMatrixSum(matrix,k,optId, baseEndWith);	

		set<Edge>::iterator it;
		for(it = matrix->getIterFRowBegin(optId); it!= matrix->getIterFRowEnd(optId); ){
			Edge temp(it->src, it->dst, 0);
			it++;
			matrix->erase(temp);
		}
		for(it = matrix->getIterFColBegin(optId); it!= matrix->getIterFColEnd(optId); ){
			Edge temp(it->src, it->dst, 0);
			it++;
			matrix->erase(temp);
		}
		cout << "The " << i << "'s node is " << optId << ", and the remaining value is "<<remaining  <<endl;


		result->push_back(pair<int, double>(optId, (double(remaining))/base));

		for(int l=0; l<k; l++){
			delete baseEndWith->at(l);
			delete baseBeginWith->at(l);
		}
		delete[] delta;

	}

	return result;
}



//------------selfTest--------------
void computeRemaining(){
	string nodeAddress = "../../Output/citation_learning/topNode/baseline/diversedRank_week12_size=100_alpha=0.9_threshold=1e-15";
	string networkAddress = "../../Dataset/citation_learning/indexNetwork.txt";
	int k=10;
	double alpha = 0.9;

	set<int>* nodeSet = inputNodeSet(nodeAddress);
	SparseMatrix* temp = inputNetwork(networkAddress);
	temp->rowNormalize();
	SparseMatrix* matrix = (*temp)*alpha;

	double base = matrixSum(matrix,k,-1, matrix->getSize()*k*2);
	cout << "base is " << base << endl;


	set<int>::iterator itInt;

	for(itInt = nodeSet->begin(); itInt != nodeSet->end(); itInt++){
		int optId = *itInt;
		set<Edge>::iterator it;
		double remaining1 = matrixSum(matrix,k,optId, matrix->getSize()*k*2);
		for(it = matrix->getIterFRowBegin(optId); it!= matrix->getIterFRowEnd(optId); ){
			Edge temp(it->src, it->dst, 0);
			it++;
			matrix->erase(temp);
		}
		for(it = matrix->getIterFColBegin(optId); it!= matrix->getIterFColEnd(optId); ){
			Edge temp(it->src, it->dst, 0);
			it++;
			matrix->erase(temp);
		}
		double remaining = matrixSum(matrix,k,-1, matrix->getSize()*k*2);
        	cout << "Afer removing " << optId << ", remaining value is " << remaining/base << " and "<<  remaining/base <<endl;
	}
}
