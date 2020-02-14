#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>

using namespace std;

void readInputData(const char* fileName, vector<double> *tim, vector<double> *tem, vector<double> *pre){
	ifstream input;
	input.open(fileName);
	vector<string> lines;
	string line;
	while(getline(input, line)){//read in lines from input data file and store in vector of strings
		lines.push_back(line);
	}
	for(int i=0; i<lines.size(); i++){//parse/split 'lines' strings by " " here into substrings
		istringstream iss(lines[i]);
		for(int j=0; j<3; j++){
			char* pEnd;
			string sub;
			double doub;
			iss >> sub;
			const char *charp = sub.c_str();
			doub = strtod(charp, &pEnd);
			if(j == 0)
				tim->push_back(doub);
			else if(j == 1)
				tem->push_back(doub);
			else if(j == 2)
				pre->push_back(doub);
		}
	}

}

void printTimes(vector<double> *tim){
	cout << "Times:" << endl;
	for(int i=0; i<tim->size(); i++){
		cout << (*tim)[i] << endl;
	}
	cout << endl;
}

void printTemps(vector<double> *tem){
	cout << "Temperatures:" << endl;
	for(int i=0; i<tem->size(); i++){
		cout << (*tem)[i] << endl;
	}
	cout << endl;
}

void printPress(vector<double> *pre){
	cout << "Pressures:" << endl;
	for(int i=0; i<pre->size(); i++){
		cout << (*pre)[i] << endl;
	}
	cout << endl;
}

void findBestFitLine(vector<double> *xData, vector<double> *yData, double *resultSlope, double *resultIntercept){
	double sumx, sumy, sumop, sumxsq;
	sumx = 0;
	sumy = 0;
	sumop = 0;
	sumxsq = 0;
	int n = xData->size();
	for(int i=0; i<xData->size(); i++){
		sumx += (*xData)[i];
		sumy += (*yData)[i];
		sumop += ((*xData)[i]*(*yData)[i]);
		sumxsq += ((*xData)[i]*(*xData)[i]);
	}

	double slope = (n*sumop - sumx*sumy)/(n*sumxsq - (sumx*sumx));
	double intercept = (sumy - (slope*sumx))/n;

	*resultSlope = slope;
	*resultIntercept = intercept;
}

void removeEntry(vector<double> *a, vector<double> *b, vector<double> *c, int index){
	a->erase(a->begin() + index);
	b->erase(b->begin() + index);
	c->erase(c->begin() + index);
}

double findStdDevofRegLine(vector<double> *xData, vector<double> *yData, double slope, double intercept){
	double n = xData->size()*1.0;
	double num = 0;
	double y, x, e;
	for(int k=0; k<n; k++){
		y = (*yData)[k];
		x = (*xData)[k];
		e = (slope*x + intercept);
//		cout << "Expected y: " << e << endl;
//		cout << "Actual y: " << y << endl;
		num += fabs(y - e);
//		cout << "numerator is: " << num << endl;
//		cout << "n is: " << n << endl;
	}
	double result = pow((num/n), 0.5);
//	cout << "result for std dev is: " << result << endl;
	return result;
}

void prune(vector<double> *xData, vector<double> *yData, double slope, double intercept, int factor, double std, vector<double> *other){
	int deleted = 0;
	double yact, yexp, x;
	int i = 0;
	while(i<(xData->size()-deleted)){
		yact = (*yData)[i];
		x = (*xData)[i];
		yexp = x*slope+intercept;
		cout << "yact = " << yact << endl;
		cout << "yexp = " << yexp << endl;
		cout << "std = " << std << endl << endl;
		if(fabs(yact-yexp)>(factor*std)){
			cout << "Choosing to remove." << endl;
			removeEntry(xData, yData, other, i);
			deleted++;
		}
		else{
			i++;
			cout << "Choosing to keep." << endl;
		}
	}

}

void writeOutputFile(vector<double> *tim, vector<double> *tem, vector<double> *pre){
	ofstream output;
	output.open("output.txt");
	for(int t=0; t<tim->size(); t++){
		output << (*tim)[t] << " " << (*tem)[t] << " " << (*pre)[t] << "\n";
	}
	output.close();
	cout << "Output file successfully written to." << endl;
}

void genRandomInputFile(){
	ofstream output;
	output.open("inputRand.txt");
	cout << "You chose to generate a random input file." << endl;
	output.close();

}

int main(int argc, char* argv[]){

	cout << endl;

	vector<double> *timesP = new vector<double>;
	vector<double> *tempsP = new vector<double>;
	vector<double> *pressP = new vector<double>;

	readInputData(argv[1], timesP, tempsP, pressP);
	printTimes(timesP);
	printTemps(tempsP);
	printPress(pressP);

	double tempSlope;
	double tempInt;
	double presSlope;
	double presInt;

	double *tempSP = &tempSlope;
	double *tempIP = &tempInt;
	double *presSP = &presSlope;
	double *presIP = &presInt;

	findBestFitLine(timesP, tempsP, tempSP, tempIP);
	findBestFitLine(timesP, pressP, presSP, presIP);

	double tempStdDev = findStdDevofRegLine(timesP, tempsP, tempSlope, tempInt);
	double presStdDev = findStdDevofRegLine(timesP, pressP, presSlope, presInt);

	cout << endl << "Best fit line for temp vs time:" << endl;
	cout << "y = " << tempSlope << "x + " << tempInt << endl;
	cout << "Standard deviation for the line is: " << tempStdDev << endl << endl;


	cout << "Best fit line for pressure vs time:" << endl;
	cout << "y = " << presSlope << "x + " << presInt << endl;
	cout << "Standard deviation for the line is: " << presStdDev << endl << endl;

	prune(timesP, tempsP, tempSlope, tempInt, 20, tempStdDev, pressP);
	prune(timesP, pressP, presSlope, presInt, 20, presStdDev, tempsP);
	findBestFitLine(timesP, tempsP, tempSP, tempIP);
	findBestFitLine(timesP, pressP, presSP, presIP);
	tempStdDev = findStdDevofRegLine(timesP, tempsP, tempSlope, tempInt);
	presStdDev = findStdDevofRegLine(timesP, pressP, presSlope, presInt);


	cout << "Recalculation after a pruning of a factor of 20 standard deviations:" << endl;
	cout << endl << "Best fit line for temp vs time:" << endl;
	cout << "y = " << tempSlope << "x + " << tempInt << endl;
	cout << "Standard deviation for the line is: " << tempStdDev << endl << endl;


	cout << "Best fit line for pressure vs time:" << endl;
	cout << "y = " << presSlope << "x + " << presInt << endl;
	cout << "Standard deviation for the line is: " << presStdDev << endl << endl;

	printTimes(timesP);
	printTemps(tempsP);
	printPress(pressP);


	prune(timesP, tempsP, tempSlope, tempInt, 20, tempStdDev, pressP);
	prune(timesP, pressP, presSlope, presInt, 20, presStdDev, tempsP);
	findBestFitLine(timesP, tempsP, tempSP, tempIP);
	findBestFitLine(timesP, pressP, presSP, presIP);
	tempStdDev = findStdDevofRegLine(timesP, tempsP, tempSlope, tempInt);
	presStdDev = findStdDevofRegLine(timesP, pressP, presSlope, presInt);

	cout << "Recalculation after a pruning of a factor of 20 standard deviations:" << endl;
	cout << endl << "Best fit line for temp vs time:" << endl;
	cout << "y = " << tempSlope << "x + " << tempInt << endl;
	cout << "Standard deviation for the line is: " << tempStdDev << endl << endl;


	cout << "Best fit line for pressure vs time:" << endl;
	cout << "y = " << presSlope << "x + " << presInt << endl;
	cout << "Standard deviation for the line is: " << presStdDev << endl << endl;

	printTimes(timesP);
	printTemps(tempsP);
	printPress(pressP);

	writeOutputFile(timesP, tempsP, pressP);


	return 0;
}
