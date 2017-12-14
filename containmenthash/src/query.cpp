#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include<fstream>
using namespace std;

int main() {
	ifstream in("../Genomes/PRJNA109269.fna");
	ofstream out("../Genomes/query.txt");
	string alpha = "ACTG";
	string txt,input;
	vector<string>str;
	vector<string>splits;

	srand(time(0));
	getline(in,txt);
	while(getline(in,txt)) {
	    if(txt[0] == '>') {
	      str.push_back(input);
	      input = "";
	    } else {
	      input += txt;
	    }
	}
	str.push_back(input);

	int split = 30;

	for(int i=0; i<str[0].size()-split+1; i=i+split) {
		splits.push_back(str[0].substr(i,split));
	}

	int len = splits.size();
	int error_cnt = 10000;

	for(int i=1;i<=error_cnt;i++) {
		int idx = rand()%len;
		string temp = splits[idx];
		temp[rand()%temp.size()] = alpha[rand() % 3];
		splits[idx] = temp;
	}

	for(int i=0; i < splits.size() ;i = i+2) {
		out << splits[i] << endl;
	}

	in.close();
	out.close();

	return 0;
}