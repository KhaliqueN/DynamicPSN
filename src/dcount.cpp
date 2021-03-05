#include <iostream>
#include <fstream>
#include<stdio.h>
#include<string.h>
#include<vector>
#include<map>
#include<utility>
#include<sstream>
#include<algorithm>
#include<set>
#include <getopt.h>
#include<time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define PROGRESS_INFO 1

using namespace std;


bool ONLY_CAUSAL_EVENTS = false; //e.g., to have temporal graphlet 1223 need that event 23 not to existed at time when 12 occurred (but it could exist some time before it);

struct TEvent
{
	int time;
	int v1;
	int v2;
};

//for sorting vector of events based on time
struct sort_events_pred {
	bool operator()(const TEvent &left, const TEvent &right) {
		return left.time < right.time;
	}
};

class TemporalNetwork
{
	//"temporal adjacency list" - for each node, list of adjacent nodes with corresponding times of contact
	map<int, vector<pair<int, int> > > adjList;
	bool isDirected;
	map<string, int> nodeIds;
	vector<TEvent> EventList;
	map<string, vector<TEvent> > nextEventsCache;//cache - for each event list of next events within deltaT

public:

	TemporalNetwork(bool isDirected, string ifname, string delimiter)
	{

		cout << "Initializing temporal network...";

		this->isDirected = isDirected;

		ifstream myfile(ifname.c_str());
		int currId = 1;
		string line;
		if (myfile.is_open()){
			string node1, node2, t_itime;
			double itime;
			int node1Id, node2Id;
			while (myfile){
				string s;
				if (!getline(myfile, s)) break;
				int i = 0, pos;

				while ((pos = s.find(delimiter)) != std::string::npos) {
					if (i == 0)
						node1 = s.substr(0, pos);
					else if (i == 1)
						node2 = s.substr(0, pos);
					s.erase(0, pos + delimiter.length());
					++i;
				}
				if (i != 2)
				{
					cout << "ERROR: Incorrect file format\n";
					return;
				}
				t_itime = s;
				itime = atoi(t_itime.c_str());

				//cout<<"node1="<<node1<<" | node2="<<node2<<endl;
				//cout<<endl;
				/*
				map<string,int>::const_iterator iter;
				for(iter=nodeIds.begin();iter!=nodeIds.end();iter++){
				cout<<"Key:"<<iter->first<<" |values: "<<iter->second<<endl;
				}
				cout<<endl;
				*/
				/////
				if (this->nodeIds.find(node1) == this->nodeIds.end())
				{
					this->nodeIds[node1] = currId++;
				}
				if (this->nodeIds.find(node2) == this->nodeIds.end())
				{
					this->nodeIds[node2] = currId++;
				}

				node1Id = this->nodeIds[node1];
				node2Id = this->nodeIds[node2];

				if (node1Id != node2Id){
					this->adjList[node1Id].push_back(std::make_pair(itime, node2Id));
					if (!this->isDirected)
						this->adjList[node2Id].push_back(std::make_pair(itime, node1Id));
				}
				else{ cout << "self loop is found. node-" << node1.c_str() << endl; }

			}
			myfile.close();
		}
		else cout << "ERROR: Unable to open file '" << ifname << "'\n";

		cout << "Network initialized: " << this->getNodeCount() << " nodes and " << this->getEventCount() << " events" << endl;
	}

	TemporalNetwork(bool isDirected, vector<string> ifnames, string delimiter)
	{

		cout << "Initializing network from " << ifnames.size() << " snapshots..." << endl;

		this->isDirected = isDirected;

		int iSnapshot = 0;
		int currId = 1;

		for (vector<string>::iterator it = ifnames.begin(); it != ifnames.end(); ++it)
		{
			ifstream myfile((*it).c_str());
			int i = 0;
			string line;

			if (myfile.is_open()){
				string node1, node2, t_itime;
				double itime;
				int node1Id, node2Id;
				while (myfile){

					string s;
					if (!getline(myfile, s)) break;
					int i = 0, pos;

					while ((pos = s.find(delimiter)) != std::string::npos) {
						node1 = s.substr(0, pos);
						s.erase(0, pos + delimiter.length());
						++i;
					}
					if (i != 1)
					{
						cout << "ERROR: Incorrect file format\n";
						return;
					}
						
					node2 = s;

					//cout<<"node1="<<node1<<" | node2="<<node2<<endl;
					//cout<<endl;

					//map<string,int>::const_iterator iter;
					//for(iter=nodeIds.begin();iter!=nodeIds.end();iter++){
					//  cout<<"Key:"<<iter->first<<" |values: "<<iter->second<<endl;
					//}
					//cout<<endl;

					if (this->nodeIds.find(node1) == this->nodeIds.end())
					{
						this->nodeIds[node1] = currId++;

					}
					if (this->nodeIds.find(node2) == this->nodeIds.end())
					{
						this->nodeIds[node2] = currId++;

					}

					node1Id = this->nodeIds[node1];
					node2Id = this->nodeIds[node2];

					//cout << "Added event: " << iSnapshot << " snapshot | " << node1Id << "(" << node1 << ") - " << node2Id << "(" << node2 << ");" << endl;

					if (node1Id != node2Id){
						this->adjList[node1Id].push_back(std::make_pair(iSnapshot, node2Id));
						if (!this->isDirected)
							this->adjList[node2Id].push_back(std::make_pair(iSnapshot, node1Id));
					}
					else{ cout << "self loop is found. node-" << node1.c_str() << endl; }

				}
				myfile.close();
				++iSnapshot;
			}
			else cout << "ERROR: Unable to open file '" << *it << "'\n";
		}
		cout << "Network initialized: " << this->getNodeCount() << " nodes and " << this->getEventCount() << " events" << endl;
	}


	int getNodeCount()
	{
		return this->nodeIds.size();
	}

	int getEventCount()
	{
		return this->getEvents().size();
	}

	string get_nodeNames(int id){
		string iname;
		map<string, int>::const_iterator it;
		for (it = nodeIds.begin(); it != nodeIds.end(); ++it){
			if (it->second == id){
				iname = it->first;
				break;
			}
		}
		return iname;
	}

	vector<TEvent> getEvents()
	{

		if (this->EventList.size() == 0)
		{

			for (map<int, vector<pair<int, int> > >::iterator outer_iter = this->adjList.begin(); outer_iter != this->adjList.end(); ++outer_iter) {
				for (vector<pair<int, int> >::iterator inner_iter = outer_iter->second.begin(); inner_iter != outer_iter->second.end(); ++inner_iter) {

					if (!this->isDirected && outer_iter->first >= inner_iter->second)
						continue;

					TEvent tempEvent;
					tempEvent.time = inner_iter->first;
					tempEvent.v1 = outer_iter->first;
					tempEvent.v2 = inner_iter->second;
					this->EventList.push_back(tempEvent);
				}
			}

			//for (vector<TEvent>::iterator iter = this->EventList.begin(); iter != this->EventList.end(); ++iter) {
			//	cout << "Event: " << (*iter).time << "| " << get_nodeNames((*iter).v1) << "-" << get_nodeNames((*iter).v2) << endl;
			//}

			//sort events based on timestamp
			std::sort(this->EventList.begin(), this->EventList.end(), sort_events_pred());

		}

		return this->EventList;
	}



	vector<TEvent> getNextEvents(TEvent prevEvent, int deltaT)
	{
		ostringstream Convert;
		Convert << prevEvent.time << "_" << prevEvent.v1 << "_" << prevEvent.v2;
		string prevEventCode = Convert.str();
		if (this->nextEventsCache.find(prevEventCode) != this->nextEventsCache.end())
			return this->nextEventsCache[prevEventCode];

		vector<TEvent> nextEvents;

		//Note: assuming events stored ordered by time

		map<int, bool> existedBefore;//whether an event existed within  (prevEvent.time - deltaT, prevEvent] - it is used to determine whether next event existed at the moment of prevEvent, or it "appeared" only after it
		for (vector<pair<int, int> >::iterator it = this->adjList[prevEvent.v1].begin(); it != this->adjList[prevEvent.v1].end(); ++it) {
			if (it->first > prevEvent.time && it->first <= prevEvent.time + deltaT && (!ONLY_CAUSAL_EVENTS || !existedBefore[it->second])){
				TEvent tempEvent;
				tempEvent.time = it->first;
				tempEvent.v1 = prevEvent.v1;
				tempEvent.v2 = it->second;
				nextEvents.push_back(tempEvent);
			}
			else if (ONLY_CAUSAL_EVENTS && it->second != prevEvent.v2 && it->first > prevEvent.time - deltaT && it->first <= prevEvent.time) {
				existedBefore[it->second] = true;
			}
			else if (it->first > prevEvent.time + deltaT)//since events are ordered
				break;

		}

		existedBefore.clear();
		for (std::vector<pair<int, int> >::iterator it = this->adjList[prevEvent.v2].begin(); it != this->adjList[prevEvent.v2].end(); ++it) {
			if (it->first > prevEvent.time && it->first <= prevEvent.time + deltaT && (this->isDirected || prevEvent.v1 != it->second) && (!ONLY_CAUSAL_EVENTS || !existedBefore[it->second])){
				TEvent tempEvent;
				tempEvent.time = it->first;
				tempEvent.v1 = prevEvent.v2;
				tempEvent.v2 = it->second;
				nextEvents.push_back(tempEvent);
			}
			else if (ONLY_CAUSAL_EVENTS && it->second != prevEvent.v1 && it->first > prevEvent.time - deltaT && it->first <= prevEvent.time) {
				existedBefore[it->second] = true;
			}
			else if (it->first > prevEvent.time + deltaT)//since events are ordered
				break;
		}

		this->nextEventsCache[prevEventCode] = nextEvents;
		return this->nextEventsCache[prevEventCode];

	}


};

struct encoding_nodeCodes{
	string encoding;
	map<int, int> nodeCodes;
};

encoding_nodeCodes encodeGraphlet(vector<TEvent> *Graphlet){
	string encoding = "";
	map<int, int> nodeCodes;
	int i = 1;
	//for all e in Graphlet do
	for (int j = 0; j < Graphlet->size(); j++){
		//u,v <- e.EndPoints()
		ostringstream con1;
		ostringstream con2;
		int u = Graphlet->at(j).v1;
		int v = Graphlet->at(j).v2;

		bool update = false;
		if (nodeCodes.find(u) == nodeCodes.end()){
			nodeCodes[u] = i;
			++i;
			update = true;
		}
		if (nodeCodes.find(v) == nodeCodes.end()){
			nodeCodes[v] = i;
			++i;
			update = true;
		}
		if (i == 4 && update == true) {
			//replace "121212...1213" with "121212...1223"
			if ((nodeCodes[u] == 1 && nodeCodes[v] == 3) || (nodeCodes[u] == 3 && nodeCodes[v] == 1)){
				int key1, key2;

				for (map<int, int>::const_iterator it = nodeCodes.begin(); it != nodeCodes.end(); ++it){
					if (it->second == 1)
						key1 = it->first;
					else if (it->second == 2)
						key2 = it->first;
				}

				nodeCodes[key1] = 2;
				nodeCodes[key2] = 1;
			}
		}
		ostringstream tmp1;
		ostringstream tmp2;
		tmp1 << min(nodeCodes[u], nodeCodes[v]);
		tmp2 << max(nodeCodes[u], nodeCodes[v]);
		encoding.append(tmp1.str());
		encoding.append(tmp2.str());


	}

	if (encoding == "1211")
	{
		for (int j = 0; j < Graphlet->size(); j++){
			cout << "(" << Graphlet->at(j).time << ", " << Graphlet->at(j).v1 << ", " << Graphlet->at(j).v2 << "); ";
		}
		cout << "\n";
	}
	encoding_nodeCodes ireturn;
	ireturn.encoding = encoding;
	ireturn.nodeCodes.insert(nodeCodes.begin(), nodeCodes.end());

	return ireturn;
}

// given encoding of a graphlet, find it encoding after addition of event: given graphlet encoding and new event -> find their encoding
encoding_nodeCodes encodeGraphletFromPrefix(encoding_nodeCodes prefixEncoding, TEvent lastEvent){

	string encoding = prefixEncoding.encoding;
	map<int, int> nodeCodes = prefixEncoding.nodeCodes;
	int i = 1;

	int maxCode = 0;
	for (map<int, int>::const_iterator it = nodeCodes.begin(); it != nodeCodes.end(); ++it){
		if (it->second > maxCode)
			maxCode = it->second;
	}

	bool update = false;
	if (nodeCodes.find(lastEvent.v1) == nodeCodes.end()){
		++maxCode;
		nodeCodes[lastEvent.v1] = maxCode;
		update = true;
	}
	if (nodeCodes.find(lastEvent.v2) == nodeCodes.end()){
		++maxCode;
		nodeCodes[lastEvent.v2] = maxCode;
		update = true;
	}

	if (maxCode == 3 && update == true) {
		//replace "121212...1213" with "121212...1223"
		if ((nodeCodes[lastEvent.v1] == 1 && nodeCodes[lastEvent.v2] == 3) || (nodeCodes[lastEvent.v1] == 3 && nodeCodes[lastEvent.v2] == 1)){
			int key1, key2;

			for (map<int, int>::const_iterator it = nodeCodes.begin(); it != nodeCodes.end(); ++it){
				if (it->second == 1)
					key1 = it->first;
				else if (it->second == 2)
					key2 = it->first;
			}

			nodeCodes[key1] = 2;
			nodeCodes[key2] = 1;
		}
	}

	ostringstream tmp1;
	ostringstream tmp2;
	tmp1 << min(nodeCodes[lastEvent.v1], nodeCodes[lastEvent.v2]);
	tmp2 << max(nodeCodes[lastEvent.v1], nodeCodes[lastEvent.v2]);
	encoding.append(tmp1.str());
	encoding.append(tmp2.str());

	encoding_nodeCodes ireturn;
	ireturn.encoding = encoding;
	ireturn.nodeCodes.insert(nodeCodes.begin(), nodeCodes.end());


	//cout << "Encoding  '" << prefixEncoding.encoding << "' extended by: " << lastEvent.time << "| " << lastEvent.v1 << ", " << lastEvent.v2 << " to '" << encoding << "'"<< endl;

	return ireturn;
}



//counting temporal graphlets with up to k_max events starting from a given prefix
void graphletCountFromPrefix(TemporalNetwork *G, TEvent lastEvent, encoding_nodeCodes prefixEN, map<string, int> *graphletCounts, map<pair<int, string>, int> *NCounts, int n_max, int k_max, int deltaT){
	encoding_nodeCodes e_n = encodeGraphletFromPrefix(prefixEN, lastEvent);

	string prefixEncoding = e_n.encoding;

	//#pragma omp critical
	(*graphletCounts)[prefixEncoding]++;

	//for all v belongs to prefNodes, do. if |prefNodes| =2 then nodePosition <- '1' else nodePosition <- nodeSet[v]

	for (map<int, int>::iterator it = e_n.nodeCodes.begin(); it != e_n.nodeCodes.end(); ++it){
		string nodePosition;
		int v = it->first;
		if (e_n.nodeCodes.size() == 2) nodePosition = "1";
		else{
			ostringstream convert;
			convert << e_n.nodeCodes[v];
			nodePosition = convert.str();
		}
		string orbit = prefixEncoding + "_" + nodePosition;

		//#pragma omp critical
		(*NCounts)[make_pair(v, orbit)]++;
	}

	/*
	cout << "Next events from (" << lastEvent.time << ", " << lastEvent.v1 << ", " << lastEvent.v2 << ") :";
	for (int m = 0; m < next_events.size(); m++){
	cout << "(" << next_events.at(m).time << ", " << next_events.at(m).v1 << ", " << next_events.at(m).v2 << "), ";
	}
	cout << "\n";
	*/

	if (prefixEncoding.length() < 2 * k_max){
		vector<TEvent> next_events = G->getNextEvents(lastEvent, deltaT);

		for (int m = 0; m < next_events.size(); m++){
			TEvent nextEvent = next_events.at(m);
			if (e_n.nodeCodes.size() < n_max || (e_n.nodeCodes.find(nextEvent.v1) != e_n.nodeCodes.end() && e_n.nodeCodes.find(nextEvent.v2) != e_n.nodeCodes.end()))
				graphletCountFromPrefix(G, nextEvent, e_n, graphletCounts, NCounts, n_max, k_max, deltaT);
		}
	}
}


int* initialize_gdvs(int size){
	int *arr = new int[size];
	for (int i = 0; i < size; i++){
		arr[i] = 0;
	}
	return arr;
}

void graphletCount(TemporalNetwork *G, int n_max, int k_max, int deltaT, ofstream *outf, ofstream *gdvf, int gdvFlag, string orbitF){




	map<string, int> counts;
	map<pair<int, string>, int> NCounts;
	vector<TEvent> events = G->getEvents();


	cout << "Counting graphlets with k=" << k_max << " and n=" << n_max;
	if (ONLY_CAUSAL_EVENTS)
		cout << " (counting only 'causal' graphlets)";
	cout << "..." << endl;

	time_t start, end;
	time(&start);




	//for all e belonging to G.events() do

	#ifdef _OPENMP
		int nProcessors = omp_get_max_threads();
		cout << nProcessors << " processors" << endl;
		omp_set_num_threads(nProcessors);
	#endif
	
	int E = events.size();
	time_t start_time = time(0);
    clock_t t_start = clock();
	#pragma omp parallel for schedule(static, 100)
	for (int i = 0; i < E; ++i) {

		#if PROGRESS_INFO
		fprintf(stderr, "\revent # %5d (timestamp # %3d) - %.1f%%  [%d min elapsed]", i, events[i].time,
		  100.0 * (float)i / E, (time(0) - start_time) / 60);
		#endif

		encoding_nodeCodes e_n;

		//	  cout << "starting from " << (*iter).time << " - " << (*iter).v1 << " - " << (*iter).v2 << "\n";

		graphletCountFromPrefix(G, events[i], e_n, &counts, &NCounts, n_max, k_max, deltaT);
	}

	time(&end);
	double diff = difftime(end, start);
    cout << endl << "Finished graphlet counting, CPU time: " << (clock() - t_start) / (double)CLOCKS_PER_SEC << " seconds; wall time:" << diff << " seconds" << endl;



	cout << "Generating output..." << endl;

	//write the dictionary of counts to the output file
	for (map<string, int>::const_iterator it = counts.begin(); it != counts.end(); ++it){
		(*outf) << it->first << "\t" << it->second << endl;
	}
	if (gdvFlag != 0){

		ifstream infile(orbitF.c_str());
		map<string, int> gdvs;
		int icounter = 0;
		while (infile)
		{
			string s;
			if (!getline(infile, s)) break;
			istringstream ss(s);
			int cc = 0;
			while (ss)
			{
				string s;
				if (!getline(ss, s, '\t')) break;
				if (cc == 0) { gdvs[s] = icounter; cc++; icounter++; }
			}
		}


		map<int, int*> node_gdvs;
		//generate a gdv row for each node
		for (map<pair<int, string>, int>::const_iterator it = NCounts.begin(); it != NCounts.end(); ++it){
			//if the node does not exist in the vector
			if (node_gdvs.find(it->first.first) == node_gdvs.end()){

				node_gdvs[it->first.first] = initialize_gdvs(gdvs.size());
				if (gdvs.find(it->first.second) != gdvs.end()){
					node_gdvs[it->first.first][gdvs[it->first.second]] = it->second;
				}
				else{ cout << "the input k_max and n_max do not match the orbit file" << endl; }
			}
			else{
				//find the orbit in gdvs
				if (gdvs.find(it->first.second) != gdvs.end()){
					node_gdvs[it->first.first][gdvs[it->first.second]] = it->second;
				}
				else{ cout << "the input k_max and n_max do not match the orbit file" << endl; }
			}
		}


		map<int, int*>::iterator it;
		for (it = node_gdvs.begin(); it != node_gdvs.end(); it++){
			// the first column is a list of node names
			(*gdvf) << G->get_nodeNames(it->first) << " ";
			for (int jj = 0; jj < gdvs.size(); jj++){
				if (jj != gdvs.size() - 1) (*gdvf) << it->second[jj] << " ";
				else (*gdvf) << it->second[jj] << "\n";
			}
		}
		gdvf->close();
	}
}

vector<string> split(const string &s, char delim) {
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}

	return elems;
}

void generate_vector_file(string vectorF, string gcountF, string graphletCountOutputF){

	cout << "Generating output in vector form..." << endl;

	ifstream infile(vectorF.c_str());
	vector<string> patterns;
	while (infile)
	{
		string s;
		if (!getline(infile, s)) break;
		s.erase(s.find_last_not_of(" \n\r\t") + 1);
		patterns.push_back(s);
	}
	if (!infile.eof())
	{
	}
	//search for the corresponding patterns stored in the vector in the output file _counts.txt
	std::stringstream kn_string;
	map<string, string> pattern_map;
	ifstream ofile(graphletCountOutputF.c_str());
	while (ofile){
		string s;
		if (!getline(ofile, s)) break;
		istringstream ss(s);
		string ipattern, icount;
		ss >> ipattern >> icount;
		pattern_map[ipattern] = icount;
	}
	//open _gcounts_k_n_t.txt file
	ofstream gfile;
	gfile.open(gcountF.c_str());
	for (int ii = 0; ii < patterns.size(); ii++){
		if (pattern_map.find(patterns.at(ii)) != pattern_map.end()){
			gfile << pattern_map[patterns.at(ii)] << endl;
		}
		else{//output 0
			gfile << "0" << endl;
		}
	}
	gfile.close();
}


int main(int argc, char *argv[]){
	//measure running time
	time_t start, end;
	time(&start);


	TemporalNetwork *G;
	int k_max, n_max;
	int deltaT;

	string outputF, inputF;
	string vectorF = "0";
	string orbitF = "0";
	string k_tmp, n_tmp, delta_tmp, tmp_separator;
	string delimiter = "\t";


	if(argc<6) {
	cout << "ERROR: Incorrect number of input parameters" << "\n";
			cout << "USAGE: " << argv[0] << "<input_file[s]> <k_max> <n_max> <delta_t> <output base filename> -d delimiter [tab by default] -v <file with temporal graphlet list> -g <file with temporal orbit list> -t <count only causal type> -h help" << endl;
	  return 1;
	  }
	else
	{
	inputF = argv[1];
	k_tmp = argv[2];
	n_tmp = argv[3];
	delta_tmp = argv[4];
	outputF = argv[5];
	}  


	int delimiterFlag = 0;
	int vectorFlag = 0;
	int gdvFlag = 0;
	int c;
	// optional argument. d stands for delimiter. v stands for vectorFile. g stands for gdvFlag. if g is not selected, then gdv files will not be generated.
	while ((c = getopt(argc, argv, "d:v:g:th")) != -1){
		switch (c){
		case 'd':
			delimiterFlag = 1;
			delimiter = optarg;
			if (delimiter == "\\t" || delimiter == "tab")
			    delimiter = "\t";
			break;
		case 'v':
			vectorFlag = 1;
			vectorF = optarg;
			break;
		case 'g':
			gdvFlag = 1;
			orbitF = optarg;
			break;
		case 't':
			ONLY_CAUSAL_EVENTS = true;
			break;
		case 'h':
			cout << "USAGE: " << argv[0] << "<input_file[s]> <k_max> <n_max> <delta_t> <output base filename> -d delimiter [tab by default] -v <file with temporal graphlet list> -g <file with temporal orbit list> -t <count only causal type>  -h help" << endl;
			return 0;
			break;
		default:
			abort();
		}
	}

	if (inputF.find(',') != string::npos){
		G = new TemporalNetwork(false, split(inputF, ','), delimiter);
	}
	else{
		G = new TemporalNetwork(false, inputF, delimiter);
	}

	k_max = atoi(k_tmp.c_str());
	n_max = atoi(n_tmp.c_str());
	deltaT = atoi(delta_tmp.c_str());

	std::stringstream kn_string;
	kn_string << k_max << "_" << n_max << "_" << deltaT;

	//graphlet count
	string graphletCountOutputF = outputF + "_d" + (ONLY_CAUSAL_EVENTS ? "c" : "") + "counts_" + kn_string.str() + ".txt";
	ofstream graphletCountOutputFile(graphletCountOutputF.c_str());

	//gdv orbits
	string gdvF = outputF + "_d" + (ONLY_CAUSAL_EVENTS ? "c" : "") + "gdv_" + kn_string.str() + ".txt";
	ofstream gdvOutputFile(gdvF.c_str());

	graphletCount(G, n_max, k_max, deltaT, &graphletCountOutputFile, &gdvOutputFile, gdvFlag, orbitF);

	//generate the output in vector form
	if (vectorFlag != 0){
		string gcountF = outputF + "_d" + (ONLY_CAUSAL_EVENTS ? "c" : "") + "gcounts_" + kn_string.str() + ".txt";
		generate_vector_file(vectorF, gcountF, graphletCountOutputF);
	}

	time(&end);
	double diff = difftime(end, start);

	//file storing running time
	string timeF = outputF + (ONLY_CAUSAL_EVENTS ? "_c" : "") + "_time_" + kn_string.str() + ".txt";
	ofstream timeFile(timeF.c_str());
	timeFile << "Running Time(sec): " << diff << endl;

	cout << "Total running time: " << diff << " seconds" << endl;

	return 0;
}
