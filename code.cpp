#include <bits/stdc++.h>

using namespace std;

#define TRACE
#ifndef ONLINE_JUDGE
#define trace(...) __f(#__VA_ARGS__, __VA_ARGS__)
template <typename Arg1>
void __f(const char* name, Arg1&& arg1){
    cerr << name << " : " << arg1 << std::endl;
}
template <typename Arg1, typename... Args>
void __f(const char* names, Arg1&& arg1, Args&&... args){
    const char* comma = strchr(names + 1, ',');cerr.write(names, comma - names) << " : " << arg1<<" | ";__f(comma+1, args...);
}
#else
#define trace(...)
#endif




const int N = 1e6 + 100;
const double EPS = 1e-10;
const int BLOCKS = 100;

typedef vector <vector <double> > v2d;

v2d dataSetSig, querySig, sim;
bool inData[N], inQuery[N];
int walkScore[N];
bool vis[N];
int n;
int nQuery, nData;
vector <double> betaValues;

struct Edge {
	int u, v;
	double w;
	Edge() {}
	Edge(int u, int v, double w) {
		this->u = u, this->v = v, this->w = w;
	}
};

void normalizeMatrix(v2d &dataSetGraph) {
	for(int j = 0; j < n; j++) {
		double sum = 0;
		for(int i = 0; i < n; i++) {
			sum += dataSetGraph[i][j];
		}
		if(sum == 0) {
			continue;
		}
		// cout << sum << endl;
		for(int i = 0; i < n; i++) {
			dataSetGraph[i][j] /= sum;
		}
	}
}

v2d mul(v2d a, v2d b) {
	v2d ret(n, vector <double> (1, 0));
	for(int i = 0; i < a.size(); i++) {
		for(int j = 0; j < b[0].size(); j++) {
			for(int k = 0; k < b.size(); k++) {
				ret[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return ret;
}

v2d mul(v2d a, double beta) {
	v2d ret(n, vector <double> (1));
	for(int i = 0; i < a.size(); i++) {
		for(int j = 0; j < a[0].size(); j++) {
			ret[i][j] = beta * a[i][j];
		}
	}
	return ret;
}

v2d add(v2d a, v2d b) {
	v2d ret(n, vector <double> (1));
	for(int i = 0; i < a.size(); i++) {
		for(int j = 0; j < a[0].size(); j++) {
			ret[i][j] = a[i][j] + b[i][j];
		}
	}
	return ret;
}

bool converged(v2d a, v2d b, double th) {
	double sum = 0;
	for(int i = 0; i < a.size(); i++) {
		for(int j = 0; j < a[0].size(); j++) {
			sum += abs(a[i][j] - b[i][j]);
		}
	}
	trace(sum, a.size(), b.size(), a[0].size(), b[0].size());
	return (sum <= th);
}

void print(v2d a) {
	for(int i = 0; i < a.size(); i++) {
		for(int j = 0; j < a[0].size(); j++) {
			cout << a[i][j] << ' ';
		}
		cout << '\n';
	}
}

v2d randomWalk(v2d &dataSetGraph, double beta) {
	v2d restart(n, vector <double> (1, 1.0 / n));
	v2d curr(n, vector <double> (1, 1.0 / n));
	int cntr = 0;
	while(true) {
		v2d tmp = mul(dataSetGraph, curr);
		tmp = mul(tmp, 1 - beta);
		v2d tmp2 = mul(restart, beta);
		v2d fin = add(tmp, tmp2);
		cntr++;
		// if(converged(fin, curr, 0.2)) {
		// 	curr = fin;
		// 	break;
		// }
		if(cntr >= 7) {
			break;
		}
		curr = fin;
	}
	cerr << "************************************************" << endl;
	// trace(cntr);
	// for(int i = 0; i < n; i++) {
	// 	cout << curr[i][0] << ' ';
	// }
	// cout << endl;
	return curr;
}


v2d betaSignature(v2d dataSetGraph) {
	v2d sig(n, vector <double>());
	for(double beta: betaValues) {
		v2d fin = randomWalk(dataSetGraph, beta);
		for(int i = 0; i < n; i++) {
			sig[i].push_back(fin[i][0]);
		}
	}
	return sig;
}

v2d betaSimilarity(v2d sigA, v2d sigB) {
	v2d sim(n, vector <double> (n, 0));
	for(int i = 0; i < n; i++) {
		// for(int j = 0; j < n; j++) {
			double sum = 0;
			for(int k = 0; k < betaValues.size(); k++) {
				sum += (sigA[i][k] - sigB[i][k]) * (sigA[i][k] - sigB[i][k]) * BLOCKS * 1.5;
			}
			sim[i][i] = 1 - sqrt(sum);
		// }
	}
	return sim;
}

int ct = 0;
bool flag;
void makeEdges(vector <Edge> &dataSetEdges, map <string, int> &mp, string inp) {
	while(cin >> inp) {
		if(inp[1] == '*') {
			break;
		}
		int u, v;
		double w;
		if(mp.find(inp) == mp.end()) {
			mp[inp] = ct++;
		}
		u = mp[inp];
		cin >> inp;
		if(mp.find(inp) == mp.end()) {
			mp[inp] = ct++;
		}
		v = mp[inp];
		cin >> inp;
		w = stod(inp);
		dataSetEdges.push_back(Edge(u, v, w));
	}
}

void makeGraph(v2d &g, vector <Edge> edges) {
	for(Edge e: edges) {
		int u, v;
		double w;
		u = e.u, v = e.v, w = e.w;
		g[u][v] = g[v][u] = w;
	}
}

int ft, md;

int dfs(int u, v2d &g, bool vis[], int d) {
	// trace(u);
	vis[u] = true;
	int ans = d;
	if(d > md) {
		md = d, ft = u;
	}
	for(int v = 0; v < n; v++) {
		if(g[u][v] != 0 and !vis[v]) {
			ans = max(dfs(v, g, vis, d + 1), ans);
		}
	}
	return ans;
}

int getRadius(v2d &g) {
	int u = -1;
	// trace(n);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(g[i][j] != 0) {
				u = i;
				trace(j);
				break;
			}
		}
		if(u != -1) {
			break;
		}
	}
	trace(u);
	cerr << "************************" << endl;
	bool vis[n];
	memset(vis, false, sizeof(vis));
	md = 0, ft = 0;
	dfs(u, g, vis, 0);
	trace(ft);
	memset(vis, false, sizeof(vis));
	md = 0;
	return dfs(ft, g, vis, 0);
}

void getNeighbour(int u, v2d &g, int d, int D, set <int> &neighbours, bool vis[]) {
	neighbours.insert(u);
	vis[u] = true;
	for(int v = 0; v < n; v++) {
		if(!vis[v] and g[u][v] != 0 and d + 1 < D) {
			getNeighbour(v, g, d + 1, D, neighbours, vis);
		}
	}
}

set <int> deltaNeighbourHood(v2d &g, int v, int r) {
	set <int> neighbours;
	bool vis[n];
	memset(vis, false, sizeof(vis));
	getNeighbour(v, g, 0, r, neighbours, vis);
	return neighbours;
}

v2d inducedSubgraph(v2d &g, int v, int r) {
	set <int> neighbours = deltaNeighbourHood(g, v, r);
	v2d tg(n, vector <double> (n, 0));
	for(int i = 0; i < n; i++) {
		if(neighbours.find(i) == neighbours.end()) {
			continue;
		}
		for(int j = 0; j < n; j++) {
			if(neighbours.find(j) != neighbours.end()) {
				tg[i][j] = tg[j][i] = g[i][j];
			}
		}
	}
	return tg;
}


void filter(set <int> &Vs, v2d queryGraph, double muv, double mus) {
	if(Vs.size() > nQuery) {
		vector <int> removed;
		for(int u: Vs) {
			double maxSim = 0;
			for(int i = 0; i < n; i++) {
				if(inQuery[i]) {
					if(i == u) {
						maxSim = max(maxSim, sim[i][i]);
					}
				}
			}
			if(maxSim < mus) {
				removed.push_back(u);
			}
		}
		for(int u: removed) {
			Vs.erase(u);
		}
	}
}

void kMatch(v2d queryGraph, v2d candidateGraph, double lambda, priority_queue <pair <int, set <int> > > &pq) {
	v2d sigQueryGraph = betaSignature(queryGraph);
	v2d sigCandidateGraph = betaSignature(candidateGraph);
	v2d bsim = betaSimilarity(queryGraph, candidateGraph);
	set <int> currNodes;
	for(int i = 0; i < n; i++) {
		if(inData[i] and inQuery[i]) {
			currNodes.insert(i);
			break;
		}
	}
	// set <int> nodes[n];
	// for(int i = 0; i < n; i++) {
	// 	nodes[i].insert(i);
	// }
	double sim[n];
	memset(sim, 0, sizeof(sim));
	// for(int i: )
	for(int t = 2; t < 2 * nQuery; t++) {
		set <int> prevNodes = currNodes;
		for(int r: prevNodes) {
			if(t % 2 == 0) {
				double currMax = 0;
				int k = -1;
				for(int p = 0; p < n; p++) {
					if(prevNodes.find(p) == prevNodes.end()) {
						if(candidateGraph[r][p] != 0 and queryGraph[r][p] != 0) {
							if(sim[r] + bsim[r][p] >= currMax) {
								currMax = sim[p] + bsim[r][p];
								k = p;
							}
						}
					}
				}
				if(k != -1) {
					sim[r] = currMax;
					currNodes.insert(k);
				}
			}
		}
	}
	trace(currNodes.size());
	pq.push({*max_element(sim, sim + n), currNodes});
}


vector <set <int> > graphMatch(v2d dataSetGraph, v2d queryGraph, double k, double muv, double mus, double lambda) {
	priority_queue <pair <int, set <int>> > pq;
	v2d sig = betaSignature(queryGraph);
	int qr = getRadius(queryGraph);
	for(int v = 0; v < n; v++) {
		if(inData[v] and inQuery[v]) {
			set <int> delta = deltaNeighbourHood(queryGraph, v, qr);
			filter(delta, queryGraph, muv, mus);
			v2d deltaGraph(n, vector <double> (n, 0));
			for(int u: delta) {
				for(int v: delta) {
					if(u != v and queryGraph[u][v] != 0) {
						deltaGraph[u][v] = deltaGraph[v][u] = queryGraph[u][v];
					}
				}
			}
			// trace(v);
			kMatch(queryGraph, deltaGraph, lambda, pq);
		}
	}
	vector <set <int> > vec;
	cout << "**********************************" << endl;
	int kcnt = 5;
	while(!pq.empty()) {
		auto x = pq.top();
		pq.pop();
		vec.push_back(x.second);
		set <int> st = x.second;
		int simValue = x.first;

		cout << simValue / st.size() << ' ';
		kcnt--;
		if(kcnt < 0) {
			break;
		}
	}
	return vec;
}

int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(NULL);
	cout.tie(NULL);
	vector <Edge> dataSetEdges, queryEdges;
	set <string> st;
	string inp;
	map <string, int> mp;
	makeEdges(dataSetEdges, mp, inp);
	flag = true;
	makeEdges(queryEdges, mp, inp);
	n = mp.size();
	cout << "number of vertices: " << n << endl;
	for(int i = 1; i <= n; i += BLOCKS) {
		betaValues.push_back((i * 1.0) / n);
	}
	cerr << "BetaValues.size: " << betaValues.size() << endl;
	v2d dataSetGraph(n, vector <double> (n, 0));
	v2d queryGraph(n, vector <double> (n, 0));
	makeGraph(dataSetGraph, dataSetEdges);
	makeGraph(queryGraph, queryEdges);

	normalizeMatrix(dataSetGraph);
	normalizeMatrix(queryGraph);

	memset(inData, false, sizeof(inData));
	memset(inQuery, false, sizeof(inQuery));

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(dataSetGraph[i][j] != 0) {
				nData += (inData[i] == false);
				inData[i] = true;
			}
			if(queryGraph[i][j] != 0) {
				nQuery += (inQuery[i] == false);
				inQuery[i] = true;
			}
		}
	}

	int qr = getRadius(queryGraph);
	int dr = getRadius(dataSetGraph);

	dataSetSig = betaSignature(dataSetGraph);
	querySig = betaSignature(queryGraph);

	sim = betaSimilarity(dataSetSig, querySig);

	double sum = 0;
	for(int i = 0; i < n; i++) {
		sum += sim[i][i];
	}
	cout << (sum / (n * 1.0)) << endl;

	vector <set <int> > vec = graphMatch(dataSetGraph, queryGraph, 0.4, 0.4, 0.4, 0.4);
}