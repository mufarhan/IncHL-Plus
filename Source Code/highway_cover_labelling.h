#ifndef HGHWAY_LABELING_H_
#define HGHWAY_LABELING_H_

#include <stdint.h>
#include <sys/time.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <thread>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <queue>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <fstream>
#include <utility>

#include "two_layer_queue.h"

//
// NOTE: Currently only unweighted and undirected graphs are supported.
//

class HighwayLabelling {
 public:
  // Constructs an index from a graph, given as a list of edges.
  HighwayLabelling(std::string filename, int k);
  HighwayLabelling();
  ~HighwayLabelling();

  void ConstructHighwayLabelling(int i, int topk[]);
  void BuildIndex(int topk[]);

  void deallocate();
  void allocate();

  void UpdateLabelling(std::string filename);
  void IncHL(int i, uint8_t *A, uint8_t *temp, int b, uint8_t d);
  bool prunable(int i, int u, uint8_t *temp, uint8_t *P);

  void SelectLandmarks_HD(int topk[]);
  int LabellingSize();

  uint8_t query(int r, int v);
  uint8_t min(uint8_t a, uint8_t b);

  void storeLabelling(std::string filename);
  void loadLabelling_Full(std::string filename, int topk[]);
  void loadLabelling_Pruned(std::string filename);

  void RemoveLandmarks(int landmarks[]);
  uint8_t UB_naive(int s, int t);
  void QueryDistance(std::string pairs, std::string output);

 private:
  int V;  // total number of vertices
  long E; // total number of edges
  int K; // total number of landmarks

  uint8_t **vertices;
  uint8_t **distances;
  uint8_t **highway;
  uint8_t *C;
  std::vector<std::vector<int> > adj;
  std::map<int, uint8_t> landmarks;

  double GetCurrentTimeSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
  }

  long GetCurrentTimeMicroSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t) 1e6 + tv.tv_usec;
  }

  long GetCurrentTimeMilliSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000LL + tv.tv_usec / 1000;
  }

  // Statistics
  double time_, time_querying_sec_;
  long time_querying_microsec_, time_querying_millisec_;
};

HighwayLabelling::HighwayLabelling() { }

HighwayLabelling::~HighwayLabelling() { }

void HighwayLabelling::deallocate() {

  for(int i = 0; i < V; i++)
    delete [] distances[i];
  delete [] distances;

  for(int i = 0; i < K; i++)
    delete [] highway[i];
  delete [] highway;
}

void HighwayLabelling::allocate() {
  // Initialization
  distances = new uint8_t*[V];
  for(int i = 0; i < V; i++) {
    distances[i] = new uint8_t[K];
    for(int j = 0; j < K; j++)
      distances[i][j] = 111;
  }

  highway = new uint8_t*[K];
  for(int i = 0; i < K; i++)
    highway[i] = new uint8_t[K];
}

HighwayLabelling::HighwayLabelling(std::string filename, int k) {
  K = k; V = 0; E = 0;

  std::ifstream ifs(filename);
  if (ifs.is_open()) {
    ifs >> V >> E;

    adj.reserve(V);

    int v, w, deg;
    while (ifs >> v >> deg) {
      adj[v].reserve(deg);
      for (int i = 0; i < deg; i++) {
        ifs >> w;
        adj[v].push_back(w);
      }
    }
    ifs.close();

    std::cout << "V : " << V << " E : " << E << std::endl << std::endl;
  } else
      std::cout << "Unable to open file" << std::endl;
}

int HighwayLabelling::LabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < K; j++) {
      if(distances[i][j] != 111)
        size++;
    }
  }
  return (V + 2 * size) / (1024 * 1024);
}

void HighwayLabelling::ConstructHighwayLabelling(int i, int topk[]) {

  uint8_t *P = new uint8_t[V];
  for(int j = 0; j < V; j++)
    P[j] = 111;

  std::queue<int> que[2];

  que[0].push(topk[i]); que[0].push(-1);
  distances[topk[i]][i] = 0; P[topk[i]] = 0; int use = 0;
  while (!que[0].empty()) {
    int u = que[use].front();
    que[use].pop();

    if(u == -1) {
      use = 1 - use;
      que[use].push(-1);
      continue;
    }

    for (int w : adj[u]) {
      if (P[w] == 111) {
        P[w] = P[u] + 1;
        if(use == 1 || landmarks.count(w) > 0)
          que[1].push(w);
        else {
          que[0].push(w);
          distances[w][i] = P[w];
        }
      }
    }
  }

  for(int j = 0; j < K; j++) {
    if(P[topk[j]] != 111) {
      highway[i][j] = P[topk[j]];
      highway[j][i] = P[topk[j]];
    }
  }

  delete [] P;
}

void HighwayLabelling::BuildIndex(int topk[]) {

  for(int i = 0; i < K; i++)
    landmarks[topk[i]] = i;

  allocate();

  // Start computing Highway Labelling (HL)
  time_ = -GetCurrentTimeSec();
  for (int i = 0; i < K; i++)
    ConstructHighwayLabelling(i, topk);
  time_ += GetCurrentTimeSec();

  std::cout << "Construction Time (sec.): " << time_ << " Labelling Size (MB): " << LabellingSize() << std::endl;
}

void HighwayLabelling::UpdateLabelling(std::string filename) {

  uint8_t *A = new uint8_t[V]; int a, b;
  uint8_t *temp = new uint8_t[V];
  for(int j = 0; j < V; j++)
    A[j] = 111;

  std::ifstream ifs(filename);

  time_ = -GetCurrentTimeSec();
  while (ifs >> a >> b) {
    adj[a].push_back(b);
    adj[b].push_back(a);

    for (int i = 0; i < K; i++) {
      int da = query(i, a);
      int db = query(i, b);

      if(da > db) {
        IncHL(i, A, temp, a, db + 1);
      } else if(da < db) {
        IncHL(i, A, temp, b, da + 1);
      }
    }
  }
  time_ += GetCurrentTimeSec();

  ifs.close();

  std::cout << "Update Time (sec.): " << time_ << " Updated Labelling Size (MB): " <<  LabellingSize() << std::endl;
  delete [] temp;
  delete [] A;
}

void HighwayLabelling::IncHL(int i, uint8_t *A, uint8_t *temp, int b, uint8_t dist) {

  // finding affected vertices
  std::queue<int> que[2];
  que[0].push(b); A[b] = dist;

  while (!que[0].empty()) {
    int u = que[0].front();
    que[0].pop();

    for (int w : adj[u]) {
      if(A[w] == 111) {
        temp[w] = query(i, w);
        if(temp[w] >= A[u] + 1)  {
          A[w] = A[u] + 1;
          que[0].push(w);
        }
      }
    }
  }

  // repairing affected vertices
  int use = 0; temp[b] = dist;
  if(prunable(i, b, temp, A)) {
    distances[b][i] = 111; que[1].push(b);
  } else {
    distances[b][i] = dist; que[0].push(b);
  }

  A[b] = 111; que[0].push(-1);
  while (!que[0].empty()) {
    int u = que[use].front();
    que[use].pop(); 

    if(u == -1) {
      use = 1 - use;
      que[use].push(-1);
      continue;
    }

    for (int w : adj[u]) {
      if(A[w] != 111) {
        temp[w] = temp[u] + 1;
        if(use == 1 || prunable(i, w, temp, A)) {
          distances[w][i] = 111;
	  que[1].push(w);
	} else {
	  distances[w][i] = temp[w];
	  que[0].push(w);
	}
	A[w] = 111;
      }
    }
  }

  while (!que[1].empty()) {
    int u = que[1].front();
    que[1].pop();

    if(u == -1)
      continue;

    for (int w : adj[u]) {
      if(A[w] != 111) {
	A[w] = 111;
	distances[w][i] = 111;

	que[1].push(w);
      }
    }
  }
}

bool HighwayLabelling::prunable(int i, int u, uint8_t *temp, uint8_t *A) {

  if(landmarks.count(u) > 0) {
    highway[i][landmarks[u]] = temp[u];
    return true;
  } else {
    for (int w : adj[u]) {
      if(A[w] == 111) {
        if(temp[w] == temp[u] - 1 && distances[w][i] == 111)
          return true;
      }
    }
  }
  return false;
}

uint8_t HighwayLabelling::query(int r, int v) {

  uint8_t m = 111;
  for(int i = 0; i < K; i++) {
    m = min(m, distances[v][i] + highway[r][i]);
  }
  return m;
}

uint8_t HighwayLabelling::min(uint8_t a, uint8_t b) {
  return (a < b) ? a : b;
}

void HighwayLabelling::RemoveLandmarks(int landmarks[]) {

  for(int i = 0; i < K; i++) {
    for (int v : adj[landmarks[i]]) {
      adj[v].erase(std::remove(adj[v].begin(), adj[v].end(), landmarks[i]), adj[v].end());
      adj[v].shrink_to_fit();
    }
    adj[landmarks[i]].clear();
    adj[landmarks[i]].shrink_to_fit();
  }
}

uint8_t HighwayLabelling::UB_naive(int s, int t) {

  uint8_t m = 111; int i, j;
  for(i = 0; i < C[s]; i++) {
    for (j = 0; j < C[t]; j++)
      m = min(m, distances[s][i] + highway[vertices[s][i]][vertices[t][j]] + distances[t][j]);
  }
  return m;
}

void HighwayLabelling::QueryDistance(std::string pairs, std::string output) {
  std::vector<TwoLayerQueue> qque; std::vector<uint8_t> qdist[2];

  qdist[0].resize(V, 111); qdist[1].resize(V, 111);
  qque.push_back(TwoLayerQueue(V)); qque.push_back(TwoLayerQueue(V));

  time_querying_millisec_ = 0; int s = 0, t = 0; int total = 0;
  std::ifstream ifs(pairs); std::ofstream ofs(output);
  while(ifs >> s >> t) { total++;

    double a = -GetCurrentTimeMilliSec();

    uint8_t dist_upper = UB_naive(s, t);
    uint8_t res = dist_upper, dis[2] = {0, 0};
    for (int dir = 0; dir < 2; dir++){
      int v = dir == 0 ? s : t;
      qque[dir].clear();
      qque[dir].push(v);
      qque[dir].next();
      qdist[dir][v] = 0;
    }

    while (!qque[0].empty() && !qque[1].empty()) {
      int use = (qque[0].size() <= qque[1].size()) ? 0 : 1;
      dis[use]++;

      if (dis[0] + dis[1] == dist_upper) {
        res = dis[0] + dis[1];
        goto LOOP_END;
      }

      while (!qque[use].empty()) {

        int v = qque[use].front();
        qque[use].pop();

        for (int w : adj[v]) {

          uint8_t &src_d = qdist[    use][w];
          uint8_t &dst_d = qdist[1 - use][w];
          if (src_d != 111) continue;
          if (dst_d != 111) {
            res = qdist[use][v] + 1 + dst_d;
            goto LOOP_END;
          } else {
            qque[use].push(w);
            qdist[use][w] = qdist[use][v] + 1;
          }
        }
      }
      qque[use].next();
    }
    LOOP_END:

    a += GetCurrentTimeMilliSec();
    time_querying_millisec_ += a;

    for (int dir = 0; dir < 2; dir++) {
      for (int v : qque[dir]) {
        qdist[dir][v] = 111;
      }
      qque[dir].clear();
    }
    ofs << s << " " << t << " " << (int) min(res, dist_upper) << "\n";
  }
  ifs.close();
  ofs.close();

  std::cout << "Average Query Time (ms) : " << (double) time_querying_millisec_ / total << std::endl;

  for(int i = 0; i < V; i++) {
    delete [] distances[i];
    delete [] vertices[i];
  }
  delete [] distances;
  delete [] vertices;
  delete [] C;

  for(int i = 0; i < K; i++)
    delete [] highway[i];
  delete [] highway;
}

void HighwayLabelling::SelectLandmarks_HD(int topk[]) {
  std::vector<std::pair<int, int> > deg(V); long sum = 0;
  for (int v = 0; v < V; v++) {
    deg[v] = std::make_pair(adj[v].size(), v);
    sum = sum + adj[v].size();
  }
  std::sort(deg.rbegin(), deg.rend());

  for (int v = 0; v < K; v++)
    topk[v] = deg[v].second;
}

void HighwayLabelling::storeLabelling(std::string filename) {
  std::ofstream ofs(std::string(filename) + std::string("_index"));

  for(int i = 0; i < V; i++) {
    uint8_t C = 0;
    for(int j = 0; j < K; j++) {
      if(distances[i][j] != 111)
        C++;
    }

    ofs.write((char*)&C, sizeof(C));
    for(uint8_t j = 0; j < K; j++) {
      if(distances[i][j] != 111) {
        ofs.write((char*)&j, sizeof(j));
        ofs.write((char*)&distances[i][j], sizeof(distances[i][j]));
      }
    }
  }
  ofs.close();

  ofs.open(std::string(filename) + std::string("_highway"));
  for(int i = 0; i < K; i++) {
    for(int j = 0; j < K; j++) {
      ofs.write((char*)&highway[i][j], sizeof(highway[i][j]));
    }
  }
  ofs.close();
}

void HighwayLabelling::loadLabelling_Full(std::string filename, int topk[]) {

  for(int i = 0; i < K; i++)
    landmarks[topk[i]] = i;

  std::ifstream ifs(std::string(filename) + std::string("_index"));

  time_ = -GetCurrentTimeSec();

  distances = new uint8_t*[V];
  for(int i = 0; i < V; i++) {
    distances[i] = new uint8_t[K];
    for(int j = 0; j < K; j++)
      distances[i][j] = 111;
  }

  uint8_t C, idx;
  for(int i = 0; i < V; i++) {
    ifs.read((char*)&C, sizeof(C));
    for(uint8_t j = 0; j < C; j++) {
      ifs.read((char*)&idx, sizeof(idx));
      ifs.read((char*)&distances[i][idx], sizeof(distances[i][idx]));
    }
  }
  ifs.close();

  ifs.open(std::string(filename) + std::string("_highway"));
  highway = new uint8_t*[K];
  for(uint8_t i = 0; i < K; i++) {
    highway[i] = new uint8_t[K];
    for(uint8_t j = 0; j < K; j++)
      ifs.read((char*)&highway[i][j], sizeof(highway[i][j]));
  }
  ifs.close();

  time_ += GetCurrentTimeSec();
}

void HighwayLabelling::loadLabelling_Pruned(std::string filename) {
  std::ifstream ifs(std::string(filename) + std::string("_index"));

  C = new uint8_t[V];
  vertices = new uint8_t*[V];
  distances = new uint8_t*[V];

  for(int i = 0; i < V; i++) {
    ifs.read((char*)&C[i], sizeof(C[i]));
    vertices[i] = new uint8_t[C[i]];
    distances[i] = new uint8_t[C[i]];
    for(uint8_t j = 0; j < C[i]; j++) {
      ifs.read((char*)&vertices[i][j], sizeof(vertices[i][j]));
      ifs.read((char*)&distances[i][j], sizeof(distances[i][j]));
    }
  }
  ifs.close();

  ifs.open(std::string(filename) + std::string("_highway"));
  highway = new uint8_t*[K];
  for(uint8_t i = 0; i < K; i++) {
    highway[i] = new uint8_t[K];
    for(uint8_t j = 0; j < K; j++)
      ifs.read((char*)&highway[i][j], sizeof(highway[i][j]));
  }
  ifs.close();
}

#endif  // PRUNED_LANDMARK_LABELING_H_
