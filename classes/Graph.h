// Original code by Gonçalo Leão
// Updated by DA 2023/2024 Team

#ifndef DA_TP_CLASSES_GRAPH
#define DA_TP_CLASSES_GRAPH

#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include <limits>
#include <algorithm>
#include "../utils/MutablePriorityQueue.h"

template<class T>
class Edge;

#define INF std::numeric_limits<double>::max()

/************************* Vertex  **************************/

template<class T>
class Vertex {
public:
    Vertex() {}

    Vertex(T in);

    Vertex(T in, double longitude, double latitude);

    bool operator<(Vertex<T> &vertex) const; // required by MutablePriorityQueue

    T getInfo() const;

    std::vector<Edge<T> *> getAdj() const;

    bool isVisited() const;

    bool isProcessing() const;

    unsigned int getIndegree() const;

    double getDist() const;

    Edge<T> *getPath() const;

    std::vector<Edge<T> *> getIncoming() const;

    double getIncFlow() const;

    double getLongitude() const;

    double getLatitude() const;

    void reverseAdj();

    void reverseIncoming();

    void setInfo(T info);

    void setVisited(bool visited);

    void setProcessing(bool processing);

    void setIndegree(unsigned int indegree);

    void setDist(double dist);

    void setPath(Edge<T> *path);

    void setIncFlow(double incFlow);

    Edge<T> *addEdge(Vertex<T> *dest, double w);

    bool removeEdge(T in);

    void removeOutgoingEdges();

    friend class MutablePriorityQueue<Vertex>;

    int queueIndex = 0;
protected:
    T info;                // info node
    double longitude;
    double latitude;
    std::vector<Edge<T> *> adj;  // outgoing edges
    std::vector<Edge<T> *> incoming; // incoming edges
    // auxiliary fields
    bool visited = false; // used by DFS, BFS, Prim ...
    bool processing = false; // used by isDAG (in addition to the visited attribute)
    unsigned int indegree; // used by topsort
    double dist = 0;
    Edge<T> *path = nullptr;
    // required by MutablePriorityQueue and UFDS
    double incFlow = 0;

    void deleteEdge(Edge<T> *edge);
};


/********************** Edge  ****************************/

template<class T>
class Edge {
public:
    Edge(Vertex<T> *orig, Vertex<T> *dest, double w);

    Vertex<T> *getDest() const;

    double getWeight() const;

    bool isSelected() const;

    Vertex<T> *getOrig() const;

    Edge<T> *getReverse() const;

    double getFlow() const;

    double getPreviousFlow() const;

    double getTempWeight() const;

    void setSelected(bool selected);

    void setReverse(Edge<T> *reverse);

    void setFlow(double flow);

    void setPreviousFlow(double previousFlow);

    void setWeight(double weight);

    void setTempWeight(double tempWeight);

protected:
    Vertex<T> *dest; // destination vertex
    double weight; // edge weight, can also be used for capacity

    // auxiliary fields

    // used for bidirectional edges
    Edge<T> *reverse = nullptr;
    Vertex<T> *orig;

    bool selected = false;
    double flow; // for flow-related problems
    double previousFlow;
    double tempWeight = 0;
};


/********************** Graph  ****************************/

template<class T>
class Graph {
public:
    ~Graph();

    /*
    * Auxiliary function to find a vertex with a given the content.
    */
    Vertex<T> *findVertex(const T &in) const;

    bool removeVertex(const T &in);

    /*
     * Adds an edge to a graph (this), given the contents of the source and
     * destination vertices and the edge weight (w).
     * Returns true if successful, and false if the source or destination vertex does not exist.
     */
    bool addEdge(const T &sourc, const T &dest, double w);

    bool removeEdge(const T &source, const T &dest);

    bool addBidirectionalEdge(const T &sourc, const T &dest, double w);

    virtual int getNumVertex() const;

    std::vector<Vertex<T> *> getVertexSet() const;

    std::vector<T> dfs() const;

    std::vector<T> dfs(const T &source) const;

    void dfsVisit(Vertex<T> *v, std::vector<T> &res) const;

    std::vector<T> bfs(const T &source) const;

    bool isDAG() const;

    bool dfsIsDAG(Vertex<T> *v) const;

    std::vector<T> topsort() const;

/*
 *  Adds a vertex with a given content or info (in) to a graph (this).
 *  Returns true if successful, and false if a vertex with that content already exists.
 */
    bool addVertex(const T &in);

protected:
    std::vector<Vertex<T> *> vertexSet;    // vertex set
    std::map<T, int> vertexIndex;    // vertex indexer (Vertex ID, Position on Vertex Set)

    double **distMatrix = nullptr;   // dist matrix for Floyd-Warshall
    int **pathMatrix = nullptr;   // path matrix for Floyd-Warshall

    int findVertexIdx(const T &in) const;
};

void deleteMatrix(int **m, int n);

void deleteMatrix(double **m, int n);


/************************* Vertex  **************************/

template<class T>
Vertex<T>::Vertex(T in) {
    this->info = in;
    this->longitude = NULL;
    this->latitude = NULL;
}

template<class T>
Vertex<T>::Vertex(T in, double longitude, double latitude) {
    this->info = in;
    this->longitude = longitude;
    this->latitude = latitude;
}

/*
 * Auxiliary function to add an outgoing edge to a vertex (this),
 * with a given destination vertex (d) and edge weight (w).
 */
template<class T>
Edge<T> *Vertex<T>::addEdge(Vertex<T> *d, double w) {
    auto newEdge = new Edge<T>(this, d, w);
    adj.push_back(newEdge);
    d->incoming.push_back(newEdge);
    return newEdge;
}

/*
 * Auxiliary function to remove an outgoing edge (with a given destination (d))
 * from a vertex (this).
 * Returns true if successful, and false if such edge does not exist.
 */
template<class T>
bool Vertex<T>::removeEdge(T in) {
    bool removedEdge = false;
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        Vertex<T> *dest = edge->getDest();
        if (dest->getInfo() == in) {
            it = adj.erase(it);
            deleteEdge(edge);
            removedEdge = true; // allows for multiple edges to connect the same pair of vertices (multigraph)
        } else {
            it++;
        }
    }
    return removedEdge;
}

/*
 * Auxiliary function to remove an outgoing edge of a vertex.
 */
template<class T>
void Vertex<T>::removeOutgoingEdges() {
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        it = adj.erase(it);
        deleteEdge(edge);
    }
}

template<class T>
bool Vertex<T>::operator<(Vertex<T> &vertex) const {
    return this->dist < vertex.dist;
}

template<class T>
T Vertex<T>::getInfo() const {
    return this->info;
}

template<class T>
std::vector<Edge<T> *> Vertex<T>::getAdj() const {
    return this->adj;
}

template<class T>
bool Vertex<T>::isVisited() const {
    return this->visited;
}

template<class T>
bool Vertex<T>::isProcessing() const {
    return this->processing;
}

template<class T>
unsigned int Vertex<T>::getIndegree() const {
    return this->indegree;
}

template<class T>
double Vertex<T>::getDist() const {
    return this->dist;
}

template<class T>
Edge<T> *Vertex<T>::getPath() const {
    return this->path;
}

template<class T>
std::vector<Edge<T> *> Vertex<T>::getIncoming() const {
    return this->incoming;
}

template<class T>
double Vertex<T>::getIncFlow() const {
    return this->incFlow;
}

template<class T>
void Vertex<T>::reverseAdj() {
    std::reverse(this->adj.begin(), this->adj.end());
}

template<class T>
void Vertex<T>::reverseIncoming() {
    std::reverse(this->incoming.begin(), this->incoming.end());
}

template<class T>
void Vertex<T>::setInfo(T in) {
    this->info = in;
}

template<class T>
void Vertex<T>::setVisited(bool visited) {
    this->visited = visited;
}

template<class T>
void Vertex<T>::setProcessing(bool processing) {
    this->processing = processing;
}

template<class T>
void Vertex<T>::setIndegree(unsigned int indegree) {
    this->indegree = indegree;
}

template<class T>
void Vertex<T>::setDist(double dist) {
    this->dist = dist;
}

template<class T>
void Vertex<T>::setPath(Edge<T> *path) {
    this->path = path;
}

template<class T>
void Vertex<T>::deleteEdge(Edge<T> *edge) {
    Vertex<T> *dest = edge->getDest();
    // Remove the corresponding edge from the incoming list
    auto it = dest->incoming.begin();
    while (it != dest->incoming.end()) {
        if ((*it)->getOrig()->getInfo() == info) {
            it = dest->incoming.erase(it);
        } else {
            it++;
        }
    }
    delete edge;
}

template<class T>
double Vertex<T>::getLongitude() const {
    return this->longitude;
}

template<class T>
double Vertex<T>::getLatitude() const {
    return this->latitude;
}


/********************** Edge  ****************************/

template<class T>
Edge<T>::Edge(Vertex<T> *orig, Vertex<T> *dest, double w): orig(orig), dest(dest), weight(w) {}

template<class T>
Vertex<T> *Edge<T>::getDest() const {
    return this->dest;
}

template<class T>
double Edge<T>::getWeight() const {
    return this->weight;
}

template<class T>
Vertex<T> *Edge<T>::getOrig() const {
    return this->orig;
}

template<class T>
Edge<T> *Edge<T>::getReverse() const {
    return this->reverse;
}

template<class T>
bool Edge<T>::isSelected() const {
    return this->selected;
}

template<class T>
double Edge<T>::getFlow() const {
    return flow;
}

template<class T>
double Edge<T>::getPreviousFlow() const {
    return previousFlow;
}

template<class T>
double Edge<T>::getTempWeight() const {
    return this->tempWeight;
}

template<class T>
void Edge<T>::setSelected(bool selected) {
    this->selected = selected;
}

template<class T>
void Edge<T>::setReverse(Edge<T> *reverse) {
    this->reverse = reverse;
}

template<class T>
void Edge<T>::setFlow(double flow) {
    this->flow = flow;
}

template<class T>
void Edge<T>::setPreviousFlow(double previousFlow) {
    this->previousFlow = previousFlow;
}

template<class T>
void Edge<T>::setWeight(double weight) {
    this->weight = weight;
}

template<class T>
void Edge<T>::setTempWeight(double tempWeight) {
    this->tempWeight = tempWeight;
}

template<class T>
void Vertex<T>::setIncFlow(double incFlow) {
    this->incFlow = incFlow;
}

/********************** Graph  ****************************/

template<class T>
int Graph<T>::getNumVertex() const {
    return vertexSet.size();
}

template<class T>
std::vector<Vertex<T> *> Graph<T>::getVertexSet() const {
    return vertexSet;
}

/*
 * Auxiliary function to find a vertex with a given content.
 */
template<class T>
Vertex<T> *Graph<T>::findVertex(const T &in) const {
    if (vertexIndex.count(in) == 0)
        return nullptr;

    auto i = vertexIndex.at(in);
    return vertexSet[i];
}

/*
 * Finds the index of the vertex with a given content.
 */
template<class T>
int Graph<T>::findVertexIdx(const T &in) const {
    for (unsigned i = 0; i < vertexSet.size(); i++)
        if (vertexSet[i]->getInfo() == in)
            return i;
    return -1;
}

/*
 *  Adds a vertex with a given content or info (in) to a graph (this).
 *  Returns true if successful, and false if a vertex with that content already exists.
 */
template<class T>
bool Graph<T>::addVertex(const T &in) {
    if (findVertex(in) != nullptr)
        return false;
    vertexSet.push_back(new Vertex<T>(in));
    vertexIndex[in] = vertexSet.size() - 1;
    return true;
}

/*
 *  Removes a vertex with a given content (in) from a graph (this), and
 *  all outgoing and incoming edges.
 *  Returns true if successful, and false if such vertex does not exist.
 */
template<class T>
bool Graph<T>::removeVertex(const T &in) {
    if (vertexIndex.count(in) == 0)
        return false;

    auto i = vertexIndex[in];
    auto v = vertexSet[i];
    v->removeOutgoingEdges();
    for (auto u: vertexSet) {
        u->removeEdge(v->getInfo());
    }
    vertexSet.erase(&vertexSet[i]);
    vertexIndex.erase(in);
    delete v;
    return true;
}

/*
 * Adds an edge to a graph (this), given the contents of the source and
 * destination vertices and the edge weight (w).
 * Returns true if successful, and false if the source or destination vertex does not exist.
 */
template<class T>
bool Graph<T>::addEdge(const T &sourc, const T &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    v1->addEdge(v2, w);
    return true;
}

/*
 * Removes an edge from a graph (this).
 * The edge is identified by the source (sourc) and destination (dest) contents.
 * Returns true if successful, and false if such edge does not exist.
 */
template<class T>
bool Graph<T>::removeEdge(const T &sourc, const T &dest) {
    Vertex<T> *srcVertex = findVertex(sourc);
    if (srcVertex == nullptr) {
        return false;
    }
    return srcVertex->removeEdge(dest);
}

template<class T>
bool Graph<T>::addBidirectionalEdge(const T &sourc, const T &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    auto e1 = v1->addEdge(v2, w);
    auto e2 = v2->addEdge(v1, w);
    e1->setReverse(e2);
    e2->setReverse(e1);
    return true;
}

/****************** DFS ********************/

/*
 * Performs a depth-first search (dfs) traversal in a graph (this).
 * Returns a vector with the contents of the vertices by dfs order.
 */
template<class T>
std::vector<T> Graph<T>::dfs() const {
    std::vector<T> res;
    for (auto v: vertexSet)
        v->setVisited(false);
    for (auto v: vertexSet)
        if (!v->isVisited())
            dfsVisit(v, res);
    return res;
}

/*
 * Performs a depth-first search (dfs) in a graph (this) from the source node.
 * Returns a vector with the contents of the vertices by dfs order.
 */
template<class T>
std::vector<T> Graph<T>::dfs(const T &source) const {
    std::vector<int> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }
    // Set that no vertex has been visited yet
    for (auto v: vertexSet) {
        v->setVisited(false);
    }
    // Perform the actual DFS using recursion
    dfsVisit(s, res);

    return res;
}

/*
 * Auxiliary function that visits a vertex (v) and its adjacent, recursively.
 * Updates a parameter with the list of visited node contents.
 */
template<class T>
void Graph<T>::dfsVisit(Vertex<T> *v, std::vector<T> &res) const {
    v->setVisited(true);
    res.push_back(v->getInfo());
    for (auto &e: v->getAdj()) {
        auto w = e->getDest();
        if (!w->isVisited()) {
            dfsVisit(w, res);
        }
    }
}

/****************** BFS ********************/
/*
 * Performs a breadth-first search (bfs) in a graph (this), starting
 * from the vertex with the given source contents (source).
 * Returns a vector with the contents of the vertices by bfs order.
 */
template<class T>
std::vector<T> Graph<T>::bfs(const T &source) const {
    std::vector<int> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }

    // Set that no vertex has been visited yet
    for (auto v: vertexSet) {
        v->setVisited(false);
    }

    // Perform the actual BFS using a queue
    std::queue<Vertex<T> *> q;
    q.push(s);
    s->setVisited(true);
    while (!q.empty()) {
        auto v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for (auto &e: v->getAdj()) {
            auto w = e->getDest();
            if (!w->isVisited()) {
                q.push(w);
                w->setVisited(true);
            }
        }
    }
    return res;
}

/****************** isDAG  ********************/
/*
 * Performs a depth-first search in a graph (this), to determine if the graph
 * is acyclic (acyclic directed graph or DAG).
 * During the search, a cycle is found if an edge connects to a vertex
 * that is being processed in the stack of recursive calls (see theoretical classes).
 * Returns true if the graph is acyclic, and false otherwise.
 */

template<class T>
bool Graph<T>::isDAG() const {
    for (auto v: vertexSet) {
        v->setVisited(false);
        v->setProcesssing(false);
    }
    for (auto v: vertexSet) {
        if (!v->isVisited()) {
            if (!dfsIsDAG(v)) return false;
        }
    }
    return true;
}

/**
 * Auxiliary function that visits a vertex (v) and its adjacent, recursively.
 * Returns false (not acyclic) if an edge to a vertex in the stack is found.
 */
template<class T>
bool Graph<T>::dfsIsDAG(Vertex<T> *v) const {
    v->setVisited(true);
    v->setProcesssing(true);
    for (auto e: v->getAdj()) {
        auto w = e->getDest();
        if (w->isProcessing()) return false;
        if (!w->isVisited()) {
            if (!dfsIsDAG(w)) return false;
        }
    }
    v->setProcesssing(false);
    return true;
}

/****************** toposort ********************/
//=============================================================================
// Exercise 1: Topological Sorting
//=============================================================================
// TODO
/*
 * Performs a topological sorting of the vertices of a graph (this).
 * Returns a vector with the contents of the vertices by topological order.
 * If the graph has cycles, returns an empty vector.
 * Follows the algorithm described in theoretical classes.
 */

template<class T>
std::vector<T> Graph<T>::topsort() const {
    std::vector<int> res;

    for (auto v: vertexSet) {
        v->setIndegree(0);
    }
    for (auto v: vertexSet) {
        for (auto e: v->getAdj()) {
            unsigned int indegree = e->getDest()->getIndegree();
            e->getDest()->setIndegree(indegree + 1);
        }
    }

    std::queue<Vertex<T> *> q;
    for (auto v: vertexSet) {
        if (v->getIndegree() == 0) {
            q.push(v);
        }
    }

    while (!q.empty()) {
        Vertex<T> *v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for (auto e: v->getAdj()) {
            auto w = e->getDest();
            w->setIndegree(w->getIndegree() - 1);
            if (w->getIndegree() == 0) {
                q.push(w);
            }
        }
    }

    if (res.size() != vertexSet.size()) {
        //std::cout << "Impossible topological ordering!" << std::endl;
        res.clear();
        return res;
    }

    return res;
}

inline void deleteMatrix(int **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete[] m[i];
        delete[] m;
    }
}

inline void deleteMatrix(double **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete[] m[i];
        delete[] m;
    }
}

template<class T>
Graph<T>::~Graph() {
    deleteMatrix(distMatrix, vertexSet.size());
    deleteMatrix(pathMatrix, vertexSet.size());
}

#endif /* DA_TP_CLASSES_GRAPH */