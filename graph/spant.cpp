/* git@drhasler - 04/19 */
#include <bits/stdc++.h>
using namespace std;

#define pq priority_queue
#define pi pair<int,int>
#define xx first
#define yy second
#define foro(a,b) for (int a=0;a<(b);a++)

/* Minimum Spanning Tree
 * of undirected graph
 * two popular methods:
 * - Prim: extend tree
 *   by adding the cheapest edge
 * - Kruskal: merge trees */

vector<pi> T[N]; // MST
struct edge { int s,d,w; }; // src dest weight

/* Prim's algorithm for dense graphs
 * O(n2) */
namespace PrimDense {
    int adj[N][N]; // weights (INF if not connected)
    bool used[N];
    pii e[N]; // edge[dst] := {src, min cst}
    int MST(int n) {
        fill(e,e+n,{0,INF});
        e[0] = {0,0}; // first we add 0
        int sum = 0;
        foro(i,n) { // we add n nodes
            int v = -1; // cheapest neighbor
            FOR(i,n) if (!used[i] && (v < 0 || e[i].yy < e[v].yy)) v = i;
            // disconnected graph, no mst
            if (e[v].yy == INF) return -1;
            // connecc
            used[v] = 1; sum += e[v].yy;
            // if not root add edge
            if (v) T[e[v].xx].pb({v,e[v].yy});
            // add neighbors
            FOR(i,n) if (adj[v][i]<e[i].yy) e[i] = {v,adj[v][i]};
        }
        return sum;
    }
}

/* Prim's aglorithm for sparse graphs
 * O(mlogm) */
namespace PrimSparse {
    bool used[N];
    vector<pi> e[N]; // edge[src] = {{dst,cst}..}
    int MST(int n) {
        int sum = 0;
        fill(used,used+n,0);
        pq<edge, vector<edge>,
            [](const edge& a, const edge& b) { return a.w>b.w; }> q;
        q.push({0,0,0});
        
        while (!q.empty()) {
            // cheapest neighbor
            edge E = q.top(); q.pop();
            if (used[E.d]) continue;
            used[E.d] = 1; sum += E.w;
            // if not root add edge
            if (E.d) T[E.s].pb({E.d,E.w});
            // add neighbors
            for (pii x:e[E.d]) q.push({E.d,x.xx,x.yy});
        }
        return sum;
    }
}


/* Kruskal's algorithm
 * O(MlogM) using DSU */

namespace Kruskal {
    int dsu[N];
    int find(int v) { return dsu[v]<0 ? v : dsu[v] = find(dsu[v]); }
    void merge(int a, int b) {
        a = find(a); b = find(b);
        if (a==b) return;
        if (dsu[a]>dsu[b]) swap(a,b);
        dsu[a] += dsu[b];
        dsu[b] = a;
    }
    vector<edge> edg;
    int MST(int n) {
        fill(dsu,dsu+n,-1);
        int sum = 0;
        sort(edg.begin(),edg.end(),
            [](const edge& a, const edge& b){ return a.w<b.w; });
        for (edge e:edg) {
            if (find(e.s)==find(e.d)) continue;
            sum += e.w;
            T[e.s].pb({e.d,e.w});
            merge(e.s,e.d);
        }
        return sum;
    }
}

/* Kirchhoff's theorem
 * for counting number of spanning trees
 * - construct Laplacian matrix of graph Lii = deg(i) Lij = adj(i,j) ? -1 : 0
 * - compute any cofactor ie det of matrix w/o 1 row & 1 col
 * proof: TODO */

int main() {
    if (0) {
        using namespace DensePrim;
        int n; cin >> n;
        foro(i,n) foro(j,n) cin >> adj[i][j];
        int cost = MST(n);
    }
}
