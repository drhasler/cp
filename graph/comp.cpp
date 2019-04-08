/* git@drhasler - 03/19 */
#include <bits/stdc++.h>
using namespace std;

#define vi vector<int>
#define pi pair<int,int>
#define xx first
#define yy second
#define foro(a,b) for (int a=0;a<(b);a++)
#define rofo(a,b) for (int a=(b)-1;a>=0;a--)

/* Components
 * summary:
 * - disjoint set union
 * - topological sort
 * - bridges, art pt, biconnected components
 * - strongly connected components
 * - lowest common ancestor */

const int N; // upper bound for number of nodes
bool seen[N];

/* disjoint set union - DSU
 * O(n α(n))
 * dsu[x] is parent or -size if self par */

int dsu[N];

void init(int n) { fill(dsu,dsu+n,-1); }

int find(int x) { return dsu[x]<0?x:dsu[x]=find(dsu[x]); }

void merge(int u, int v) {
    u = find(u); v = find(v);
    if (u==v) return;
    // biggest is most negative
    if (dsu[v]<dsu[u]) swap(u,v);
    dsu[u] += dsu[v];
    dsu[v] = u;
}


/* topological sort
 * O(n+m)
 * doesnt check for backedges
 * patial sort or check if total order */

int top[N];

void topo_aux(int v, int &idx) {
    used[v] = 1;
    // pass parent for cycle check
    for (int u:e[v]) if (!used[u]) topo_aux(u,idx);
    top[--idx] = v;
}

void toposort(int n) {
    int idx = n;
    foro(i,n) if (!used[i]) topo_aux(i,idx);
}


/* bridges - articulation points - biconnected components
 * O(n+m) offline
 * undirected graph
 * biconnected components are separated by bridges
 * using a map, we can delete bridges from the original graph
 * and run traversal to get biconnected components */

namespace OfflineBridges {
    int timer;
    int tin[N]; // time of entry
    int low[N]; // lowest time of entry for this component
    vi artPts;
    vector<pi> bridges;

    void DFS(int v,int p) { // p = v for root
        seen[v] = 1;
        tin[i] = low[i] = timer++;
        int cnt = 0; // cnt of children comp (for artPts)
        bool isArt = 0;
        for (int u:e[v]) {
            if (!seen[u]) {
                dfs(u,v);
                low[v] = min(low[v],low[u]); // unroll
                if (low[u]>tin[v]) bridges.pb({v,u}); // new bcc
                ++cnt;
                if (low[u]>=tin[v] && p!=v) isArt = 1;
            } else if (v != p) { // back edge
                low[v] = min(low[v],tin[u]); // unroll
            }
        }
        // is linked to more than one component
        if (p == v ? cnt > 1 : isArt) artPts.pb(v);
    }

    int BCC(int n) { // number of biconnected components
        int cnt = 0; // number of connected components
        fill(seen,seen+n,0);
        FOR (i,n) if (!seen[i]) {
            ++cnt;
            DFS(i,i);
        }
        return cnt + bridges.size();
    }
}

/* bridges online
 * amortized O(logN) per edge
 * using LCA */

namespace OnlineBridges {
    int bridges; // we could also be store them in a map
    int it; // for naive LCA
    int lit[N]; // last it
    int par[N]; // parent in block-cut forest
    int dsu_cc[N]; // bcc to cc
    int dsu_bcc[N]; // vertex to cc
    // usual find and merge

    void init(int n) {
        fill(par,par+n,-1);
        iota(dsu_cc,dsu_cc+n,0);
        iota(dsu_bcc,dsu_bcc+n,0);
    }

    void reroot(int v) { // makes parent a child
        int p = par[v];
        if (p<0) return;
        reroot(p);
        par[p] = v;
    }

    void compress_loop(int a,int b) { // amortized O(1)
        int lca;
        ++it;
        stack<int> sa,sb;
        while (1) {
            if (a>=0) {
                sa.push(a);
                if (lit[a]==it) { lca = a; break; }
                lit[a] = it;
                a = par[a];
            }
            if (b>=0) {
                sb.push(b);
                if (lit[b]==it) { lca = b; break; }
                lit[b] = it;
                b = par[b];
            }
        }
        while (sa.top()!=lca) sa.pop();
        sa.pop();
        while (!sa.empty()) {
            --bridges;
            merge_bcc(lca,sa.top());
            sa.pop();
        }
        while (sb.top()!=lca) sb.pop();
        sb.pop();
        while (!sb.empty()) {
            --bridges;
            merge_bcc(lca,sb.top());
            sb.pop();
        }
    }

    void add_edge(int a,int b) {
        a = find_bcc(a);
        b = find_bcc(b);
        if (a == b) return;
        // now we consider two block-cut trees vertices
        // get their root
        int ca = find_cc(a);
        int cb = find_cc(b);

        if (ca != cb) { // merge trees
            ++bridges;
            // we will reroot the smaller tree
            if (dsu_cc_size[ca] > dsu_cc_size[cb]) {
                swap(a, b);
                swap(ca, cb);
            }
            reroot(a);
            par[a] = dsu_cc[a] = b;
            dsu_cc_size[cb] += dsu_cc_size[a];
        } else {
            compress_loop(a,b);
        }
    }
}


/* Strongly Connected Components
 * O(n+m)
 * directed graph w/ partial topo sort */

vector<int> scc[N];

void SCC_aux(int v,int c) {
    scc[c].pb(v);
    used[v] = 1;
    for (int u:e[v]) if (!used[u]) SCC_aux(u,c);
}

void build_SCC(int n) {
    toposort(n);
    fill(used,used+n,0);
    int cnt=0;
    rofo (i,n) {
        int v = top[i];
        if (!used[v]) SCC_aux(v,cnt++);
    }
}


/* online lowest common ancestor - LCA
 * O(logN) query time
 * N power of 2 */

int T[4*N]; // segment tree of euler tour of idx
int toT[N]; // vertex to pos in T
int fromT[N]; // from value in T to vertex
int T_idx = 2*N; // index in tree
int e_idx; // index in euler tour

void euler_dfs(int v,int p) {
    int vi = e_idx++;
    toT[v] = T_idx;
    fromT[vi] = v;
    T[T_idx++] = vi;
    for (int u:e[v]) if (u!=p) {
        LCA_aux(u,v);
        T[T_idx++] = vi; // eulers return
    }
}

void LCA_build() {
    LCA_aux(1,1); // explore tree from root
    ROF (i,2*N) T[i] = min(T[2*i],T[2*i+1]); // T[0] doesnt matter
}

int LCA(int u,int v) {
    // vertices to positions
    u = toT[u]; v = toT[v];
    if (u>v) swap(u,v);
    v++; // exclusive range mininum query
    int ans = 1e9; // inf
    while (u<v) {
        if (u&1) ans = min(ans,T[u++]);
        if (v&1) ans = min(ans,T[--v]);
        u/=2; v/=2;
    }
    return fromT[ans]; // real vertex
}


/* Offline Lowest Common Ancestor
 * O(q α(n) + n)
 * requires dsu init */
namespace OLCA {

    int anc[N]; // anc[find(v)] := ancestor of the set v belongs to
    int ans[Q]; // ans[q_idx]:= answer
    vector<pi> queries[N]; // pair<other,q_idx> - for both v

    void OLCA(int v,int p) { // call LCA(root,root)
        anc[v] = v; // or iota(anc,anc+n,0);
        for (int u:e[v]) if (u!=p) {
            LCA(u,v); // post-order
            merge(u,v);
            anc[find(v)] = v;
        }
        seen[v] = 1;
        for (pii q:queries[v]) {
            if (seen[q.xx]) {
                ans[q.yy] = anc[find(q.xx)];
            }
        }
    }
}

