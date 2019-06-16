**COMP3927**
- [Asymptotic analysis](#asymptotic-analysis)
- [Graphs](#graphs)
- [Greedy](#greedy)
- [Dynamic programming](#dynamic-programming)
- [Network flow](#network-flow)
- [NP](#np)
- [prove templete](#prove-templete)

# Asymptotic analysis
## symbols
&theta;(n), O(n), &Omega;(n)
## exp vs poly
exponential time: c<sup>n</sup>.  
polynomial time: a*n<sup>b</sup>. (log is also in poly time)   
a,b,c are constant greater than 1  

factorial time > exponential time > polynomial time  
(n! > c<sup>n</sup> > n<sup>c</sup>)

# Graphs
## basics
### representation
  - adiacency matrix  
  - adjacency list  
  
### path
a path in an undirected graph G = (V,E)is a sequence of vertices that each consecutive pair is joined by an edge in E  
path is simple if all vertices are distinct.  
### connectness
an undirected graph is connected if for every pair of vertices, there is path between them.  
### cycle
A cycle is a path v1,v2...,vk in which  
  1. v1 = vk
  2. k > 2
  3. the first k-1 vertices are all distinct.  
  
**tree**: an undirected graph that is connected and does not contain a cycle.  
In a tree, there is n-1 edges.
## transitive closure
### definition
The transitive closure of a graph G = (V,E) is a new graph G’ that has the same set of V, but there is an edge between every connect pair of vertices in G.
### properties
the transitive closure of a connected graph is a complete graph.  
### algorithm
Use BFS to find all connected component in a graph, make them complete.  
## BFS
### applications
  1. calculate the shortest path length from a single *s*.
  2. find the set of vertices connected to a single source *s*.
  3. find transitive closure.
  
### time
  O(n+m) *considering all egdes and vertices.*
### space
  - adjacency list: &theta;(|V|+|E|).  
  - adjacency matrix: &theta;(|V|<sup>2</sup>).
  
## DFS
### application
  1. check connectness. If the DFS result is a single tree, then the graph is connected.
  2. find all paths from a vertex to any others.

### time
  O(n+m) *considering all egdes and vertices.*
## bipartite graph
### definition
The nodes in a graph can be colored red or blue such that every edge has one red and one blue end.
### application
  1. stable marriage
  2. job scheduling
  
### test bipartiteness
A graph G is bipartite if and only if it contains no odd length cycle.
#### algorithm to test bipartiteness
run BFS start from an arbitrary vertex s, if no edge of G joins two vertices of the same layer, G is bipartite. *n edge of G joins two nodes of the same layer and G contains an odd-length cycle, hence G is not bipartite.*  

time: O(n+m) *BFS*
## cut edge
### definition
the removal of a cut edge disconnects the graph.
### algorithm
Use DFS to construct a DFS tree, iteratively remove an edge in the DFS tree and check if the graph is still connected (using another DFS).

time: O(nm) *DFS+n*DFS*
# Greedy
The best way to say greedy won't work is giving ONE counter example.
## applications
### interval scheduling
find maximum subset of mutually compatible jobs.
#### greedy approach
sorted by earliest finish time.
#### prove
Exchange argument
### interval partitioning
the minimum number of classrooms to schedule all lectures so that no two occur at the same time in the same room.  

**depth**: the max number of intervals in any given time.  
Number of classrooms needed >= depth.

#### greedy approach
sorted by starting time
#### prove
greedy stay ahead
### Scheduling to minimizing lateness
give a list of jobs (process time and due time) to a person, schedule all jobs to minimize maximum lateness.  

**idle time**: the time that the person isn’t doing any job.  
There exists an optimal schedule with no idle time. (greedy)  

**inversion**: for a pair of jobs i and k, the due time for i < k but k is scheduled before i.

#### greedy approach
sorted by earlist due time.  
  - greed yschedule has no inversions.
  - Swap two adjacent, inverted jobs reduces the number of inversions by one and does not increase the max lateness.
## Minimum Spanning Tree
### definition
Given a connected graph G with edge weights, a MST is a subset of the edges such that the subset is a spanning tree whose sum of edge weights is minimized.  
such that, a tree is a MST then:  
  1. connected
  2. no cycle
  3. minimum sum of edge weights

**cayley's theorem**: there are n<sup>n-2</sup> spanning trees of K<sub>n</sub>.  
**cycle**: set of edges that form a cycle  
**cut**: a subset of vertices
**cutset**: the subset of edges in a cut with exactly one endpoint in the cut.
### properties
  - **simplifying assumption**: All edge costs are distinct.
  - **cut property**:  Let S be any subset of nodes, and let e be the min cost edge with exactly one endpoint in S, then MST contains e.
  - **cycle property**: Let C be any cycle, and let f be the max cost edge belonging to C. then the MST does not contain f.
### prim
#### process
  1. initialize S = any node
  2. Apply cut property
  3. Add min cost edge in cutset corresponding to S, add the new vertex to S
#### properties
  - If edge costs are distinct, then MST is unique.
  - At every step of Prim’s, current tree is contained in every MST.
#### time complexity
  - use priority queue: O(mlogn)
  - use array: O(n<sup>2</sup>)
### kruskal
#### process
  1. start with empty graph T, consider edges in ascending order of weight
  2. If adding e to T creates a cycle, discard e according to cycle property.
  3. Otherwise, insert e = (u, v) into T according to cut property where S = set of nodes in u's connected component.


**Lexicographic tirebreaking**: if the simplifying assumption property cannot be guaranteed, the edge with smaller number is less.

### shortest paths
in a shortest path problem, we are given:  
  - directed graph G = (V,E)
  - source vertex s, destination vertex t
  - length l<sub>e</sub> = length of edge e

the target is to find the shortest directed path from s to t.  
if the length of all egdes are unit, can be solved using **BFS**.
#### algorithm
Can be computed using Dijkstra.  
time complexity is O(m+nlogn) using priority queue.

# Divide and conquer
  1. **Divide**: Break up problem into several parts.
  2. **Conquer**: Solve each part recursively.
  3. **Combine**: solutions to sub-problems into overall solution.
## application
general application:
  - binary search: time O(logn)
  - merge sort: time O(nlogn)
  
complex application:
  - counting inversions: time O(nlogn)
  - closest pair of points: preprocessing time O(nlogn) *sort* + running time O(nlogn).
  - multiplication: O(n<sup>2</sup>)
  
## unrolling
### geometric series
  - 1 + x + x<sup>2</sup> + ... = 1/(1-x) for x < 1
  - 1 + x + x<sup>2</sup> + ... + x<sup>n</sup> = (1-x<sup>n+1</sup>)/(1-x) for x != 1
  
### process
draw the recursion tree in terms of non recursion time in every level. and sum them all together.

## master theorem
for a recursion time expression:  
T(n) = a * T(n/b) + n<sup>d</sup> * log<sup>k</sup>n, where a>0, b>0, d>=0 and k>=0.  
f(n) = n<sup>d</sup> * log<sup>k</sup>n
### case 1: d < log<sub>b</sub>a
T(n) = &theta;(n<sup>log<sub>b</sub>a</sup>)
### case 2: d = log<sub>b</sub>a
T(n) = &theta;(n<sup>log<sub>b</sub>a</sup>log<sup>k+1</sup>n)
### case 3: d > log<sub>b</sub>a
And satisfy the **regularity condition** that a * f(n/b) <= c * f(n) for some constant c < 1  
T(n) = &theta;(f(n))

# Dynamic programming
divide and conquer breaks a problem into sub problems and solve each of them **independently**.  
dynamic programming breaks a problem into a series of **overlapping** sub problems.
## key steps of DP
  1. define subproblems
  2. find recurrence relating subproblems
  3. solve the base cases.
  
## approaches
  1. top down: start from the biggest problem. If a smaller problem hasn't been found, go find it, otherwise use it directly. The technique is called **memoization**. top down is usually implemented in recursion.
  2. bottom up: start from the smallest problem. Construct all sequential problems using the small problems that have been solved. bottom up is usually implemented in loops, and its easier to use.
  
## application 
  - 1D dp
    - weighted interval scheduling: preprocessing O(nlogn) + running O(n). *unweighted version can be solved by greedy*
    - segmented least squares: time O(n<sup>3</sup>), space O(n<sup>2</sup>)
    - maximum-sum contiguous subarray (MCS): time O(n)
    - longest increasing subsequence: time O(nlogn)
  - 2D dp
    - knapsack: time O(nW)
    - shortest path with negative edges: **Bellman-Ford** algo, time O(mn), space O(n)
    - longest common subsequence: time O(nm), space O(nm)
  - dp over intervals
    - RNA secondary stucture: time O(n<sup>3</sup>)
    
# Network flow
## basics
### definition
  – Abstraction for material flowing through the edges.  
  – G = (V, E): a directed graph with no parallel edges.  
  – Two distinguished nodes: s = source, t = sink.  
  – The source has no incoming edges and the sink has no outgoing edges.  
  – c(e) = capacity of edge e.  
  
### s-t flow
A s-t flow satisfies:
  - for each edge e, 0 <= f(e) <= c(e)   *capacity*
  - e is saturated if f(e) = c(e)
  - for each vertex v, the flow goes into v equals the flow leaves v.  *conservation*
  
**value of a flow**: is calculated by summing all the flow leaves the source s.
## max flow problem
### definition
A flow is max flow if:
  - the flow saturates every edge OR  
  - the flow saturates every edge out of s
  
A **s-t cut** is a partition (A,B) of V such that s is in A and t is in B    
The **capacity** of a (A,B) cut, cap(A,B) is the sum of the capacities of all edges come out of A.  
**min s-t cut problem**: find a s-t cut of minimum capacity. The min cut problem leads to max flow problem as **min cut = max flow**  
**flow value lemma**: Let f be any flow, and let (A, B) be any s-t cut. Then, the **net flow** sent across the cut is equal to the amount
leaving s. such that v(f) = f<sup>out</sup>(A) - f<sup>in</sup>(A)  
**weak duality**: Let f be any flow, and let (A, B) be any s-t cut. Then the value of the flow is at most the capacity of the cut, such that v(f) <= cap(A,B)  
## process
use a greedy approach:  
  – Start with f(e) = 0, for all edge e in E  
  – Find an s-t path P where each edge has f(e) < c(e).  
  – **Augment flow along path P.**  
  – Repeat until you get stuck.  
  
## residual graph
original edge has capacity c(e) and flow f(e).  
the residual edge:
  - "undo" flow sent
  - e = (u,v) and e' = (v,u)
  - calculate residual capacity c<sub>f</sub>(e):
    - = c(e) - f(e) if e in E
    - = f(e) if e' in E

residual graph G' is the graph consists of residual edges.
## augment path algorithm
P = a simple s-t path in G'  
bottleneck(P,f) = minimum residual capacity of any edge on P with respect to the current flow f.  
minus the bottleneck value to all forward edges in P, add the bottleneck value to all backward edges in P.  
**augmenting path theorem**: Flow f is a max flow if and only if there are no augmenting paths in the residual graph.  

**Max-flow min-cut theorem**: The value of the max flow is equal to the value of the min cut.  

**Integrality theorem**: If all capacities are integers then every flow value f(e) and every residual capacities c'(e) remains an integer throughout the algorithm.  
## max flow algorithm
  1. Ford-Fulkerson: time O(mF)  
  2. Ford-Fulkerson with max bottleneck augmenting: time O(m<sup>2</sup>F)  
  3. Ford-Fulkerson with shortest augmenting path: time O(nm<sup>2</sup>)  
  
If F is extremely large choose 2, otherwise 1 with a general augmenting path strategy runs faster.
## reduce problem to max flow problem
There are following steps:
  1. **max flow formulation**. remember to explicitly state:
    - the directed graph with vertices set and edge set
    - edge capacities
    - source s and sink t
    - the edges come out of s and the edges come to t.
  2. **use a max flow algo**. there are [three choices](#max-flow-algorithm), choose the one runs fast.
  3. **extract solution from max flow**. explicitly state the property relationship between both problems.
  
Pove the correctness of the reduction from both directions, such that:
  - max flow value -> the problem solution  
  (integrality theorem may help if the capacities are all unit) and
  - the problem solution -> max flow value
  
## circulation with supply and demand
  - directed graph G = (V,E)
  - edge capacities c(e) for e in E
  - node supply and demands d(v) for v in V
    - demand if d(v) > 0; supply if d(v) < 0; transshipment if d(v) = 0
    
**circulation**: A circulation is a function that satisfies:
  - for each e in E, 0 <= f(e) <= c(e)  *capacity*
  - for each v in V, the sum of flow come in v - the sum of flow come out v = d(v)  *conservation*
  
### circulation problem
given (V,E,c,d), does there exist a circulation?  
**necessary condition**: sum of supplies = sum of demands, otherwise there is no circulation.  

solve the circulation by formulating to max flow:
  - add new source s and sink t
  – For each v with d(v) < 0, add edge (s, v) with capacity -d(v).
  – For each v with d(v) > 0, add edge (v, t) with capacity d(v).
  
**The graph has circulation if and only if max flow of the formulated network flow saturates all edges leaving s and entering t**  
## survey design
  – Design survey asking n1 consumers about n2 products.  
  – Can only survey consumer i about a product j if they own it.  
  – Ask consumer i between ci and ci' questions.  
  – Ask between pj and pj' consumers about product j.  
  
# NP
## polynomial reduction (Cook)
problem X **polynomial reduces to** problem Y, denoted X &leq;<sub>p</sub> Y, if arbitrary instances of problem X can be solved using:  
  – Polynomial number of standard computational steps, plus  
  – **Polynomial number** of calls to an oracle/black box that solves problem Y.  
  
## purpose
  1. Design algorithms: If X &leq;<sub>p</sub> Y and Y can be solved in polynomial-time, then X **can** also be solved in polynomial time.  
  2. Establish intractability: If X &leq;<sub>p</sub> Y and X cannot be solved in polynomial-time, then Y **cannot** be solved in polynomial time.
  
**transitivity**: if X &leq;<sub>p</sub> Y and Y &leq;<sub>p</sub> Z, then X &leq;<sub>p</sub> Z.  
**equivalence**: if X &leq;<sub>p</sub> Y and Y &leq;<sub>p</sub> X, then X &equiv;<sub>p</sub> Y.

## decision vs optimization vs search
### decision problem
Does there **exist** an object satisfying some given properties?  

### optimization problem
What is the **size** of the biggest/smallest such object?  

### search problem
**Find** a biggest/smallest such object  

**decision &equiv;<sub>p</sub> optimization &equiv;<sub>p</sub> search**
## reduction process
there are such steps to reduce problem X to Y:
  1. Construct a polynomial-time algorithm that transforms instance I of X to an instance f(I) of Y.  
  2. Prove correctness, which boils down to showing answer for I is Yes iff answer for f(I) is Yes.  
  
## P and NP and NP-complete
### P
Decision problems for which there is a deterministic poly-time algorithm.
### NP
Decision problems for which YES-instances have a poly-time certifier.  

**certificate**: a solution instance of the problem.  

**certifier**: an algorithm to check the correctness of a certificate, that the result is yes if and only if the the certificate is a valid solution.  

#### prove NP
  - give an arbitrary certificate  
  - state a certifier that checks the certificate is 'yes' or 'no'.  
  - prove the certifier is in poly-time.  
  
#### NP vs P
P is a subset of NP, but it is not sure P = NP or P != NP.
### NP-complete
#### polynomial transformation (Karp)
problem X **polynomial transforms to** problem Y, denoted X &leq;<sub>p</sub> Y, if arbitrary instances of problem X can be solved using:  
  – Polynomial number of standard computational steps, plus  
  – **one** call to an oracle/black box that solves problem Y.  
  
#### NP-hardness
There are a set of [NP-complete problems](#np-complete-problems-need-to-know) out there, pick one NP-complete problem X. For a problem Y, we can prove Y's NP-hardness by X &leq;<sub>p</sub>Y (use Karp reduction).  

**NP-hard problem is not necessary in NP.**
#### prove NP-completeness
if we want to prove problem Y is NP-complete.
  1. [prove Y is in NP](#prove-np)
  2. [prove Y is NP-hard](#np-hardness) use karp reduction
  
#### NP-complete problems need to know
  - 3-SAT: Given a formula in Conjunctive Normal Form, where each clause has exactly 3 literals, is there a satisfying truth assignment?  
  - 3-Colour: Given a graph, is it possible to colour the vertices using at most 3 colours, such that no pair of adjacent vertices have the same colour?  
  - Hamiltonian-Cycle: Given a graph, is there a simple cycle which visits every vertex of the graph?  
  - Independent-Set: Given a graph G = (V,E) and an integer k, is there a subset of vertices S ⊆ V such that |S| ≥ k and for each edge at most one of its endpoints is in S?  
  - Knapsack: Given a set of items with weights and values, an integer capacity C, and an integer k, is it possible to select a subset of items with total weight no greater than C, but total value at least k?  
  - Subset-Sum: Given a set of integers S and an integer k, is there a subset of S which sums to k?  
  - Partition: Given a set of integers S, is there a partition of S into two subsets S1 and S2 such that the sum of integers in S1 equals the sum of integers in S2?  
  - Set-Cover: Given a set of elements U, a collection of subsets of it, and an integer k, does there exist a collection of k of these subsets whose union covers U?  
  - Traveling-Salesman-Problem: Given a set of n cities and a pairwise distance function d(u,v), is there a tour of length ≤ D?  
  - Vertex-Cover: Given a graph G = (V,E) and integer k, is there a subset of vertices S ⊆ V such that |S| ≤ k and for each edge at least one of its endpoints is in S?  

# prove templete
## exchange argument (for greedy)
  1. Define your greedy solution.  
  2. Compare solutions. If X<sub>greedy</sub> !=; X<sub>opt</sub>, then they must differ in some specific way.  
  3. Exchange Pieces. Transform X<sub>opt</sub> to a solution that is “closer” to X<sub>greedy</sub> and prove cost doesn’t increase.  
  4. Iterate. By iteratively exchanging pieces one can turn X<sub>opt</sub> into X<sub>greedy</sub> without impacting the quality of the solution.  
### exchange argument templete
  - Assume the solution produced by the greedy algorithm is g1, g2,... and there exists an optimal solution i1, i2,...  
  - If they are identical, then the prove is done. Otherwise, they must differ in some specific way.  
  - consider base case (when the idx = 1). If g1 is not i1, and by the definition of geedy, g1 ... than i1, which means g1 will not contradict i2,i3,... Therefore, we can swap i1 for g1 and the solution is still feasible and optimal.  
  - By inductive hypothesis, assume 1->idx-1 g1 = i1, g2 = i2... are all optimal.  
  - consider inductive case (when idx > 1). if g idx is not i idx, and by the definition of greedy, g idx ... than i idx, which means g idx will not contradict i idx, i idx+1... Therefore, we can swap i idx for g idx and the solution is still feasible and optimal  
  - Continuely iterating the process until the optimal solution is entirelly transformed into the greedy solution wihout increasing any potential cost.  
  
## induction
  1. Define induction invariant n. such that, the variable for induction.
  2. prove the **base case** n = 1 is correct.
  3. explicitly state **inductive hypothesis**: suppose the theorem holds for all instances of size < n. *n must be integer*!
  4. prove the **inductive step**: such that give 2 and 3, for n > 1, the instance holds the theorem.
  5. Therefore, we can deduce by induction that the theorem holds for all instance for size >= 1.
  
## max flow reduction
  1. problem -> according flow value:  
    - satisfy capacity constraint  
    - flow value add up equals to the problem solution
   
  2. according flow value -> problem:  
    - use integrality theorem, to say the flow value of each edge is integral  
    - consider the specific arrangement
    
  3. In summary, given any arrangement of ..., we can find a flow of the same size and vice versa. Therefore the optimal solution of the problem is equivalent to the maximum flow of the formulated flow network.
  
## prove NP-hardness
  1. Yes for the NP hard problem -> Yes for problem
  2. Yes for problem -> Yes for the NP hard problem

  
