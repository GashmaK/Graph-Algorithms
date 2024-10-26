# Graph Algorithms

This repository contains implementations of various graph algorithms, including Dijkstra's algorithm, Bellman-Ford algorithm, and Johnson's algorithm. The following sections provide the necessary theoretical background and explanations of the code.

## Table of Contents
- [1. Graphs](#1-graphs)
- [2. Dijkstra's Algorithm](#2-dijkstras-algorithm)
- [3. Bellman-Ford Algorithm](#3-bellman-ford-algorithm)
- [4. Johnson's Algorithm](#4-johnsons-algorithm)
- [5. Code Explanation](#5-code-explanation)

## 1. Graphs

**Definition**: A graph is a data structure consisting of vertices (or nodes) and edges (or connections) between them. Graphs can be directed (edges have direction) or undirected (edges do not have direction).

**Graph Representation**:
- **Adjacency Matrix**: A two-dimensional array where the element `graph[i][j]` represents the weight of the edge between vertices `i` and `j`.
- **Adjacency List**: An array of lists where each element of the array represents a vertex and contains a list of its neighbors.

**Types of Graphs**:
- **Weighted Graphs**: Edges have weights (cost).
- **Unweighted Graphs**: Edges do not have weights.
- **Cyclic Graphs**: Contain cycles (a path that returns to the starting vertex).
- **Acyclic Graphs**: Do not contain cycles.

## 2. Dijkstra's Algorithm

**Description**: Dijkstra's algorithm is used to find the shortest paths from a single source vertex to all other vertices in a graph with non-negative edge weights.

**Working Principle**:
1. **Initialization**: Set the distance to the source vertex to 0 and all other distances to infinity.
2. **Vertex Selection**: At each step, select the vertex with the minimum distance.
3. **Distance Update**: For each neighbor of the selected vertex, update the distances.
4. **Completion**: The process continues until all vertices have been processed.

**Complexity**: O((V + E) log V) when using a priority queue.

## 3. Bellman-Ford Algorithm

**Description**: The Bellman-Ford algorithm is used to find the shortest paths from a single source vertex to all other vertices in a graph that may contain negative edge weights.

**Working Principle**:
1. **Initialization**: Set the distance to the source vertex to 0 and all other distances to infinity.
2. **Edge Relaxation**: For each edge in the graph, update the distances to neighboring vertices for `V-1` iterations.
3. **Negative Cycle Check**: If a distance can be updated again, the graph contains a negative cycle.

**Complexity**: O(V * E).

## 4. Johnson's Algorithm

**Description**: Johnson's algorithm is used to find the shortest paths between all pairs of vertices in a graph that may contain negative edge weights but does not contain negative cycles.

**Working Principle**:
1. **Add a Dummy Vertex**: Create a new vertex and connect it to each existing vertex with zero weight edges.
2. **Run Bellman-Ford**: Use it to find the shortest distances from the dummy vertex.
3. **Transform Edge Weights**: Update the edge weights to eliminate negative weights.
4. **Run Dijkstra's Algorithm**: For each vertex, use Dijkstra's algorithm to find the shortest paths to all other vertices.

**Complexity**: O(V^2 * log V + V * E).

## 5. Code Explanation

### Reading the Graph

The function `readGraph` reads a graph from a binary file and constructs an adjacency matrix. It handles the conversion of weights, setting weights to infinity where appropriate.

### Bellman-Ford Algorithm

The `bellmanFord` function implements the Bellman-Ford algorithm. It initializes distances, relaxes edges, and checks for negative cycles.

### Dijkstra's Algorithm

The `dijkstra` function implements Dijkstra's algorithm using a priority queue to efficiently find the shortest paths from a source vertex.

### Johnson's Algorithm

The `johnson` function implements Johnson's algorithm. It adds a dummy vertex, runs the Bellman-Ford algorithm, transforms edge weights, and then applies Dijkstra's algorithm for each vertex.

### Graph Analysis

The `analyzeGraph` function computes the diameter, radius, and identifies central and peripheral vertices based on the distance matrix.

### Main Function

The `main` function orchestrates the reading of the graph, execution of algorithms, and output of results to a specified file.

---

