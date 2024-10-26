#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <queue>
#include <cstring>

// Define infinity as the maximum value of int
#define INF std::numeric_limits<int>::max()

using namespace std;

// Read the graph from a binary file
vector<vector<int>> readGraph(const string& filename) {

    cout << "entered readGraph " << endl;

    ifstream infile(filename, ios::binary);  // Open the file for reading in binary mode
    if (!infile) {   // Check if the file is successfully opened
        cerr << "Error opening the file." << endl;
        exit(1);
    }

    int16_t size;  // Size of the graph
    infile.read(reinterpret_cast<char*>(&size), sizeof(size));  // Read the size of the graph from the file

    vector<vector<int>> graph(size, vector<int>(size));  // Create an adjacency matrix of size x size
    for (int i = 0; i < size; ++i) {   // Loop through all vertices of the graph
        for (int j = 0; j < size; ++j) {   // Loop through all edges for each vertex
            int16_t weight;  // Weight of the edge
            infile.read(reinterpret_cast<char*>(&weight), sizeof(weight));  // Read the weight of the edge from the file
            graph[i][j] = (weight == 0 && i != j) ? INF : weight;  // If weight is 0 and it's not a loop, set to infinity, otherwise set to weight
        }
    }
    infile.close();
    return graph;
}

// Bellman-Ford algorithm
bool bellmanFord(const vector<vector<int>>& graph, int src, vector<int>& dist) {

    cout << "entered Bellmanford " << endl;

    int V = graph.size();  // Get the number of vertices in the graph
    dist.assign(V, INF); // Initialize distances to all vertices as infinity
    dist[src] = 0; // Distance to the source vertex is 0

    for (int i = 0; i < V - 1; ++i) {
        for (int u = 0; u < V; ++u) {  // Loop through all vertices of the graph
            for (int v = 0; v < V; ++v) {  // Loop through all neighbors of vertex u
                if (graph[u][v] != INF && dist[u] != INF && dist[u] + graph[u][v] < dist[v]) {
                    // Update distance to vertex v if a shorter path is found
                    dist[v] = dist[u] + graph[u][v];
                }
            }
        }
    }

    // Check for negative cycles
    for (int u = 0; u < V; ++u) {  // Loop through all vertices
        for (int v = 0; v < V; ++v) {  // Loop through all neighbors
            if (graph[u][v] != INF && dist[u] != INF && dist[u] + graph[u][v] < dist[v]) {
                return true;  // Negative cycle detected
            }
        }
    }
    return false;  // No negative cycles
}

// Dijkstra's algorithm for a single source
vector<int> dijkstra(const vector<vector<int>>& graph, int src) {

    cout << "entered Dijkstra " << endl;

    int V = graph.size();   // Get the number of vertices in the graph
    vector<int> dist(V, INF);  // Initialize distances to all vertices as infinity
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;  // Create a priority queue to store pairs (distance, vertex)

    dist[src] = 0;  // Distance to the source vertex is 0
    pq.push({ 0, src });   // Add the source vertex to the queue

    while (!pq.empty()) {    // While the queue is not empty
        int u = pq.top().second;   // Extract the vertex with the minimum distance
        pq.pop();   // Remove this vertex from the queue

        for (int v = 0; v < V; ++v) {  // Loop through all neighbors
            if (graph[u][v] != INF && dist[u] != INF && dist[u] + graph[u][v] < dist[v]) {
                // Update distance to vertex v if a shorter path is found
                dist[v] = dist[u] + graph[u][v];
                pq.push({ dist[v], v });   // Add the updated distance to the queue
            }
        }
    }

    return dist;
}

// Johnson's algorithm
vector<vector<int>> johnson(const vector<vector<int>>& graph) {

    cout << "entered Johnson " << endl;

    int V = graph.size();
    vector<vector<int>> modifiedGraph(V, vector<int>(V));  // Create a modified adjacency matrix
    vector<vector<int>> dist(V, vector<int>(V, INF));  // Create a distance matrix and initialize it with infinities

    // Add a dummy vertex and run Bellman-Ford
    vector<int> h;
    if (bellmanFord(graph, 0, h)) {
        cerr << "Graph contains a negative cycle, cannot apply Johnson's algorithm." << endl;
        exit(1);
    }

    // Transform edge weights
    for (int u = 0; u < V; ++u) {  // Loop through all vertices
        for (int v = 0; v < V; ++v) {
            if (graph[u][v] != INF) {
                // Transform edge weights to eliminate negative weights
                modifiedGraph[u][v] = graph[u][v] + h[u] - h[v];
            }
            else {
                modifiedGraph[u][v] = INF;  // If there is no edge, set to infinity
            }
        }
    }

    // Dijkstra's algorithm for each vertex
    for (int src = 0; src < V; ++src) {
        vector<int> d = dijkstra(modifiedGraph, src);
        for (int v = 0; v < V; ++v) {
            if (d[v] != INF) {
                // Restore original weights and save results in the distance matrix
                dist[src][v] = d[v] + h[v] - h[src];
            }
        }
    }

    return dist;
}

// Analyze the graph to compute diameter, radius, central and peripheral vertices
void analyzeGraph(const vector<vector<int>>& dist, int& diameter, int& radius, vector<int>& centralVertices, vector<int>& peripheralVertices) {

    cout << "entered analyzeGraph " << endl;

    int V = dist.size();
    vector<int> eccentricity(V, 0);  // Vector to store the eccentricities of each vertex

    diameter = 0;
    radius = INF;

    for (int u = 0; u < V; ++u) {   // Loop through all vertices
        for (int v = 0; v < V; ++v) {
            if (dist[u][v] != INF) {
                // Update the eccentricity of vertex u
                eccentricity[u] = max(eccentricity[u], dist[u][v]);
            }
        }
        diameter = max(diameter, eccentricity[u]);  // Update the diameter of the graph with the largest eccentricity
        radius = min(radius, eccentricity[u]); // Update the radius of the graph with the smallest eccentricity
    }

    for (int u = 0; u < V; ++u) {
        if (eccentricity[u] == radius) {
            // If eccentricity equals the radius, add the vertex to central vertices
            centralVertices.push_back(u);
        }
        if (eccentricity[u] == diameter) {
            // If eccentricity equals the diameter, add the vertex to peripheral vertices
            peripheralVertices.push_back(u);
        }
    }
}

// Main function of the program
int main(int argc, char* argv[]) {

    cout << "entered main " << endl;

    // Check if enough command line arguments are provided
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " inputfile [-o outputfile]" << endl;
        return 1;
    }

    string inputFile = argv[1];  // Read the input file name from command line arguments
    string outputFile = "output.txt"; // Set the default output file name

    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            outputFile = argv[++i];
        }
    }

    // Read the graph
    vector<vector<int>> graph = readGraph(inputFile);

    // Bellman-Ford algorithm
    vector<int> dist;
    bool hasNegativeCycle = bellmanFord(graph, 0, dist);

    ofstream outfile(outputFile); // Open the output file to write results

    if (hasNegativeCycle) {
        outfile << "Graph contains edges with negative weight." << endl;
    } else {
        outfile << "Graph does not contain edges with negative weight." << endl;
        outfile << "Shortest paths lengths:" << endl;

        for (int v = 0; v < dist.size(); ++v) {
            outfile << "0 - " << v << ": ";
            if (dist[v] == INF) {
                outfile << "∞" << endl;  // Write "∞" for infinity
            } else {
                outfile << dist[v] << endl; // Write the distance value
            }
        }

        // Johnson's algorithm
        vector<vector<int>> allPairsDist = johnson(graph);

        // Analyze the distance matrix to get diameter, radius, and central/peripheral vertices
        int diameter, radius;
        vector<int> centralVertices, peripheralVertices;
        analyzeGraph(allPairsDist, diameter, radius, centralVertices, peripheralVertices);

        // Output results
        outfile << "Diameter of the graph: " << diameter << endl;
        outfile << "Radius of the graph: " << radius << endl;

        outfile << "Central vertices: ";
        for (int v : centralVertices) {
            outfile << v << " ";
        }
        outfile << endl;

        outfile << "Peripheral vertices: ";
        for (int v : peripheralVertices) {
            outfile << v << " ";
        }
        outfile << endl;
    }

    outfile.close();
    return 0;
}
