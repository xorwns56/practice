using System;
using System.Collections.Generic;
public class Solution {
    int[] visited;
    List<List<int>> graph;
    public int solution(int n, int[,] edge) {
        visited = new int[n];
        graph = new List<List<int>>();
        for(int i = 0; i< n; i++){
            visited[i] = 0;
            graph.Add(new List<int>());
        }
        for(int i = 0; i<edge.GetLength(0); i++){
            graph[edge[i,0]-1].Add(edge[i,1]-1);
            graph[edge[i,1]-1].Add(edge[i,0]-1);
        }
        return BFS(0);
    }
    
    int BFS(int first){
        Queue<int> q = new Queue<int>();
        int max = 1;
        int count = 1;
        visited[first] = 1;
        q.Enqueue(first);
        while(q.Count > 0){
            int v = q.Dequeue();
            for(int i = 0; i< graph[v].Count; i++){
                if(visited[graph[v][i]] == 0){
                    if(max < visited[v] + 1){
                        max = visited[v] + 1;
                        count = 1;
                    }
                    else if(max == visited[v] + 1) count++;
                    visited[graph[v][i]] = visited[v] + 1;
                    q.Enqueue(graph[v][i]);
                }
            }
        }
        return count;
    }
}