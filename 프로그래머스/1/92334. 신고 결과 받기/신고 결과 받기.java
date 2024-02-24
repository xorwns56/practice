import java.util.*;
class Solution {
    public int[] solution(String[] id_list, String[] report, int k) {
        HashMap<String, Integer> map = new HashMap<>();
        boolean[][] graph = new boolean[id_list.length][id_list.length];
        int[] reported = new int[id_list.length];
        for(int i = 0; i < id_list.length; i++) map.put(id_list[i], i);
        for(int i = 0; i < report.length; i++){
            String[] sp = report[i].split("\\s");
            int from = map.get(sp[0]);
            int to = map.get(sp[1]);
            if(!graph[from][to]) reported[to]++;
            graph[from][to] = true;
        }
        int[] answer = new int[id_list.length];
        for(int i = 0; i < id_list.length; i++){
            if(reported[i] >= k){
                for(int j = 0; j < id_list.length; j++){
                    if(graph[j][i]) answer[j]++;
                }
            }
        }
        return answer;
    }
}