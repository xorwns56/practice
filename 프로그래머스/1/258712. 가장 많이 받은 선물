import java.util.*;
class Solution {
    public int solution(String[] friends, String[] gifts) {
        HashMap<String, Integer> map = new HashMap();
        int[][] graph = new int[friends.length][friends.length];
        int[] count = new int[friends.length];
        for(int i = 0; i < friends.length; i++) map.put(friends[i], i);
        for(int i = 0; i < gifts.length; i++){
            String[] sp = gifts[i].split("\\s");
            int from = map.get(sp[0]);
            int to = map.get(sp[1]);
            graph[from][to]++;
            count[from]++;
            count[to]--;
        }
        int max_gift = 0;
        for(int i = 0; i < friends.length; i++){
            int gift = 0;
            for(int j = 0; j < friends.length; j++){
                if(i == j) continue;
                if(graph[i][j] - graph[j][i] > 0) gift++;
                else if(graph[i][j] - graph[j][i] == 0 && count[i] > count[j]) gift++;
            }
            if(max_gift < gift) max_gift = gift;
        }
        return max_gift;
    }
}
