import java.util.*;
class Solution {
    public String[] solution(String[] players, String[] callings) {
        HashMap<String, Integer> map = new HashMap<>();
        for(int i = 0; i < players.length; i++) map.put(players[i], i);
        for(int i = 0; i < callings.length; i++){
            int curr = map.get(callings[i]);
            players[curr] = players[curr - 1];
            players[curr - 1] = callings[i];
            map.put(players[curr], curr);
            map.put(players[curr - 1], curr - 1);
        }
        String[] answer = new String[players.length];
        for(Map.Entry<String, Integer> entry : map.entrySet()) answer[entry.getValue()] = entry.getKey();
        return answer;
    }
}