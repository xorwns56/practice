import java.util.*;
class Solution {
    public int[] solution(String today, String[] terms, String[] privacies) {
        HashMap<Character, Integer> map = new HashMap<>();
        for(int i = 0; i < terms.length; i++) map.put(terms[i].charAt(0), Integer.parseInt(terms[i].substring(2)));
        String[] sp = today.split("[.]");
        int today_y = Integer.parseInt(sp[0]);
        int today_m = Integer.parseInt(sp[1]);
        int today_d = Integer.parseInt(sp[2]);
        List<Integer> list = new ArrayList<>();
        for(int i = 0; i < privacies.length; i++){
            sp = privacies[i].split("[. ]");
            int d = Integer.parseInt(sp[2]) - 1;
            int m = Integer.parseInt(sp[1]) + map.get(sp[3].charAt(0));
            if(d == 0){
                m--;
                d = 28;
            }
            int y = Integer.parseInt(sp[0]) + m / 12;
            m = m % 12;
            if(m == 0){
                y--;
                m = 12;
            }
            if(today_y == y){
                if(today_m == m){
                    if(today_d <= d) continue;
                }else if(today_m < m) continue;
            }else if(today_y < y) continue;
            list.add(i + 1);
        }
        int[] answer = new int[list.size()];
        for(int i = 0; i < answer.length; i++) answer[i] = list.get(i);
        return answer;
    }
}