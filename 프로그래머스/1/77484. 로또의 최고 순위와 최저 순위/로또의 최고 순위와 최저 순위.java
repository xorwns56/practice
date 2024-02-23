import java.util.*;
class Solution {
    public int[] solution(int[] lottos, int[] win_nums) {
        HashSet<Integer> set = new HashSet<>();
        int[] answer = new int[2];
        for(int i = 0; i < win_nums.length; i++) set.add(win_nums[i]);
        for(int i = 0; i < lottos.length; i++){
            if(lottos[i] == 0) answer[0]++;
            else if(set.contains(lottos[i])){
                answer[0]++;
                answer[1]++;
            }
        }
        answer[0] = answer[0] <= 1 ? 6 : 7 - answer[0];
        answer[1] = answer[1] <= 1 ? 6 : 7 - answer[1];
        return answer;
    }
}