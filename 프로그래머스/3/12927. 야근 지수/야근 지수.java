import java.util.*;
class Solution {
    public long solution(int n, int[] works) {
        long answer = 0;
        Arrays.sort(works);
        int idx = works.length - 1;
        while(n > 0){
            if(works[idx] == 0) break;
            works[idx]--;
            if(idx == 0 || works[idx] >= works[idx - 1]) idx = works.length - 1;
            else if(works[idx] < works[idx - 1]) idx--;
            n--;
        }
        for(int i = 0; i < works.length; i++) answer += works[i] * works[i];
        return answer;
    }
}