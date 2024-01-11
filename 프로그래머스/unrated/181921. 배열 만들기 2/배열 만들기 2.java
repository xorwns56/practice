import java.util.*;
class Solution {
    public int[] solution(int l, int r) {
        List<Integer> list = new ArrayList<>();
        for(int i = 1;;i++){
            int tmp = i;
            int digit = 0;
            int rslt = 0;
            while(tmp != 0){
                rslt += (tmp & 1) * Math.pow(10, digit) * 5;
                tmp >>= 1;
                digit++;
            }
            if(r < rslt) break;
            else if(l <= rslt) list.add(rslt);
        }
        if(list.size() == 0) list.add(-1);
        int[] answer = new int[list.size()];
        for(int i = 0; i < list.size(); i++) answer[i] = list.get(i);
        return answer;
    }
}