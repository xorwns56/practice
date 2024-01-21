import java.util.*;
import java.util.stream.IntStream;
class Solution {
    public int solution(int i, int j, int k) {
        int answer = 0;
        for(int n = i; n <= j; n++){
            int curr = n;
            while(curr > 0){
                if(curr % 10 == k) answer++;
                curr /= 10;
            }
        }
        return answer;
    }
}