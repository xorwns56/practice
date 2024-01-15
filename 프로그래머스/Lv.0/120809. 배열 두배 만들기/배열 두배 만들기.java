import java.util.*;
class Solution {
    public int[] solution(int[] numbers) {
        int[] answer = Arrays.copyOf(numbers, numbers.length);
        for(int i = 0; i < answer.length; i++) answer[i] *= 2;
        return answer;
    }
}