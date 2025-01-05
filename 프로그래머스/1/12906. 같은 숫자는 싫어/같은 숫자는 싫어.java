import java.util.*;

public class Solution {
    public int[] solution(int []arr) {
        Stack<Integer> stk = new Stack<Integer>();
        for(int i = 0; i < arr.length; i++){
            if(stk.isEmpty() || stk.peek() != arr[i]) stk.push(arr[i]);
        }
        int[] answer = new int[stk.size()];
        for(int i = 0; i < answer.length; i++) answer[answer.length - 1 - i] = stk.pop();
        return answer;
    }
}