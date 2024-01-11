import java.util.*;
class Solution {
    public int[] solution(int[] arr) {
        Stack<Integer> stk = new Stack<>();
        for(int i = 0; i < arr.length; i++){
            if(stk.isEmpty() || stk.peek() != arr[i]) stk.push(arr[i]);
            else if(stk.peek() == arr[i]) stk.pop();
        }
        if(stk.isEmpty()) stk.push(-1);
        
        int[] answer = new int[stk.size()];
        for(int i = answer.length - 1; i >= 0; i--) answer[i] = stk.pop();
        return answer;
    }
}